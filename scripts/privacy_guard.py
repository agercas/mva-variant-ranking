#!/usr/bin/env python3
"""Block individual-level genomic data and secrets from entering git history.

This is the enforcement arm of the project's privacy rules. It is deliberately
noisy: a false positive costs you one ``--allow`` annotation, a false negative
costs a child their privacy.

Usage
-----
    python scripts/privacy_guard.py            # scan git-staged files
    python scripts/privacy_guard.py --tree     # scan the whole working tree
    python scripts/privacy_guard.py --paths a b

``--paths`` accepts directories and walks them recursively. 🔴 It did not always:
``scan_path`` returns early for anything that is not a file, so a directory
argument used to scan **nothing** and report ``OK - 1 file(s) scanned``. That is
the worst possible failure mode for this script - a green check that means the
check did not run - and an earlier script passed exactly such a directory.
See ``expand_paths``.

Exit codes: 0 = clean, 1 = violation found, 2 = usage error.
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent

# ---------------------------------------------------------------------------
# Rules
# ---------------------------------------------------------------------------

#: Extensions that may never be committed, whatever their content.
FORBIDDEN_SUFFIXES: frozenset[str] = frozenset(
    {
        ".vcf",
        ".gvcf",
        ".bcf",
        ".bam",
        ".bai",
        ".cram",
        ".crai",
        ".sam",
        ".fastq",
        ".fq",
        ".fa",
        ".fasta",
        ".fai",
        ".2bit",
        ".gtc",
        ".idat",
        ".cel",
        ".ped",
        ".fam",
        ".bim",
        ".phenopacket",
    }
)

#: Directories whose contents may never be committed.
FORBIDDEN_DIRS: tuple[str, ...] = ("data/private", "data/cache", "secrets")

#: Content signatures. (name, compiled regex, human explanation)
CONTENT_SIGNATURES: list[tuple[str, re.Pattern[str], str]] = [
    (
        "vcf-header",
        re.compile(r"^##fileformat=VCFv", re.MULTILINE),
        "VCF file header — individual-level variant calls",
    ),
    (
        "vcf-column-header",
        re.compile(r"^#CHROM\s+POS\s+ID\s+REF\s+ALT", re.MULTILINE),
        "VCF column header — individual-level variant calls",
    ),
    (
        "sam-header",
        re.compile(r"^@HD\s+VN:\d", re.MULTILINE),
        "SAM/BAM header — aligned sequencing reads",
    ),
    (
        "fastq-record",
        re.compile(r"^@[!-~]+\n[ACGTNacgtn]{40,}\n\+", re.MULTILINE),
        "FASTQ record — raw sequencing reads",
    ),
    (
        "long-read-sequence",
        re.compile(r"[ACGTN]{200,}"),
        "long nucleotide run — possible embedded sequence data",
    ),
    (
        "vcf-genotype-field",
        re.compile(r"\bGT:AD:DP:GQ(:PL|:PGT|:PID)?\b"),
        "VCF FORMAT/genotype field — individual-level genotypes",
    ),
    (
        "hf-token",
        re.compile(r"\bhf_[A-Za-z0-9]{30,}\b"),
        "Hugging Face access token",
    ),
    (
        "aws-key",
        re.compile(r"\bAKIA[0-9A-Z]{16}\b"),
        "AWS access key id",
    ),
    (
        "anthropic-key",
        re.compile(r"\bsk-ant-[A-Za-z0-9_\-]{20,}"),
        "Anthropic API key",
    ),
    (
        "generic-bearer",
        re.compile(r"(?i)\b(api[_-]?key|secret|password|token)\s*[=:]\s*['\"][^'\"\s]{16,}['\"]"),
        "hard-coded credential",
    ),
]

#: A table of many variant rows is individual-level data even as CSV/TSV.
VARIANT_ROW = re.compile(
    r"^(chr)?(\d{1,2}|[XYM]|MT)[\t,;]\d{3,}[\t,;][ACGTN]+[\t,;][ACGTN]+",
    re.MULTILINE | re.IGNORECASE,
)
#: Above this many variant-shaped rows in one file we call it a variant table.
VARIANT_ROW_LIMIT = 25

#: Escape hatch for reviewed files (e.g. documentation showing a format example).
ALLOW_MARKER = "privacy-guard: allow"

#: Files we never scan for content (binary-ish or self-referential).
SKIP_NAMES = {"privacy_guard.py"}
SKIP_SUFFIXES = {".png", ".jpg", ".jpeg", ".gif", ".pdf", ".mp4", ".webm", ".ico", ".woff", ".woff2"}

MAX_TRACKED_BYTES = 5 * 1024 * 1024  # 5 MB: nothing legitimate here is bigger


@dataclass(frozen=True)
class Violation:
    path: str
    rule: str
    detail: str

    def __str__(self) -> str:  # pragma: no cover - formatting only
        return f"  {self.path}\n      [{self.rule}] {self.detail}"


# ---------------------------------------------------------------------------
# Scanning
# ---------------------------------------------------------------------------


def _git(*args: str) -> list[str]:
    """Run a git command and return its non-empty output lines.

    🔴 Raises on failure rather than returning an empty list. Swallowing the
    error makes a git that could not run indistinguishable from a repository
    with no files -- and this script reports the second case as success.

    Raises:
        RuntimeError: git exited non-zero, or is not installed.
    """
    try:
        out = subprocess.run(["git", *args], cwd=REPO_ROOT, capture_output=True, text=True, check=False)
    except OSError as exc:  # git not installed
        raise RuntimeError(f"cannot run git: {exc}") from exc
    if out.returncode != 0:
        raise RuntimeError(f"git {' '.join(args)} failed: {out.stderr.strip() or out.returncode}")
    return [line for line in out.stdout.splitlines() if line.strip()]


def in_git_work_tree() -> bool:
    """Whether REPO_ROOT is inside a git working tree."""
    try:
        return _git("rev-parse", "--is-inside-work-tree") == ["true"]
    except RuntimeError:
        return False


def staged_paths() -> list[Path]:
    return [REPO_ROOT / p for p in _git("diff", "--cached", "--name-only", "--diff-filter=ACMR")]


def tracked_and_untracked_paths() -> list[Path]:
    """Every file in the working tree, git-aware when possible.

    🔴 Falls back to walking the filesystem when this is not a git checkout.
    It previously relied on ``git ls-files`` alone, whose failure was swallowed
    into an empty list -- so running ``--tree`` outside a git checkout (an
    unpacked archive, say) printed ``OK - 0 file(s) scanned``. A guard that reports success having read nothing
    is worse than no guard, because it is believed.
    """
    if not in_git_work_tree():
        return expand_paths([REPO_ROOT])
    seen = dict.fromkeys(_git("ls-files") + _git("ls-files", "--others", "--exclude-standard"))
    return [REPO_ROOT / p for p in seen]


#: Directories never worth walking into when a caller names a parent directory.
WALK_SKIP_DIRS = frozenset({".git", ".venv", "__pycache__", ".pytest_cache", ".ruff_cache", ".mypy_cache"})


def expand_paths(paths: list[Path]) -> list[Path]:
    """Expand directories to the files inside them, recursively.

    ``scan_path`` returns early for anything that is not a file. Without this
    expansion a directory argument scans nothing and the script reports success
    -- a green check that means the check never ran. An earlier caller passed
    the submission directory this way, so it was never actually scanned.

    Args:
        paths: Files and/or directories as given on the command line.

    Returns:
        Every file reachable from ``paths``, de-duplicated, in a stable order.
        Non-existent paths are passed through unchanged so that ``scan_path``
        -- not this helper -- decides what to do with them.
    """
    out: list[Path] = []
    for path in paths:
        if path.is_dir():
            for child in sorted(path.rglob("*")):
                if not child.is_file():
                    continue
                if WALK_SKIP_DIRS & set(child.relative_to(path).parts):
                    continue
                out.append(child)
        else:
            out.append(path)
    return list(dict.fromkeys(out))


def _rel(path: Path) -> str:
    try:
        return str(path.relative_to(REPO_ROOT))
    except ValueError:
        return str(path)


def scan_path(path: Path) -> list[Violation]:
    """Return every violation found in ``path``."""
    rel = _rel(path)
    violations: list[Violation] = []

    if not path.is_file():
        return violations

    # A .gz/.bgz wrapper hides the real suffix, so look one level deeper.
    suffixes = [s.lower() for s in path.suffixes]
    effective = {s for s in suffixes if s not in {".gz", ".bgz", ".zip", ".bz2", ".xz"}}
    if effective & FORBIDDEN_SUFFIXES:
        violations.append(
            Violation(
                rel,
                "forbidden-extension",
                f"{sorted(effective & FORBIDDEN_SUFFIXES)} is a genomic data format",
            )
        )

    if any(rel == d or rel.startswith(f"{d}/") for d in FORBIDDEN_DIRS):
        violations.append(Violation(rel, "quarantined-directory", "this directory must never be committed"))

    if path.name in SKIP_NAMES or path.suffix.lower() in SKIP_SUFFIXES:
        return violations

    size = path.stat().st_size
    if size > MAX_TRACKED_BYTES:
        violations.append(
            Violation(
                rel,
                "oversized-file",
                f"{size / 1e6:.1f} MB exceeds the {MAX_TRACKED_BYTES / 1e6:.0f} MB limit",
            )
        )

    try:
        text = path.read_text(encoding="utf-8", errors="strict")
    except (UnicodeDecodeError, OSError):
        # Binary content is suspicious anywhere in this repo.
        if size > 64 * 1024:
            violations.append(Violation(rel, "binary-blob", "unreadable binary file larger than 64 KB"))
        return violations

    if ALLOW_MARKER in text:
        return [v for v in violations if v.rule == "quarantined-directory"]

    for name, pattern, explanation in CONTENT_SIGNATURES:
        if pattern.search(text):
            violations.append(Violation(rel, name, explanation))

    n_variant_rows = len(VARIANT_ROW.findall(text))
    if n_variant_rows > VARIANT_ROW_LIMIT:
        violations.append(
            Violation(
                rel,
                "variant-table",
                f"{n_variant_rows} variant-shaped rows — looks like an individual-level variant table",
            )
        )

    return violations


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--tree", action="store_true", help="scan the whole working tree")
    group.add_argument("--paths", nargs="+", help="scan specific paths")
    args = parser.parse_args(argv)

    if args.paths:
        paths = expand_paths([Path(p).resolve() for p in args.paths])
        scope = "given paths"
    elif args.tree:
        paths = tracked_and_untracked_paths()
        scope = "working tree"
    else:
        paths = staged_paths()
        scope = "staged changes"

    # 🔴 Scanning nothing is a failure, never a pass. Staged mode is the one
    # one exception: a pre-commit hook with nothing staged has nothing to do.
    if not paths and (args.tree or args.paths):
        print(
            f"privacy-guard: ERROR — nothing to scan in {scope}. "
            "Refusing to report success for a scan that did not happen.",
            file=sys.stderr,
        )
        return 2

    violations = [v for p in paths for v in scan_path(p)]

    if not violations:
        print(f"privacy-guard: OK — {len(paths)} file(s) scanned in {scope}.")
        return 0

    print(f"privacy-guard: BLOCKED — {len(violations)} violation(s) in {scope}.\n", file=sys.stderr)
    for violation in violations:
        print(violation, file=sys.stderr)
    print(
        "\nThis repository must never contain individual-level genomic data or secrets.\n"
        "  • Genomic data belongs in data/private/ (git-ignored).\n"
        "  • Secrets belong in .env (git-ignored).\n"
        f"  • A reviewed false positive can be annotated with '{ALLOW_MARKER}' in the file.\n"
        "See the Privacy section of the README.",
        file=sys.stderr,
    )
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
