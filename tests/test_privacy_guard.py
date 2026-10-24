"""The privacy guard is the last line of defence — verify it actually catches things."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

from privacy_guard import (  # noqa: E402
    _git,
    expand_paths,
    in_git_work_tree,
    main,
    scan_path,
    tracked_and_untracked_paths,
)


def _rules(path: Path) -> set[str]:
    return {violation.rule for violation in scan_path(path)}


def test_catches_vcf_content(tmp_path: Path) -> None:
    # Renamed to .txt on purpose: name-based blocking alone is not enough.
    path = tmp_path / "notes.txt"
    path.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
        "15\t40163000\t.\tC\tT\t50\tPASS\t.\tGT:AD:DP:GQ:PGT:PID\t0/1:10,10:20:99:0|1:40163000\n"
    )
    rules = _rules(path)
    assert "vcf-header" in rules
    assert "vcf-column-header" in rules
    assert "vcf-genotype-field" in rules


def test_catches_fastq_content(tmp_path: Path) -> None:
    path = tmp_path / "reads.txt"
    path.write_text("@READ1\n" + "ACGT" * 20 + "\n+\n" + "I" * 80 + "\n")
    assert _rules(path) & {"fastq-record", "long-read-sequence"}


def test_catches_forbidden_extension(tmp_path: Path) -> None:
    path = tmp_path / "sample.vcf.gz"
    path.write_text("harmless looking text")
    assert "forbidden-extension" in _rules(path)


def test_catches_variant_table(tmp_path: Path) -> None:
    path = tmp_path / "candidates.csv"
    rows = "\n".join(f"chr15,{40000000 + i},A,G" for i in range(60))
    path.write_text("chrom,pos,ref,alt\n" + rows)
    assert "variant-table" in _rules(path)


@pytest.mark.parametrize(
    ("secret", "rule"),
    [
        ('token = "hf_' + "a" * 34 + '"', "hf-token"),
        ('key = "AKIA' + "A" * 16 + '"', "aws-key"),
    ],
)
def test_catches_secrets(tmp_path: Path, secret: str, rule: str) -> None:
    path = tmp_path / "creds.py"
    path.write_text(secret + "\n")
    assert rule in _rules(path)


def test_allows_ordinary_documentation(tmp_path: Path) -> None:
    path = tmp_path / "notes.md"
    path.write_text(
        "# Notes\n\nBUB1B is on chromosome 15q15.1 and encodes BubR1.\n"
        "A handful of coordinates in prose is fine: chr15:40163000.\n"
    )
    assert not scan_path(path)


def test_allow_marker_waives_content_rules(tmp_path: Path) -> None:
    path = tmp_path / "format-example.md"
    path.write_text("privacy-guard: allow\n\n##fileformat=VCFv4.2\n")
    assert not scan_path(path)


def test_repository_tree_is_clean() -> None:
    """The repo itself must always pass. This test is the regression guard."""
    sys.path.insert(0, str(ROOT / "scripts"))
    from privacy_guard import scan_path as scan
    from privacy_guard import tracked_and_untracked_paths

    violations = [v for p in tracked_and_untracked_paths() for v in scan(p)]
    assert not violations, "privacy violations in the working tree:\n" + "\n".join(str(v) for v in violations)


# ---------------------------------------------------------------------------
# 🔴 Directory arguments
#
# ``scan_path`` returns early for anything that is not a file, so before
# ``expand_paths`` existed a directory argument scanned NOTHING and the script
# reported "OK - 1 file(s) scanned". An earlier caller passed exactly such a
# directory, so the submission's privacy check was a green light that never ran. These tests exist so that cannot come back.
# ---------------------------------------------------------------------------


def test_scan_path_alone_is_blind_to_directories(tmp_path: Path) -> None:
    """Document the sharp edge that made the bug possible."""
    (tmp_path / "leak.txt").write_text("##fileformat=VCFv4.2\n")
    assert scan_path(tmp_path) == []  # the directory itself: nothing
    assert _rules(tmp_path / "leak.txt")  # the file inside it: caught


def test_expand_paths_walks_directories_recursively(tmp_path: Path) -> None:
    nested = tmp_path / "evidence" / "deep"
    nested.mkdir(parents=True)
    (tmp_path / "top.md").write_text("fine\n")
    (nested / "buried.txt").write_text("also fine\n")

    expanded = expand_paths([tmp_path])

    assert set(expanded) == {tmp_path / "top.md", nested / "buried.txt"}


def test_expand_paths_finds_a_leak_buried_in_a_subdirectory(tmp_path: Path) -> None:
    """The regression: a leak one level down must not be invisible."""
    nested = tmp_path / "evidence"
    nested.mkdir()
    (nested / "dump.txt").write_text(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
    )

    violations = [v for path in expand_paths([tmp_path]) for v in scan_path(path)]

    assert {"vcf-header", "vcf-column-header"} <= {v.rule for v in violations}


def test_expand_paths_skips_vcs_and_cache_directories(tmp_path: Path) -> None:
    (tmp_path / ".git").mkdir()
    (tmp_path / ".git" / "objectish").write_text("##fileformat=VCFv4.2\n")
    (tmp_path / "__pycache__").mkdir()
    (tmp_path / "__pycache__" / "x.pyc").write_text("noise\n")
    (tmp_path / "real.md").write_text("content\n")

    assert expand_paths([tmp_path]) == [tmp_path / "real.md"]


def test_expand_paths_passes_files_and_missing_paths_through(tmp_path: Path) -> None:
    real = tmp_path / "a.md"
    real.write_text("x\n")
    missing = tmp_path / "nope.md"

    assert expand_paths([real, missing]) == [real, missing]


def test_expand_paths_deduplicates(tmp_path: Path) -> None:
    f = tmp_path / "a.md"
    f.write_text("x\n")

    assert expand_paths([tmp_path, f]) == [f]


# Where the submission artefacts live.
SUBMISSION_DIRS = ("submission",)


def test_submission_directory_is_actually_scanned() -> None:
    """Scripts invoke the guard as `--paths <dir>`; prove that reaches files.

    The regression: `scan_path()` returned early for anything that was not a
    file, so `--paths submission` reported `OK - 1 file(s) scanned`
    having read nothing at all.
    """
    present = [ROOT / name for name in SUBMISSION_DIRS if (ROOT / name).is_dir()]
    if not present:
        pytest.skip("no submission directory in this checkout")

    for directory in present:
        expanded = expand_paths([directory])
        assert len(expanded) > 1, f"{directory} expanded to nothing"
        assert all(path.is_file() for path in expanded)
        assert not [v for path in expanded for v in scan_path(path)]


# ---------------------------------------------------------------------------
# 🔴 A scan that reads nothing must never report success.
#
# Two ways that happened. `scan_path()` returned early for directories, so
# `--paths <dir>` scanned nothing (fixed; covered above). And `_git()` swallowed
# git's exit code into an empty list, so `--tree` outside a git checkout (an
# unpacked archive, say) printed `OK - 0 file(s) scanned`.
# ---------------------------------------------------------------------------


def test_git_failure_is_loud_not_empty() -> None:
    """A git error must raise, not read as 'no files'."""
    with pytest.raises(RuntimeError):
        _git("ls-files", "--not-a-real-git-flag")


def test_tree_mode_walks_the_filesystem_outside_a_git_checkout(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A copy without git metadata (an unpacked archive) must still be scanned."""
    (tmp_path / "a.md").write_text("content\n")
    (tmp_path / "pkg").mkdir()
    (tmp_path / "pkg" / "b.py").write_text("x = 1\n")

    import privacy_guard

    monkeypatch.setattr(privacy_guard, "REPO_ROOT", tmp_path)
    assert not in_git_work_tree()

    found = tracked_and_untracked_paths()
    assert {p.name for p in found} == {"a.md", "b.py"}


def test_tree_mode_refuses_to_pass_when_it_scanned_nothing(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    import privacy_guard

    monkeypatch.setattr(privacy_guard, "REPO_ROOT", tmp_path)
    monkeypatch.setattr(privacy_guard, "tracked_and_untracked_paths", list)

    assert main(["--tree"]) == 2
    assert "did not happen" in capsys.readouterr().err


def test_paths_mode_refuses_to_pass_when_it_scanned_nothing(
    capsys: pytest.CaptureFixture[str],
) -> None:
    import privacy_guard

    original = privacy_guard.expand_paths
    try:
        privacy_guard.expand_paths = lambda paths: []
        assert main(["--paths", "whatever"]) == 2
    finally:
        privacy_guard.expand_paths = original
    assert "did not happen" in capsys.readouterr().err
