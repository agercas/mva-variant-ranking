#!/usr/bin/env python3
"""Turn a raw result file into the version that is safe to publish.

Several checks print more about the child than the report needs: the genes where
the child has a homozygous change, rsIDs and genotypes of common benign variants,
every mitochondrial variant (which together give the maternal lineage), and the
child's phenotype terms. Only what an argument turns on should be published: the
two variants of the diagnosis, the secondary finding, and counts.

This script reads a raw result file on stdin and writes the public view on stdout.
It fails (exit 2) if the file does not look like the kind it was told to expect,
so a changed output format cannot silently pass unredacted text through.

Usage::

    python3 scripts/biallelic_ranking.py ... | python3 scripts/public_view.py --kind ranking

Kinds: ``ranking``, ``differential``, ``expanded``, ``mito``, ``phenotype``, ``completeness``.
"""

from __future__ import annotations

import argparse
import re
import sys

#: Genes that may be named in the public view.
REPORTED_GENES = {"BUB1B", "LZTR1"}

#: Competing-diagnosis genes the report names. Any other gene in the differential list is
#: counted but not named, because the list as a whole would show the child's clinical picture.
NAMED_DIFFERENTIALS = {
    "HRAS",
    "CDKN1C",
    "PTPN11",
    "SOS1",
    "RAF1",
    "KRAS",
    "NRAS",
    "BRAF",
    "MAP2K1",
    "SHOC2",
    "CBL",
    "RIT1",
    "CEP57",
    "TRIP13",
    "BUB1",
}


class Unexpected(Exception):
    """The input did not look like the kind of file the view was written for."""


def ranking(text: str) -> str:
    """Name only BUB1B in the ranking tables; other genes become '(withheld)'."""
    out, in_table, seen_table = [], False, False
    for line in text.splitlines():
        if line.lstrip().startswith("rank  gene"):
            in_table, seen_table = True, True
            out.append(line)
            continue
        if in_table:
            m = re.match(r"^(\s*\d+\s+)(\S+)(\s+)(.*)$", line)
            if m and not line.startswith("-"):
                gene = m.group(2)
                if gene not in REPORTED_GENES:
                    gene = "(withheld)"
                out.append(f"{m.group(1)}{gene:<15}{m.group(4)}")
                continue
            if not line.strip():
                in_table = False
        out.append(line)
    if not seen_table:
        raise Unexpected("no ranking table found")
    return "\n".join(out) + "\n"


def differential(text: str) -> str:
    """Keep gene-level conclusions; drop rsIDs, genotypes and phenotype hints."""
    lines = text.splitlines()
    out: list[str] = []
    i = 0
    seen_part2 = False
    block_head = re.compile(r"^([A-Z0-9]+) [—-] (.*)$")
    hidden = {"none": 0, "benign": 0, "unclassified": 0, "pathogenic": 0}
    while i < len(lines):
        line = lines[i]
        # Part 1: drop the per-variant carrier lines, keep BUB1B and LZTR1 and the count.
        if line.startswith("  ") and re.match(r"^\s+[A-Z0-9]+\s+\d+:\d+ ", line):
            if line.split()[0] in REPORTED_GENES:
                out.append(re.sub(r"(DP=\d+).*$", r"\1", line).rstrip())
            i += 1
            continue
        m = block_head.match(line)
        if m and "PART" not in line:
            seen_part2 = True
            gene, desc = m.group(1), m.group(2)
            desc = re.sub(r"\s*-\s*blind rank \d+.*$", "", desc)
            desc = re.sub(r"\s*\([^)]*(nephrocalcinosis|growth)[^)]*\)", "", desc)
            block = []
            i += 1
            while i < len(lines) and lines[i].strip():
                block.append(lines[i])
                i += 1
            if gene in REPORTED_GENES:
                out.append(f"{gene} — {desc}")
                out.extend(re.sub(r"\s+GT=\S+ DP=\d+", "", b) for b in block)
            else:
                variants = [b for b in block if re.search(r"\sGT=\S+", b)]
                clinvar = [(re.search(r"ClinVar=(\S+)", b) or [None, "-"])[1] for b in variants]
                benign = sum("enign" in c for c in clinvar)
                pathogenic = sum("athogenic" in c and "enign" not in c for c in clinvar)
                unclassified = len(clinvar) - benign - pathogenic
                excluded = any("EXCLUDED" in b for b in block)
                if not variants and not excluded:
                    raise Unexpected(f"unrecognised block for {gene}")
                if not variants:
                    verdict = "no damaging variants"
                else:
                    verdict = (
                        f"{len(variants)} damaging-annotated variant(s): {benign} ClinVar benign or likely "
                        f"benign, {unclassified} not classified, {pathogenic} pathogenic or likely pathogenic"
                    )
                if gene in NAMED_DIFFERENTIALS:
                    out.append(f"{gene} — {desc}")
                    out.append(f"   {verdict}" + ("  -> EXCLUDED" if excluded else ""))
                else:
                    if not variants:
                        hidden["none"] += 1
                    elif pathogenic:
                        hidden["pathogenic"] += 1
                    elif unclassified:
                        hidden["unclassified"] += 1
                    else:
                        hidden["benign"] += 1
                    if "HIDDEN_SUMMARY" not in out:
                        out.append("HIDDEN_SUMMARY")
                    i += 1
                    continue
            out.append("")
            i += 1  # the blank line that ended the block
            continue
        out.append(line)
        i += 1
    if not seen_part2:
        raise Unexpected("no differential blocks found")
    total = sum(hidden.values())
    summary = (
        f"{total} further named competing diagnoses (genes not named): {hidden['none']} have no damaging "
        f"variant; {hidden['benign']} carry only variants that ClinVar calls benign or likely benign; "
        f"{hidden['unclassified']} carry benign variants plus one that is not classified; "
        f"{hidden['pathogenic']} carry a pathogenic or likely pathogenic variant.\n"
    )
    result = "\n".join(out).rstrip() + "\n"
    return result.replace("HIDDEN_SUMMARY\n", summary + "\n") if total else result


def mito(text: str) -> str:
    """Keep coverage, the noise floor and the verdict; drop every per-position variant list."""
    if "EMPIRICAL NOISE FLOOR" not in text:
        raise Unexpected("not a mitochondrial screen")
    position_row = re.compile(r"^\s+\d+ [ACGT]>[ACGT]\s+[\d.]+%\s+depth\s+\d+")
    out: list[str] = []
    skip_next = False
    for line in text.splitlines():
        if skip_next:
            skip_next = False
            if "continuum and are reported" in line:
                continue
        if position_row.match(line) or re.match(r"^\s+\d+ [ACGT]>[ACGT] at [\d.]+% ", line):
            continue
        if "sit INSIDE the noise" in line:
            out.append("  For completeness, some P/LP positions sit INSIDE the noise continuum and are")
            out.append("  reported as artefacts, not findings (positions withheld).")
            skip_next = True
            continue
        out.append(line)
    result = "\n".join(out).rstrip() + "\n"
    if re.search(r"^\s+\d+ [ACGT]>[ACGT]\s+[\d.]+%\s+depth\s+\d+\s*<-", result, re.M):
        raise Unexpected("a per-position variant row is still present")
    return result


#: Diseases the report names in the phenotype ranking. The full top of the list would show
#: which clinical features the child has, so every other disease is dropped from the tables.
NAMED_DISEASES = ("Mosaic variegated", "MVA", "Costello", "Beckwith", "Li-Fraumeni")


def keep_named_disease_rows(text: str) -> str:
    """Drop ranking rows for diseases the report does not name; say how many were dropped."""
    row = re.compile(r"^\s+\d+\s+[\d.]+\s+(\S.*)$")
    out: list[str] = []
    dropped = 0
    for line in text.splitlines():
        m = row.match(line)
        if m and not any(name in m.group(1) for name in NAMED_DISEASES):
            dropped += 1
            continue
        if dropped:
            out.append("  (other diseases withheld)")
            dropped = 0
        out.append(line)
    return "\n".join(out).rstrip() + "\n"


def phenotype(text: str) -> str:
    """Replace the table of the child's terms with the counts the report uses."""
    if "the child presents with" not in text:
        raise Unexpected("not a phenotype reverse check")
    out, skip = [], False
    for line in text.splitlines():
        if line.strip().startswith("the child presents with"):
            out.append("  (the child's feature names are withheld; counts only)")
            skip = True
            continue
        if skip:
            if not line.strip():
                skip = False
                out.append(line)
            continue
        out.append(line)
    return keep_named_disease_rows("\n".join(out))


def completeness(text: str) -> str:
    """Check that no rsID or gene-level detail is in the audit text; pass it through otherwise."""
    if "Which inheritance models" not in text:
        raise Unexpected("not a completeness audit")
    if re.search(r"\brs\d{4,}\b", text) or re.search(r"\b\d{1,2}:\d{4,}\b", text):
        raise Unexpected("an rsID or coordinate is present in the completeness audit")
    return text


def expanded(text: str) -> str:
    """Keep the per-gene counts; show carried P/LP alleles only for the reported gene."""
    if "ClinVar P/LP records inside the panel spans" not in text:
        raise Unexpected("not an expanded differential")
    out: list[str] = []
    gene = None
    in_carried = False
    for line in text.splitlines():
        if line.startswith("=== ClinVar Pathogenic/Likely-pathogenic alleles carried"):
            in_carried = True
            out.append(line)
            continue
        if not in_carried:
            out.append(line)
            continue
        m = re.match(r"^([A-Z0-9]+):$", line)
        if m:
            gene = m.group(1)
            out.append(line if gene in REPORTED_GENES else f"{gene}: (alleles withheld)")
            continue
        if line.startswith("   ") and gene is not None:
            if gene in REPORTED_GENES:
                out.append(re.sub(r"(DP=\d+).*$", r"\1", line))
            continue
        out.append(line)
    return "\n".join(out).rstrip() + "\n"


KINDS = {
    "ranking": ranking,
    "differential": differential,
    "mito": mito,
    "phenotype": phenotype,
    "completeness": completeness,
    "expanded": expanded,
}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--kind", choices=sorted(KINDS), required=True)
    args = ap.parse_args()
    try:
        sys.stdout.write(KINDS[args.kind](sys.stdin.read()))
    except Unexpected as err:
        print(f"public_view: {err}; refusing to pass the text through", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
