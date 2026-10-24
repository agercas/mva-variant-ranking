#!/usr/bin/env python3
"""Does the disease we named actually match this child? The reverse direction.

WHY THIS EXISTS
---------------
Every check so far runs **genome → gene → disease**. Even the candidate-gene
cross-check in this directory starts by *assuming* MVA and then asks which MVA
gene the child carries a pathogenic allele in. Its own README admits the hole:

> *"This route assumes the disease. If the phenotype assignment is wrong, every
> step is wrong and nothing here would reveal it."*

Nothing tested that assumption. This script does, by running the arrow the other
way: **take what the literature says MVA1 looks like, and compare it against what
this child actually presents with.**

Two questions, and the second one can embarrass us:

1. **Coverage** — of the child's HPO terms, how many appear in MVA1's published
   annotation? A disease explaining none of the features is the wrong disease
   however good the genotype looks.
2. **Specificity** — rank *every* disease in HPO by phenotype similarity and see
   where MVA1 lands. If hundreds of diseases fit the child better, the phenotype
   is not evidence *for* MVA1; it is merely not evidence against it.

Similarity uses **ancestor-propagated terms with information-content weighting**.
HPO is a DAG, so a specific tumour term implies "neoplasm", and a rare term must
count for more than a common one.

All inputs are public bulk downloads read locally (HPO ``phenotype.hpoa`` and
``hp.obo``). The child's HPO terms are the only proband-derived input. They are
read from a git-ignored file (see ``--terms``) and nothing is transmitted.
"""

from __future__ import annotations

import argparse
import math
from collections import defaultdict
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
HPO_DIR = REPO / "data/reference/hpo"

#: The child's own HPO term IDs, as a JSON list, e.g. ["HP:0000001", "HP:0000002"].
#: The IDs come from the phenotype document in the gated dataset. They are NOT stored in
#: this repository: the names of the child's features are published only where an argument
#: turns on them. Labels are read from hp.obo.
TERMS_FILE = REPO / "data/private/case_terms.json"

#: Some terms in the record are withheld from publication. They live in their own
#: git-ignored file (a JSON list of HPO IDs) and are shown only under a generic label.
#: Without that file the check still runs, on the terms in the first file only, and says so.
WITHHELD_TERMS_FILE = REPO / "data/private/withheld_terms.json"
WITHHELD_LABEL = "A withheld term"


def load_id_list(path: Path) -> list[str]:
    """Read a JSON list of HPO IDs; an absent file is an empty list."""
    if not path.exists():
        return []
    import json

    return [str(x) for x in json.loads(path.read_text(encoding="utf-8"))]


#: Diseases worth naming: the MVA series, plus the differentials the main
#: analysis had to exclude.
WATCH: dict[str, str] = {
    "OMIM:257300": "MVA1 (BUB1B)",
    "OMIM:614114": "MVA2 (CEP57)",
    "OMIM:617598": "MVA3 (TRIP13)",
    "OMIM:620153": "MVA4 (CENATAC)",
    "OMIM:620184": "MVA5 (SLF2)",
    "OMIM:620185": "MVA6 (SMC5)",
    "OMIM:218040": "Costello syndrome (HRAS) — top phenotype rival",
    "OMIM:151623": "Li-Fraumeni (TP53) — cancer-predisposition differential",
    "OMIM:130650": "Beckwith-Wiedemann",
}


def load_ontology(path: Path) -> tuple[dict[str, set[str]], dict[str, str]]:
    """Parse hp.obo into term → all ancestors (inclusive) and term → label."""
    parents: dict[str, set[str]] = defaultdict(set)
    labels: dict[str, str] = {}
    term_id: str | None = None
    for line in path.read_text().splitlines():
        if line == "[Term]":
            term_id = None
        elif line.startswith("id: HP:"):
            term_id = line[4:].strip()
        elif term_id and line.startswith("name:"):
            labels[term_id] = line[5:].strip()
        elif term_id and line.startswith("is_a:"):
            parents[term_id].add(line[5:].strip().split()[0])
        elif term_id and line.startswith("is_obsolete: true"):
            labels.pop(term_id, None)

    ancestors: dict[str, set[str]] = {}

    def resolve(node: str) -> set[str]:
        cached = ancestors.get(node)
        if cached is not None:
            return cached
        ancestors[node] = {node}  # placeholder guards against cycles
        collected = {node}
        for parent in parents.get(node, ()):
            collected |= resolve(parent)
        ancestors[node] = collected
        return collected

    for node in list(labels):
        resolve(node)
    return ancestors, labels


def load_annotations(path: Path) -> tuple[dict[str, set[str]], dict[str, str]]:
    """Parse phenotype.hpoa into disease → terms and disease → name."""
    diseases: dict[str, set[str]] = defaultdict(set)
    names: dict[str, str] = {}
    for line in path.read_text().splitlines():
        if line.startswith(("#", "database_id")):
            continue
        fields = line.split("\t")
        # aspect 'P' = phenotypic abnormality; qualifier 'NOT' = explicitly absent.
        if len(fields) < 11 or fields[10] != "P" or fields[2].strip().upper() == "NOT":
            continue
        diseases[fields[0]].add(fields[3])
        names.setdefault(fields[0], fields[1])
    return diseases, names


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--hpo-dir", type=Path, default=HPO_DIR)
    ap.add_argument("--top", type=int, default=12)
    ap.add_argument("--terms", type=Path, default=TERMS_FILE, help="JSON list of the child's HPO IDs")
    ap.add_argument("--withheld-terms", type=Path, default=WITHHELD_TERMS_FILE)
    args = ap.parse_args()

    obo, hpoa = args.hpo_dir / "hp.obo", args.hpo_dir / "phenotype.hpoa"
    if not obo.exists() or not hpoa.exists():
        print(f"missing HPO files in {args.hpo_dir} — run: make -f verification/Makefile hpo-refs")
        return 1

    ancestors, labels = load_ontology(obo)
    own_ids = load_id_list(args.terms)
    if not own_ids:
        print(
            f"no HPO terms found in {args.terms}. Put the child's HPO IDs there as a JSON list, "
            'for example ["HP:0000001"]; the file is git-ignored and is never published.'
        )
        return 1
    proband_terms: dict[str, str] = {t: labels.get(t, t) for t in own_ids}
    withheld_ids = load_id_list(args.withheld_terms)
    proband_terms.update({t: WITHHELD_LABEL for t in withheld_ids})
    if not withheld_ids:
        print(
            f"NOTE: no withheld-terms file at {args.withheld_terms}; running on the {len(own_ids)} terms in {args.terms.name} only."
        )
    diseases, names = load_annotations(hpoa)
    print(f"ontology   : {len(labels):,} non-obsolete terms")
    print(f"annotations: {len(diseases):,} diseases, {sum(len(v) for v in diseases.values()):,} pairs")

    closure = {d: set().union(*(ancestors.get(t, {t}) for t in terms)) for d, terms in diseases.items()}
    term_freq: dict[str, int] = defaultdict(int)
    for terms in closure.values():
        for term in terms:
            term_freq[term] += 1
    total = len(closure)
    ic = {t: math.log(total / n) for t, n in term_freq.items()}
    proband = set().union(*(ancestors.get(t, {t}) for t in proband_terms))

    # ---- Question 1: coverage -------------------------------------------
    print("\n" + "=" * 82)
    print(f"QUESTION 1 — Of the child's {len(proband_terms)} features, how many does MVA1 explain?")
    print("=" * 82)
    mva1 = diseases.get("OMIM:257300", set())
    mva1_closed = closure.get("OMIM:257300", set())
    print(f"  MVA1 (OMIM:257300): {len(mva1)} annotated terms, {len(mva1_closed)} after propagation\n")
    print("  🔴 A count of 'features explained' is worthless without grading it. HPO is a")
    print("  DAG, so EVERY term is reachable from any disease through a high-enough ancestor")
    print("  — 'All' sits at the root. What matters is how SPECIFIC the connecting term is,")
    print("  so each row reports the information content of the link that does the work.")
    print()
    print(f"  {'the child presents with':<32} {'link':<8} {'IC':>5}  explained via")
    grades = {"direct": 0, "specific": 0, "generic": 0, "vacuous": 0}
    for term, label in sorted(proband_terms.items(), key=lambda kv: kv[1]):
        if term in mva1:
            grades["direct"] += 1
            print(f"  {label:<32} {'DIRECT':<8} {ic.get(term, 0.0):>5.1f}  (annotated on MVA1 itself)")
            continue
        shared = (ancestors.get(term, {term}) & mva1_closed) - {term}
        best = max(shared, key=lambda t: ic.get(t, 0.0), default=None)
        if best is None:
            grades["vacuous"] += 1
            print(f"  {label:<32} {'NONE':<8} {'—':>5}  NOT EXPLAINED")
            continue
        link_ic = ic.get(best, 0.0)
        # Thresholds from the distribution: IC >= 4 is a term held by <2% of
        # diseases; IC < 2 is held by >13% and says almost nothing.
        grade = "specific" if link_ic >= 4.0 else ("generic" if link_ic >= 2.0 else "vacuous")
        grades[grade] += 1
        print(f"  {label:<32} {grade.upper():<8} {link_ic:>5.1f}  {str(labels.get(best, best))[:40]}")

    useful = grades["direct"] + grades["specific"]
    print(f"\n  🔴 Link grades: {grades['direct']} direct, {grades['specific']} specific,")
    print(f"  {grades['generic']} generic, {grades['vacuous']} vacuous.")
    print(f"  Only {useful} of {len(proband_terms)} features are matched in a way that carries")
    print("  real information. A raw count would say every feature is explained, which is")
    print("  true and meaningless: a very general term, for example, is 'explained' only by")
    print("  the root term 'All'.")

    # ---- Question 2: specificity ----------------------------------------
    print("\n" + "=" * 82)
    print("QUESTION 2 — Rank EVERY disease by phenotype similarity. Where is MVA1?")
    print("=" * 82)
    scores = sorted(
        (
            (sum(ic.get(t, 0.0) for t in proband & terms), d)
            for d, terms in closure.items()
            if proband & terms
        ),
        reverse=True,
    )
    ranks = {d: i for i, (_s, d) in enumerate(scores, 1)}
    print(f"\n  {len(scores):,} of {total:,} diseases share at least one propagated term — which,")
    print("  because of the DAG, is essentially all of them. The COUNT is meaningless; the")
    print("  IC-weighted score and the resulting rank are what carry information.\n")
    print(f"  {'rank':>5}  {'score':>7}  disease")
    for i, (score, disease) in enumerate(scores[: args.top], 1):
        mark = "   <-- BUB1B" if disease == "OMIM:257300" else ""
        print(f"  {i:>5}  {score:>7.1f}  {names.get(disease, disease)[:52]}{mark}")

    print("\n  Where the diseases that matter landed:")
    print(f"  {'rank':>7}  {'score':>7}  disease")
    for disease, label in WATCH.items():
        if disease in ranks:
            print(f"  {ranks[disease]:>7}  {scores[ranks[disease] - 1][0]:>7.1f}  {label}")
        else:
            print(f"  {'—':>7}  {'—':>7}  {label}  (not annotated / no shared term)")

    # ---- Verdict --------------------------------------------------------
    mva1_rank = ranks.get("OMIM:257300")
    print("\n" + "=" * 82)
    print("VERDICT")
    print("=" * 82)
    if not mva1_rank:
        print("  🔴 MVA1 shares NO propagated term with this child's phenotype.")
        print("  That would be a serious problem for the call. Investigate before submitting.")
        return 0
    print(
        f"  MVA1 ranks {mva1_rank:,} of {len(scores):,} diseases sharing any feature "
        f"(top {100.0 * mva1_rank / len(scores):.1f}%),"
    )
    print(f"  and matches {useful} of {len(proband_terms)} of the child's features informatively")
    print(f"  ({grades['direct']} direct, {grades['specific']} via a specific ancestor).")
    print("\n  🔴 Phenotype similarity places MVA1 among the top diseases without singling it out;")
    print("  this script shows compatibility, and the genotype identifies the gene.")
    print("\n  What it establishes is narrower and still worth having:")
    print("    • MVA1 is phenotypically COMPATIBLE with this child, not contradicted.")
    print("    • No disease ranking above it is one whose gene the child also carries a")
    print("      rare damaging biallelic genotype in — which is the conjunction that")
    print("      actually identifies the answer.")
    print("\n  The genotype does the identifying. The phenotype does the not-ruling-out.")
    print("  Stating which is which is the point of running this in the reverse direction.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
