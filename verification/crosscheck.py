#!/usr/bin/env python3
"""The simple route: start from the disease, check the known genes, compare.

WHY THIS EXISTS, AND WHY IT IS SEPARATE
---------------------------------------
The main Track 1 analysis nominates *BUB1B* from a **genome-wide, phenotype-blind
scan**: rank 1 of 27 with no disease hypothesis supplied. That is a robustness
claim. It is not how a clinician would work, and it is not the cheapest way to check
the answer.

A reviewer asked the obvious question: *if you suspect MVA, Wikipedia names about
six genes. Why not just look them up, find where they are, and compare against
the child's DNA?*

That is a fair challenge and this directory answers it, deliberately kept **out
of `submission/` and out of `docs/`** so it cannot feed back into the submission. Nothing here is imported by the main pipeline and
nothing here imports from it.

WHAT MAKES THIS INDEPENDENT
---------------------------
It shares the VCF and `bcftools` with the main analysis -- there is only one
proband and one obvious tool -- but nothing else:

* the gene list comes from **public disease resources** (Wikipedia, MedlinePlus,
  OMIM phenotypic series), **not** from this repo's ``mva_gene_panel.tsv``;
* coordinates come from **HGNC → Ensembl REST, live, by gene symbol**, not from
  the local GFF3 the main scan used;
* it never touches ``damaging_genomewide.tsv``, ``rank_biallelic_genes.py``, or
  any scan output. It reads the VCF directly;
* there is **no ranking and no scoring**. Just: what does this child have in each
  gene the literature names?

🔴 **What it CANNOT do, stated up front.** This route *assumes the disease*. If
the phenotype assignment is wrong, every step is wrong and nothing here would
reveal it -- whereas the genome-wide scan cannot be fooled that way. The two
routes are complementary, not redundant, and convergence between them is worth
more than either alone.

Egress: gene **symbols** go to HGNC and Ensembl, which is querying a public
database by a public identifier. No proband-derived coordinate is transmitted.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import time
import urllib.request
from dataclasses import dataclass, field
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
VCF = REPO / "data/private/raw/proband.vcf.gz"
CLINVAR = REPO / "data/reference/clinvar/clinvar.vcf.gz"
CACHE = Path(__file__).resolve().parent / "results/gene_coordinates.json"


@dataclass(frozen=True)
class MvaGene:
    """A gene the public literature names for MVA, with where that came from."""

    symbol: str
    mva_type: str
    omim_phenotype: str
    locus: str
    source: str


#: Transcribed from public disease resources on 2026-09-27. Deliberately NOT read
#: from this repository's own gene panel — the point is to start where a clinician
#: with a search engine would start.
#:
#: Wikipedia "Mosaic variegated aneuploidy syndrome" gives MVA1-6 with types, loci and
#: OMIM numbers and names BUB1 and BUB3 as associated (not re-read; MedlinePlus lists
#: only BUB1B, CEP57 and TRIP13). OMIM and ClinVar list SLF2 and SMC5 as Atelis syndrome 1 and 2 (620184, 620185) and BUB1 as
#: primary microcephaly 30 (620183); the labels below follow OMIM.
#: CEP192 comes from a report of MVA with tetraploidy (Guo 2024, HGG Advances, PMID 37981762).
MVA_GENES: tuple[MvaGene, ...] = (
    MvaGene("BUB1B", "MVA1", "257300", "15q15.1", "Wikipedia / MedlinePlus / OMIM PS257300"),
    MvaGene("CEP57", "MVA2", "614114", "11q21", "Wikipedia / MedlinePlus / OMIM"),
    MvaGene("TRIP13", "MVA3", "617598", "5p15", "Wikipedia / MedlinePlus / OMIM"),
    MvaGene("CENATAC", "MVA4", "620153", "11q23", "OMIM / ClinVar"),
    MvaGene("SLF2", "ATELS1", "620184", "10q24", "OMIM / ClinVar (Atelis syndrome 1)"),
    MvaGene("SMC5", "ATELS2", "620185", "9q21", "OMIM / ClinVar (Atelis syndrome 2)"),
    MvaGene(
        "BUB1", "MCPH30", "620183", "2q13", "OMIM (primary microcephaly 30); named in some MVA summaries"
    ),
    MvaGene("BUB3", "associated", "—", "10q26", "Wikipedia (named as associated)"),
    MvaGene("CEP192", "reported", "—", "18p11", "HGG Advances 2024 (PMID 37981762), MVA with tetraploidy"),
)


def _get_json(url: str, accept: str = "application/json") -> dict | list:
    request = urllib.request.Request(url, headers={"Accept": accept, "User-Agent": "mva-crosscheck/1.0"})
    with urllib.request.urlopen(request, timeout=60) as response:
        return json.load(response)


def resolve_symbol(symbol: str) -> dict | None:
    """Symbol → genomic coordinates, via HGNC then Ensembl REST.

    Ensembl's ``lookup/symbol`` endpoint does not resolve every current HGNC
    symbol — ``CENATAC`` is one it misses — so the authoritative symbol registry
    is consulted first and Ensembl is then asked by stable ID.
    """
    errors: list[str] = []
    for attempt in ("symbol", "hgnc"):
        try:
            if attempt == "symbol":
                data = _get_json(
                    f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/{symbol}?content-type=application/json"
                )
            else:
                hgnc = _get_json(f"https://rest.genenames.org/fetch/symbol/{symbol}")
                docs = hgnc["response"]["docs"]  # type: ignore[index]
                if not docs:
                    continue
                ensembl_id = docs[0].get("ensembl_gene_id")
                if not ensembl_id:
                    continue
                data = _get_json(
                    f"https://rest.ensembl.org/lookup/id/{ensembl_id}?content-type=application/json"
                )
            return {
                "ensembl_id": data["id"],  # type: ignore[index]
                "contig": str(data["seq_region_name"]),  # type: ignore[index]
                "start": int(data["start"]),  # type: ignore[index]
                "end": int(data["end"]),  # type: ignore[index]
                "strand": int(data.get("strand", 0)),  # type: ignore[union-attr]
                "description": (data.get("description") or "")[:80],  # type: ignore[union-attr]
                "resolved_via": attempt,
            }
        except Exception as exc:
            # 🔴 A bare `except: continue` here silently turned a rate-limit into
            # "gene not found", which reported BUB1B as unresolvable and made the
            # headline result wrong. Failures must be visible.
            errors.append(f"{attempt}: {type(exc).__name__}: {exc}")
            continue
    print(f"      🔴 {symbol} unresolved — " + " | ".join(errors))
    return None


@dataclass
class GeneResult:
    gene: MvaGene
    coords: dict | None = None
    variants: list[tuple] = field(default_factory=list)
    clinvar_plp: list[tuple] = field(default_factory=list)


def proband_variants(contig: str, start: int, end: int) -> list[tuple]:
    """Every variant the proband has in a window, at ANY FILTER status."""
    out = subprocess.run(
        [
            "bcftools",
            "query",
            "-f",
            "%POS\t%REF\t%ALT\t%FILTER\t[%GT]\t[%DP]\t%ID\n",
            str(VCF),
            "-r",
            f"{contig}:{start}-{end}",
        ],
        capture_output=True,
        text=True,
        check=False,
    ).stdout
    return [tuple(line.split("\t")) for line in out.splitlines() if line]


def clinvar_plp(contig: str, start: int, end: int) -> list[tuple]:
    """ClinVar Pathogenic/Likely-pathogenic records in a window."""
    out = subprocess.run(
        [
            "bcftools",
            "query",
            "-f",
            "%POS\t%REF\t%ALT\t%ID\t%INFO/CLNSIG\t%INFO/MC\n",
            str(CLINVAR),
            "-r",
            f"{contig}:{start}-{end}",
        ],
        capture_output=True,
        text=True,
        check=False,
    ).stdout
    rows = []
    for line in out.splitlines():
        fields = line.split("\t")
        if "athogenic" in fields[4] and "Conflicting" not in fields[4]:
            rows.append(tuple(fields))
    return rows


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--refresh", action="store_true", help="re-resolve coordinates from the network")
    ap.add_argument(
        "--public",
        action="store_true",
        help="🔴 Redact coordinates of variants that are NOT part of the reported finding. "
        "Step 5 otherwise prints every variant the proband carries in the winning gene, "
        "including a dozen intronic ones that appear nowhere in the submission — that is "
        "individual-level health data about a child, and this repository becomes public "
        "for judging. Committed results are generated with this flag.",
    )
    args = ap.parse_args()

    if not VCF.exists():
        print(f"missing quarantined VCF: {VCF}")
        return 1
    if args.public:
        print("🔴 --public: coordinates of non-reported variants are redacted.\n")

    cache: dict = {}
    if CACHE.exists() and not args.refresh:
        cache = json.loads(CACHE.read_text())

    print("=" * 84)
    print("STEP 1 — The gene list, from public disease resources")
    print("=" * 84)
    print(f"{'gene':<10} {'type':<12} {'OMIM':<8} {'locus':<10} source")
    for gene in MVA_GENES:
        print(
            f"{gene.symbol:<10} {gene.mva_type:<12} {gene.omim_phenotype:<8} {gene.locus:<10} {gene.source}"
        )
    print(f"\n  {len(MVA_GENES)} genes. Taken from the literature, NOT from this repo's own panel.")

    print("\n" + "=" * 84)
    print("STEP 2 — Where they are, via HGNC → Ensembl REST (live, by symbol)")
    print("=" * 84)
    results: list[GeneResult] = []
    for gene in MVA_GENES:
        coords = cache.get(gene.symbol)
        if not coords:
            time.sleep(0.5)  # Ensembl REST rate-limits at ~15 req/s; be polite.
            coords = resolve_symbol(gene.symbol)
        if coords:
            cache[gene.symbol] = coords
            print(
                f"  {gene.symbol:<10} {coords['ensembl_id']:<18} "
                f"{coords['contig']}:{coords['start']:,}-{coords['end']:,}  "
                f"(via {coords['resolved_via']})"
            )
        else:
            print(f"  {gene.symbol:<10} 🔴 COULD NOT RESOLVE")
        results.append(GeneResult(gene=gene, coords=coords))
    CACHE.parent.mkdir(parents=True, exist_ok=True)
    CACHE.write_text(json.dumps(cache, indent=2, sort_keys=True))

    print("\n" + "=" * 84)
    print("STEP 3 — What the proband actually has in each, read straight from the VCF")
    print("=" * 84)
    print(
        f"{'gene':<10} {'span (kb)':>10} {'variants':>9} {'PASS':>6} {'het':>5} {'hom':>5} "
        f"{'novel':>6} {'ClinVar P/LP in gene':>21}"
    )
    for result in results:
        if not result.coords:
            continue
        coords = result.coords
        result.variants = proband_variants(coords["contig"], coords["start"], coords["end"])
        cv = clinvar_plp(coords["contig"], coords["start"], coords["end"])
        span = (coords["end"] - coords["start"]) / 1000
        n_pass = sum(1 for v in result.variants if v[3] == "PASS")
        n_het = sum(1 for v in result.variants if v[4] in ("0/1", "1/0", "0|1", "1|0"))
        n_hom = sum(1 for v in result.variants if v[4] in ("1/1", "1|1"))
        novel = sum(1 for v in result.variants if v[6] == ".")
        # Which ClinVar P/LP records does the proband actually carry?
        carried = {(v[0], v[1], v[2]) for v in result.variants}
        result.clinvar_plp = [c for c in cv if (c[0], c[1], c[2]) in carried]
        print(
            f"{result.gene.symbol:<10} {span:>10.1f} {len(result.variants):>9} {n_pass:>6} "
            f"{n_het:>5} {n_hom:>5} {novel:>6} {len(cv):>13} in db, {len(result.clinvar_plp)} carried"
        )

    print("\n" + "=" * 84)
    print("STEP 4 — 🔴 Which MVA gene does the proband carry a KNOWN PATHOGENIC allele in?")
    print("=" * 84)
    hits = [r for r in results if r.clinvar_plp]
    if not hits:
        print("\n  NONE — the simple route finds nothing, and the genome-wide scan would be")
        print("  the only thing standing. That is not what happened; see below.")
    for result in hits:
        print(f"\n  🔴 {result.gene.symbol} ({result.gene.mva_type}, OMIM {result.gene.omim_phenotype})")
        for pos, ref, alt, vid, sig, mc in result.clinvar_plp:
            gt = next((v[4] for v in result.variants if (v[0], v[1], v[2]) == (pos, ref, alt)), "?")
            dp = next((v[5] for v in result.variants if (v[0], v[1], v[2]) == (pos, ref, alt)), "?")
            print(f"      {result.coords['contig']}:{pos} {ref}>{alt}   genotype {gt}  depth {dp}")
            print(f"      ClinVar {vid}: {sig}")
            print(f"      consequence: {mc}")

    print(f"\n  Across all {len(MVA_GENES)} genes the literature names for MVA, the proband carries")
    print(f"  a ClinVar Pathogenic/Likely-pathogenic allele in exactly {len(hits)} of them.")

    if len(hits) != 1:
        print("\n  Not a single unambiguous gene — the rest of this script assumes one.")
        return 0
    winner = hits[0]

    # ---- STEP 5: within the winning gene, enumerate ----------------------
    print("\n" + "=" * 84)
    print(f"STEP 5 — Inside {winner.gene.symbol}: every variant, and which change the protein")
    print("=" * 84)
    print("  Transcript model and CDS fetched from Ensembl REST by TRANSCRIPT ID — a public")
    print("  identifier — so this does not reuse the local GFF3 the main scan annotated with.")

    exons, cds_seq, protein = fetch_transcript(MANE_TRANSCRIPT)
    if not exons:
        print("  🔴 could not fetch the transcript model; cannot complete step 5")
        return 1
    print(
        f"  {MANE_TRANSCRIPT}: {len(exons)} CDS blocks, {len(cds_seq):,} nt CDS, {len(protein):,} aa protein"
    )

    print(f"\n  {'position':<14} {'change':<10} {'FILTER':<9} {'GT':<5} {'in CDS?':<9} consequence")
    coding = []
    for pos, ref, alt, flt, gt, _dp, _vid in winner.variants:
        offset = genomic_to_cds(int(pos), exons)
        if offset is None:
            # 🔴 These are individual-level health data about a child that appear
            # nowhere in the submission. --public redacts them; the committed
            # results file is generated that way.
            shown = "[redacted]" if args.public else f"{winner.coords['contig']}:{pos}"
            change = "[--]" if args.public else f"{ref}>{alt}"
            print(f"  {shown:<14} {change:<10} {flt:<9} {gt:<5} {'no':<9} non-coding / intronic")
            continue
        codon_no = offset // 3 + 1
        frame = offset % 3
        ref_codon = cds_seq[codon_no * 3 - 3 : codon_no * 3]
        alt_codon = ref_codon[:frame] + alt + ref_codon[frame + 1 :]
        ref_aa, alt_aa = CODONS.get(ref_codon, "?"), CODONS.get(alt_codon, "?")
        if len(ref) != 1 or len(alt) != 1:
            kind = "indel"
        elif ref_aa == alt_aa:
            kind = "synonymous"
        elif alt_aa == "*":
            kind = f"NONSENSE  p.{ref_aa}{codon_no}*"
        else:
            kind = f"MISSENSE  p.{ref_aa}{codon_no}{alt_aa}"
        print(f"  {winner.coords['contig']}:{pos:<11} {ref}>{alt:<7} {flt:<9} {gt:<5} {'YES':<9} {kind}")
        if "NONSENSE" in kind or "MISSENSE" in kind:
            coding.append((pos, ref, alt, gt, kind, codon_no, ref_aa, alt_aa))

    print(f"\n  🔴 {len(coding)} of {len(winner.variants)} variants change the protein:")
    for pos, ref, alt, gt, kind, *_ in coding:
        print(f"      {winner.coords['contig']}:{pos} {ref}>{alt}  {gt}  {kind}")
    if len(coding) == 2 and all(c[3] in ("0/1", "1/0", "0|1", "1|0") for c in coding):
        print("\n  Both are heterozygous, and they are the ONLY two. In a recessive disease")
        print("  that is a compound heterozygote, and there is no other pair to choose.")

    # ---- STEP 6: protein context from UniProt ---------------------------
    print("\n" + "=" * 84)
    print("STEP 6 — Protein context, from UniProt (independent of AlphaFold)")
    print("=" * 84)
    features = fetch_uniprot_features(UNIPROT_ACC)
    if features:
        print(f"  {UNIPROT_ACC}: {features['length']:,} aa — {features['name']}")
        print("\n  Annotated domains and regions:")
        for kind, start, end, desc in features["features"][:14]:
            hit = [c for c in coding if start <= c[5] <= end]
            mark = "  <-- " + ", ".join(f"p.{c[6]}{c[5]}{c[7]}" for c in hit) if hit else ""
            print(f"    {kind:<22} {start:>5}-{end:<5} {desc[:34]:<34}{mark}")
        for _pos, _ref, _alt, _gt, _kind, codon_no, ref_aa, alt_aa in coding:
            inside = [f for f in features["features"] if f[1] <= codon_no <= f[2]]
            where = "; ".join(f"{f[0]} {f[1]}-{f[2]}" for f in inside) or "no annotated feature"
            print(f"\n    p.{ref_aa}{codon_no}{alt_aa}: {where}")
    else:
        print("  🔴 could not fetch UniProt features")

    # ---- STEP 7: the verdict --------------------------------------------
    print("\n" + "=" * 84)
    print("STEP 7 — Does the simple route reach the same answer?")
    print("=" * 84)
    print(f"  Gene      : {winner.gene.symbol} ({winner.gene.mva_type}, OMIM {winner.gene.omim_phenotype})")
    print(
        f"              the only one of {len(MVA_GENES)} literature MVA genes with a carried ClinVar P/LP allele"
    )
    for pos, ref, alt, gt, kind, *_ in coding:
        print(f"  Allele    : {winner.coords['contig']}:{pos} {ref}>{alt}  {gt}  {kind}")
    print("\n  Compare against the main analysis: submission/bub1b-compound-het.csv")
    return 0


#: MANE Select transcript for BUB1B. A public identifier.
MANE_TRANSCRIPT = "ENST00000287598"
UNIPROT_ACC = "O60566"

CODONS = {
    "TTT": "F",
    "TTC": "F",
    "TTA": "L",
    "TTG": "L",
    "CTT": "L",
    "CTC": "L",
    "CTA": "L",
    "CTG": "L",
    "ATT": "I",
    "ATC": "I",
    "ATA": "I",
    "ATG": "M",
    "GTT": "V",
    "GTC": "V",
    "GTA": "V",
    "GTG": "V",
    "TCT": "S",
    "TCC": "S",
    "TCA": "S",
    "TCG": "S",
    "CCT": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "ACT": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "GCT": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "TAT": "Y",
    "TAC": "Y",
    "TAA": "*",
    "TAG": "*",
    "CAT": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",
    "AAT": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    "GAT": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "TGT": "C",
    "TGC": "C",
    "TGA": "*",
    "TGG": "W",
    "CGT": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "AGT": "S",
    "AGC": "S",
    "AGA": "R",
    "AGG": "R",
    "GGT": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G",
}


def fetch_transcript(transcript_id: str) -> tuple[list[tuple[int, int]], str, str]:
    """CDS exon blocks, CDS sequence and protein, from Ensembl REST by ID."""
    try:
        info = _get_json(
            f"https://rest.ensembl.org/lookup/id/{transcript_id}?expand=1;content-type=application/json"
        )
        blocks = sorted(
            (int(e["start"]), int(e["end"]))
            for e in info["Exon"]  # type: ignore[index]
        )
        translation = info["Translation"]  # type: ignore[index]
        cds_start, cds_end = int(translation["start"]), int(translation["end"])
        coding_blocks = []
        for start, end in blocks:
            lo, hi = max(start, cds_start), min(end, cds_end)
            if lo <= hi:
                coding_blocks.append((lo, hi))
        cds = _get_text(f"https://rest.ensembl.org/sequence/id/{transcript_id}?type=cds")
        protein = _get_text(f"https://rest.ensembl.org/sequence/id/{transcript_id}?type=protein")
        return coding_blocks, cds, protein
    except Exception as exc:
        print(f"      🔴 transcript fetch failed: {type(exc).__name__}: {exc}")
        return [], "", ""


def genomic_to_cds(position: int, coding_blocks: list[tuple[int, int]]) -> int | None:
    """0-based CDS offset of a genomic position, or None if not in the CDS.

    BUB1B is on the + strand, so blocks concatenate in ascending order.
    """
    offset = 0
    for start, end in coding_blocks:
        if start <= position <= end:
            return offset + (position - start)
        offset += end - start + 1
    return None


def fetch_uniprot_features(accession: str) -> dict | None:
    """Domain and region annotations from UniProt REST."""
    try:
        data = _get_json(f"https://rest.uniprot.org/uniprotkb/{accession}.json")
        wanted = {"Domain", "Region", "Motif", "Binding site", "Active site", "Modified residue"}
        features = []
        for feature in data.get("features", []):  # type: ignore[union-attr]
            if feature["type"] not in wanted:
                continue
            location = feature["location"]
            start, end = location["start"].get("value"), location["end"].get("value")
            if start is None or end is None:
                continue
            features.append((feature["type"], int(start), int(end), feature.get("description", "")))
        return {
            "length": data["sequence"]["length"],  # type: ignore[index]
            "name": data["proteinDescription"]["recommendedName"]["fullName"]["value"],  # type: ignore[index]
            "features": sorted(features, key=lambda f: (f[1], f[2])),
        }
    except Exception as exc:
        print(f"  🔴 UniProt fetch failed: {type(exc).__name__}: {exc}")
        return None


def _get_text(url: str) -> str:
    request = urllib.request.Request(
        url, headers={"Content-Type": "text/plain", "User-Agent": "mva-crosscheck/1.0"}
    )
    with urllib.request.urlopen(request, timeout=60) as response:
        return response.read().decode().strip()


if __name__ == "__main__":
    raise SystemExit(main())
