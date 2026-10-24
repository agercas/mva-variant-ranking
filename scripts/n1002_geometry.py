#!/usr/bin/env python3
"""Where does residue N1002 sit in the AlphaFold model of BubR1?

The report describes p.Asn1002Lys as a change to a mostly buried asparagine that makes
three hydrogen bonds to the protein backbone and lies far from the ATP site. This script
measures those things on the public AlphaFold model (UniProt O60566, version 6), so they
can be rechecked. No data from the child is used.

It reports:

* relative solvent accessibility (Shrake-Rupley, 1.4 A probe; maximum accessible surface
  area for Asn 195 A^2, Tien et al. 2013), and how that compares with the other Asn
  residues in the kinase domain (766-1050)
* the model confidence (pLDDT) at residue 1002
* distances to K795 (the lysine that binds ATP in active kinases) and to T1008 (a
  phosphorylation site)
* every polar contact of N1002 within 3.3 A to atoms outside residues 1001-1003

Usage::

    python3 scripts/n1002_geometry.py --download data/reference/alphafold

Needs Biopython (``pip install biopython``). An AlphaFold model is a prediction and not a
measurement; the numbers describe the model.
"""

from __future__ import annotations

import argparse
import statistics
import urllib.request
from pathlib import Path

from Bio.PDB import NeighborSearch, PDBParser, ShrakeRupley

MODEL = "AF-O60566-F1-model_v6.pdb"
URL = "https://alphafold.ebi.ac.uk/files/" + MODEL
MAX_ASA_ASN = 195.0  # Tien et al. 2013, theoretical maximum accessible surface area, A^2
SITE = 1002


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dir", type=Path, default=Path("data/reference/alphafold"))
    ap.add_argument("--download", type=Path, help="fetch the public model into this directory first")
    args = ap.parse_args()
    if args.download:
        args.download.mkdir(parents=True, exist_ok=True)
        req = urllib.request.Request(URL, headers={"User-Agent": "Mozilla/5.0 research-script"})
        (args.download / MODEL).write_bytes(urllib.request.urlopen(req, timeout=120).read())
        args.dir = args.download

    model = PDBParser(QUIET=True).get_structure("af", args.dir / MODEL)[0]
    chain = model["A"]
    ShrakeRupley(probe_radius=1.4, n_points=960).compute(model, level="R")
    res = chain[SITE]
    assert res.get_resname() == "ASN", f"residue {SITE} is {res.get_resname()}, expected ASN"

    def dist(a, b) -> float:
        return float(((a.coord - b.coord) ** 2).sum() ** 0.5)

    print(f"AlphaFold model {MODEL} (public; a prediction, not a measurement)\n")
    rsa = res.sasa / MAX_ASA_ASN
    domain = [r.sasa / MAX_ASA_ASN for r in chain if r.get_resname() == "ASN" and 766 <= r.id[1] <= 1050]
    print(f"N{SITE}: accessible surface {res.sasa:.1f} A^2, relative solvent accessibility {rsa:.2f}")
    print(
        f"  the {len(domain)} Asn residues in 766-1050 have median relative accessibility "
        f"{statistics.median(domain):.2f}; {sum(v > rsa for v in domain)} of {len(domain)} are more exposed than N{SITE}"
    )
    print(f"  model confidence (pLDDT) at residue {SITE}: {res['CA'].get_bfactor():.1f}\n")

    k795, t1008 = chain[795], chain[1008]
    print(
        f"distance from N{SITE} side chain (ND2) to the K795 side chain (NZ): {dist(res['ND2'], k795['NZ']):.1f} A"
    )
    print(f"distance between the alpha carbons of N{SITE} and T1008: {dist(res['CA'], t1008['CA']):.1f} A\n")

    polar = [a for r in chain for a in r if a.element in ("N", "O")]
    search = NeighborSearch(polar)
    print(f"polar contacts of N{SITE} within 3.3 A (outside residues {SITE - 1}-{SITE + 1}):")
    contacts = 0
    for atom in res:
        if atom.element not in ("N", "O"):
            continue
        for other in search.search(atom.coord, 3.3):
            partner = other.get_parent()
            if SITE - 1 <= partner.id[1] <= SITE + 1:
                continue
            contacts += 1
            print(
                f"  {atom.get_id():<4} -- {partner.get_resname()} {partner.id[1]} {other.get_id():<4} "
                f"{dist(atom, other):.2f} A"
            )
    backbone = sum(
        1
        for atom in res
        if atom.element in ("N", "O")
        for other in search.search(atom.coord, 3.3)
        if other.get_parent().id[1] not in (SITE - 1, SITE, SITE + 1) and other.get_id() in ("N", "O")
    )
    print(
        f"{contacts} contacts, {backbone} of them to backbone atoms. Lysine has no side-chain oxygen, so it"
    )
    print("cannot accept the hydrogen bond that the asparagine side-chain oxygen takes from a backbone N-H.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
