#!/usr/bin/env python3
"""Pick the most dissimilar (lowest-BLOSUM62) substitution at MHC-II anchor positions.

For a 9-mer MHC class II binding core, the P1/P4/P6/P9 pockets carry most of the
binding energy. This script reports, for each anchor, the amino acid with the
LOWEST BLOSUM62 score against the wild-type residue -- i.e. the most chemically
dissimilar replacement -- and builds the combined quadruple mutant.

Caveat worth stating in any write-up: BLOSUM62 scores evolutionary
interchangeability, not MHC binding. A low-BLOSUM substitution is a reasonable
first guess at "disrupt this anchor", but it is not the same as "lowest predicted
binding". Re-run NetMHCIIpan on the output to confirm the epitope is actually gone.

Examples
--------
  python worst_anchor_substitution.py --core SGGGLVQPG
  python worst_anchor_substitution.py --cores-file cores.txt --exclude CP --tsv out.tsv
"""
import argparse
import csv
import sys
from pathlib import Path

AA = "ARNDCQEGHILKMFPSTWYV"

# Standard BLOSUM62, row/column order as in AA above.
_B62_ROWS = """
 4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0
-1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3
-2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3
-2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3
 0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1
-1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2
-1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2
 0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3
-2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3
-1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3
-1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1
-1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2
-1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1
-2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1
-1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2
 1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2
 0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0
-3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3
-2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1
 0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4
"""

BLOSUM62 = {}
for _i, _line in enumerate(_B62_ROWS.strip().splitlines()):
    for _j, _val in enumerate(_line.split()):
        BLOSUM62[(AA[_i], AA[_j])] = int(_val)


def worst_substitutions(residue, allowed):
    """Return (min_score, [amino acids achieving it]) for `residue`."""
    scored = [(BLOSUM62[(residue, a)], a) for a in allowed if a != residue]
    if not scored:
        raise ValueError(f"No replacement available for {residue} with the current --exclude set")
    low = min(s for s, _ in scored)
    return low, sorted(a for s, a in scored if s == low)


def process_core(core, anchors, allowed):
    core = core.strip().upper()
    if len(core) != 9:
        raise ValueError(f"Core must be 9 residues, got {len(core)}: {core!r}")
    bad = sorted(set(core) - set(AA))
    if bad:
        raise ValueError(f"Core {core!r} contains non-standard residue(s): {','.join(bad)}")

    rows = []
    combined = list(core)
    for p in anchors:
        wt = core[p - 1]
        low, picks = worst_substitutions(wt, allowed)
        chosen = picks[0]
        combined[p - 1] = chosen
        rows.append({
            "core": core,
            "anchor": p,
            "wt": wt,
            "chosen": chosen,
            "blosum62": low,
            "self_score": BLOSUM62[(wt, wt)],
            "ties": "".join(picks),
            "single_mutant": core[:p - 1] + chosen + core[p:],
        })
    return rows, "".join(combined)


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--core", help="A single 9-mer binding core, e.g. SGGGLVQPG")
    src.add_argument("--cores-file",
                     help="File of cores, one per line; first whitespace field is used, "
                          "'#' comments and blank lines ignored")
    ap.add_argument("--anchors", default="1,4,6,9",
                    help="1-based anchor positions within the core (default: 1,4,6,9)")
    ap.add_argument("--exclude", default="C",
                    help="Residues never used as replacements (default: C, to protect "
                         "disulfides). Use '' to allow all 20; 'CP' also blocks proline.")
    ap.add_argument("--tsv", help="Also write results to this TSV file")
    args = ap.parse_args()

    try:
        anchors = [int(x) for x in args.anchors.split(",") if x.strip()]
    except ValueError:
        ap.error("--anchors must be comma-separated integers, e.g. 1,4,6,9")
    if not anchors or any(p < 1 or p > 9 for p in anchors):
        ap.error("--anchors values must be between 1 and 9")

    excluded = set(args.exclude.upper().replace(",", "").replace(" ", ""))
    if excluded - set(AA):
        ap.error("--exclude must contain only standard one-letter amino acid codes")
    allowed = [a for a in AA if a not in excluded]
    if not allowed:
        ap.error("--exclude leaves no replacement residues")

    if args.core:
        cores = [args.core]
    else:
        cores = []
        for line in Path(args.cores_file).read_text().splitlines():
            line = line.strip()
            if line and not line.startswith("#"):
                cores.append(line.split()[0])
        if not cores:
            ap.error(f"No cores found in {args.cores_file}")

    all_rows = []
    for core in cores:
        try:
            rows, combined = process_core(core, anchors, allowed)
        except ValueError as exc:
            print(f"SKIPPED {core}: {exc}", file=sys.stderr)
            continue
        all_rows.extend(rows)

        print(f"\nCore {rows[0]['core']}   (excluded as replacements: "
              f"{''.join(sorted(excluded)) or 'none'})")
        print(f"  {'pos':>4}  {'WT':^3}  {'->':^3}  {'B62':>4}  {'self':>4}  "
              f"{'ties':<8}  single mutant")
        for r in rows:
            print(f"  {'P' + str(r['anchor']):>4}  {r['wt']:^3}  {r['chosen']:^3}  "
                  f"{r['blosum62']:>4}  {r['self_score']:>4}  {r['ties']:<8}  "
                  f"{r['single_mutant']}")
        print(f"  all-{len(anchors)} combined mutant: {combined}")

    if args.tsv and all_rows:
        with open(args.tsv, "w", newline="") as fh:
            w = csv.DictWriter(fh, delimiter="\t", fieldnames=list(all_rows[0].keys()))
            w.writeheader()
            w.writerows(all_rows)
        print(f"\nWrote {len(all_rows)} rows to {args.tsv}")


if __name__ == "__main__":
    main()
