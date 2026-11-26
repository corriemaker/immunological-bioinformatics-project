#!/usr/bin/env python3
"""Build paired VH/VL FASTA files for ABodyBuilder from a mutant chain FASTA.

ABodyBuilder models the Fv, so each input needs a heavy and a light variable
domain together. The mutant FASTA holds full chains (VH+CH1, VL+CL) and each
record varies only one of them, so for every selected mutant this script trims
the mutated chain to its variable domain and pairs it with the WILD-TYPE partner
domain taken from the same file.

Variable-domain boundaries are the variable-region sequences SAbDab publishes for
4G3Y: heavy 1-120 (ends ...WGQGTTLTVSS), light 1-107 (ends ...FGSGTNLEVK). They
are verified against the WT records at run time rather than trusted blindly.

Output: one FASTA per mutant containing records 'H' and 'L'.

Examples
--------
  python make_abodybuilder_input.py --fasta mutants_pssm.fasta \
      --map mutants_pssm.tsv --records M08,M09,M18,M51 --outdir abb_input

  python make_abodybuilder_input.py --fasta mutants_pssm.fasta \
      --map mutants_pssm.tsv --summary step6_pssm_summary.tsv --outdir abb_input
"""
import argparse
import csv
import re
import sys
from pathlib import Path

HEAVY_CHAIN = "4G3Y_2_Chain"
LIGHT_CHAIN = "4G3Y_1_Chain"
VH_END, VL_END = 120, 107
VH_TAIL, VL_TAIL = "WGQGTTLTVSS", "FGSGTNLEVK"

# SAbDab IMGT CDRs, as sequential positions in prot_AB.fasta
CDRS = {
    HEAVY_CHAIN: [("CDR-H1", 26, 33), ("CDR-H2", 51, 60), ("CDR-H3", 99, 109)],
    LIGHT_CHAIN: [("CDR-L1", 27, 32), ("CDR-L2", 50, 52), ("CDR-L3", 89, 97)],
}


def read_fasta(path):
    seqs, name, parts = {}, None, []
    for raw in Path(path).read_text().splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name:
                seqs[name] = "".join(parts).upper()
            name, parts = line[1:].split()[0], []
        else:
            parts.append(line)
    if name:
        seqs[name] = "".join(parts).upper()
    return seqs


def write_fasta(path, records, width=60):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as fh:
        for header, seq in records:
            fh.write(f">{header}\n")
            for i in range(0, len(seq), width):
                fh.write(seq[i:i + width] + "\n")


def region_of(chain, pos):
    for name, lo, hi in CDRS.get(chain, []):
        if lo <= pos <= hi:
            return name
    return "framework"


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fasta", required=True, help="mutants_pssm.fasta")
    ap.add_argument("--map", dest="mapfile", required=True, help="mutants_pssm.tsv")
    ap.add_argument("--records", help="Comma-separated record ids, e.g. M08,M09,M18,M51")
    ap.add_argument("--summary",
                    help="step6 summary TSV; selects CLEARED rows with no new binders")
    ap.add_argument("--include-cdr", action="store_true",
                    help="With --summary, keep hits inside a CDR (excluded by default)")
    ap.add_argument("--outdir", default="abodybuilder_input")
    args = ap.parse_args()

    seqs = read_fasta(args.fasta)
    with open(args.mapfile, newline="") as fh:
        mapping = {r["record_id"]: r for r in csv.DictReader(fh, delimiter="\t")}

    # wild-type variable domains
    wt = {}
    for chain, end, tail in ((HEAVY_CHAIN, VH_END, VH_TAIL), (LIGHT_CHAIN, VL_END, VL_TAIL)):
        wt_id = next((r["wt_record_id"] for r in mapping.values()
                      if r["chain"] == chain and r.get("wt_record_id")), None)
        if wt_id is None or wt_id not in seqs:
            sys.exit(f"No wild-type record for {chain} in {args.fasta}")
        domain = seqs[wt_id][:end]
        if not domain.endswith(tail):
            sys.exit(f"{chain}: residues 1-{end} do not end in {tail!r} -- "
                     "the variable-domain boundary is wrong for this FASTA")
        wt[chain] = domain
    print(f"WT domains verified: VH {len(wt[HEAVY_CHAIN])} aa, VL {len(wt[LIGHT_CHAIN])} aa")

    # which records
    if args.records:
        wanted = [r.strip() for r in args.records.split(",") if r.strip()]
    elif args.summary:
        wanted = []
        with open(args.summary, newline="") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                if row["verdict"] != "CLEARED" or int(row.get("n_new") or 0):
                    continue
                m = re.search(r"\d+", row["substitutions"])
                if not m:
                    continue
                pos, chain = int(m.group()), row["chain"]
                end = VH_END if chain == HEAVY_CHAIN else VL_END
                if pos > end:
                    continue
                if not args.include_cdr and region_of(chain, pos) != "framework":
                    continue
                wanted.append(row["record_id"])
    else:
        ap.error("give --records or --summary")

    outdir, written = Path(args.outdir), 0
    for rec in wanted:
        entry = mapping.get(rec)
        if entry is None:
            print(f"{rec}: not in the map file", file=sys.stderr)
            continue
        if rec not in seqs:
            print(f"{rec}: not in the FASTA", file=sys.stderr)
            continue

        chain = entry["chain"]
        end = VH_END if chain == HEAVY_CHAIN else VL_END
        mutated = seqs[rec][:end]

        pos = int(re.search(r"\d+", entry["substitutions"]).group())
        if pos > end:
            print(f"{rec}: substitution at {pos} is outside the variable domain; "
                  "skipped", file=sys.stderr)
            continue

        if chain == HEAVY_CHAIN:
            records = [("H", mutated), ("L", wt[LIGHT_CHAIN])]
        else:
            records = [("H", wt[HEAVY_CHAIN]), ("L", mutated)]

        # the substitution must actually be present, and nothing else may differ
        parent = wt[chain]
        diffs = [i + 1 for i, (a, b) in enumerate(zip(parent, mutated)) if a != b]
        if diffs != [pos]:
            print(f"{rec}: expected one change at {pos}, found {diffs}", file=sys.stderr)
            continue

        label = entry.get("label") or rec
        short = label.split("__")[-1]
        name = f"{rec}_{short}"
        write_fasta(outdir / f"{name}.fasta", records)
        print(f"  {name:<24} {chain:<13} {entry['substitutions']:<10} "
              f"{region_of(chain, pos):<10} H={len(records[0][1])} L={len(records[1][1])}")
        written += 1

    print(f"\nWrote {written} paired FASTA file(s) to {outdir}/")
    if written:
        print("Each holds 'H' and 'L' variable domains -- upload one per ABodyBuilder job.")


if __name__ == "__main__":
    main()
