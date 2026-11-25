#!/usr/bin/env python3
"""Pick the most dissimilar (lowest-BLOSUM62) substitution at MHC-II anchor
positions, and optionally write the mutated full-length chains as FASTA.

For a 9-mer MHC class II binding core, the P1/P4/P6/P9 pockets carry most of the
binding energy. This script reports, for each anchor, the amino acid with the
LOWEST BLOSUM62 score against the wild-type residue -- i.e. the most chemically
dissimilar replacement -- builds the combined quadruple mutant, and (with
--fasta) locates each core in the parent chain and writes new FASTA records
carrying the substitutions.

Caveat worth stating in any write-up: BLOSUM62 scores evolutionary
interchangeability, not MHC binding. A low-BLOSUM substitution is a reasonable
first guess at "disrupt this anchor", but it is not the same as "lowest predicted
binding", and it can create a new register rather than removing the epitope.
Re-run NetMHCIIpan on the output FASTA to confirm the rank actually moved.

Examples
--------
  # cores only, no FASTA
  python worst_anchor_substitution.py --core SGGGLVQPG

  # mutate the parent chains, one multi-record FASTA for a single NetMHCIIpan run
  python worst_anchor_substitution.py --cores-file cores.txt \
      --fasta rcsb_pdb_4G3Y.fasta --out-fasta mutants.fasta --map mutants.tsv

  # one .txt FASTA file per mutation instead
  python worst_anchor_substitution.py --cores-file cores.txt \
      --fasta rcsb_pdb_4G3Y.fasta --split --outdir mutant_fastas --ext txt
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


# --------------------------------------------------------------------------- #
# FASTA handling
# --------------------------------------------------------------------------- #

def read_fasta(path):
    """Return {identifier: sequence}. Identifier is the first whitespace token."""
    sequences, name, parts = {}, None, []
    for raw in Path(path).read_text().splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name is not None:
                sequences[name] = "".join(parts).upper()
            name, parts = line[1:].split()[0], []
            if name in sequences:
                raise ValueError(f"Duplicate FASTA identifier: {name}")
        else:
            if name is None:
                raise ValueError("FASTA must start with a '>' header")
            parts.append(line)
    if name is not None:
        sequences[name] = "".join(parts).upper()
    if not sequences:
        raise ValueError("Empty FASTA")
    return sequences


def write_fasta(path, records, width=60):
    """records: list of (header, sequence)."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as fh:
        for header, seq in records:
            fh.write(f">{header}\n")
            for i in range(0, len(seq), width):
                fh.write(seq[i:i + width] + "\n")


def locate_core(sequences, core, chain_filter=None):
    """Find `core` in exactly one chain. Returns (chain, zero-based start).

    Raises ValueError if the core is absent, appears in more than one chain,
    or appears more than once within a chain.
    """
    hits = []
    for chain, seq in sequences.items():
        if chain_filter and chain != chain_filter:
            continue
        first = seq.find(core)
        if first == -1:
            continue
        if seq.find(core, first + 1) != -1:
            raise ValueError(f"core {core} occurs more than once in chain {chain}")
        hits.append((chain, first))
    if not hits:
        where = f"chain {chain_filter}" if chain_filter else "any chain"
        raise ValueError(f"core {core} not found in {where}")
    if len(hits) > 1:
        raise ValueError(
            f"core {core} found in several chains ({', '.join(c for c, _ in hits)}); "
            "use --chain to disambiguate")
    return hits[0]


# --------------------------------------------------------------------------- #
# Substitution logic
# --------------------------------------------------------------------------- #

def worst_substitutions(residue, allowed):
    """Return (min_score, [amino acids achieving it]) for `residue`."""
    scored = [(BLOSUM62[(residue, a)], a) for a in allowed if a != residue]
    if not scored:
        raise ValueError(
            f"no replacement available for {residue} with the current --exclude set")
    low = min(s for s, _ in scored)
    return low, sorted(a for s, a in scored if s == low)


def process_core(core, anchors, allowed):
    """Return (per-anchor rows, combined mutant core)."""
    core = core.strip().upper()
    if len(core) != 9:
        raise ValueError(f"core must be 9 residues, got {len(core)}: {core!r}")
    bad = sorted(set(core) - set(AA))
    if bad:
        raise ValueError(f"core {core!r} contains non-standard residue(s): {','.join(bad)}")

    rows, combined = [], list(core)
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
            "single_mutant_core": core[:p - 1] + chosen + core[p:],
        })
    return rows, "".join(combined)


# --------------------------------------------------------------------------- #

def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--core", help="A single 9-mer binding core, e.g. SGGGLVQPG")
    src.add_argument("--cores-file",
                     help="Whitespace- or tab-separated file of binder rows. By default "
                          "the first field of each line is taken as the core; use "
                          "--core-column for a different layout. A header line is skipped "
                          "automatically. '#' comments and blank lines are ignored.")

    ap.add_argument("--core-column", type=int, default=1,
                    help="1-based column holding the 9-mer core in --cores-file "
                         "(default: 1)")
    ap.add_argument("--chain-column", type=int,
                    help="1-based column holding the FASTA identifier for each row. "
                         "If given, each core is mutated in that chain instead of being "
                         "searched for across the whole FASTA.")

    ap.add_argument("--anchors", default="1,4,6,9",
                    help="1-based anchor positions within the core (default: 1,4,6,9)")
    ap.add_argument("--exclude", default="C",
                    help="Residues never INTRODUCED as replacements (default: C, to avoid "
                         "creating free thiols). Note this does not stop a wild-type C "
                         "from being mutated. Use 'CP' to also block proline, '' for all 20.")

    ap.add_argument("--fasta",
                    help="Parent FASTA (e.g. rcsb_pdb_4G3Y.fasta). Enables FASTA output.")
    ap.add_argument("--chain",
                    help="Restrict core lookup to this FASTA identifier")
    ap.add_argument("--mutants", choices=["single", "combined", "both"], default="both",
                    help="Which mutated chains to emit (default: both)")
    ap.add_argument("--id-style", choices=["short", "descriptive"], default="short",
                    help="FASTA record naming. 'short' (default) writes WT1/WT2 and "
                         "M01, M02, ... because NetMHCIIpan truncates the Identity "
                         "field at 15 characters, which silently merges descriptive "
                         "names into one record. 'descriptive' writes "
                         "<chain>__<substitutions>; only safe if the FASTA is not "
                         "going through NetMHCIIpan.")
    ap.add_argument("--id-prefix", default="M",
                    help="Prefix for short mutant ids (default: M)")
    ap.add_argument("--no-wt", action="store_true",
                    help="Omit the unmutated parent chains from the FASTA output "
                         "(they are included by default as a NetMHCIIpan baseline)")

    ap.add_argument("--out-fasta", default="mutants.fasta",
                    help="Multi-record output FASTA (default: mutants.fasta)")
    ap.add_argument("--split", action="store_true",
                    help="Write one FASTA file per mutant instead of one combined file")
    ap.add_argument("--outdir", default="mutant_fastas",
                    help="Directory for --split output (default: mutant_fastas)")
    ap.add_argument("--ext", default="fasta",
                    help="File extension for --split output, e.g. 'txt' (default: fasta)")

    ap.add_argument("--map", dest="mapfile",
                    help="Write a TSV mapping every record back to core/anchor/position")
    ap.add_argument("--tsv", help="Write the per-anchor BLOSUM62 table to this TSV")

    args = ap.parse_args()

    try:
        anchors = [int(x) for x in args.anchors.split(",") if x.strip()]
    except ValueError:
        ap.error("--anchors must be comma-separated integers, e.g. 1,4,6,9")
    if not anchors or any(p < 1 or p > 9 for p in anchors):
        ap.error("--anchors values must be between 1 and 9")
    anchors = sorted(dict.fromkeys(anchors))

    excluded = set(args.exclude.upper().replace(",", "").replace(" ", ""))
    if excluded - set(AA):
        ap.error("--exclude must contain only standard one-letter amino acid codes")
    allowed = [a for a in AA if a not in excluded]
    if not allowed:
        ap.error("--exclude leaves no replacement residues")

    if args.core_column < 1 or (args.chain_column is not None and args.chain_column < 1):
        ap.error("--core-column and --chain-column are 1-based")

    def looks_like_core(text):
        return len(text) == 9 and not (set(text.upper()) - set(AA))

    pairs = []  # (core, chain_or_None), in file order
    if args.core:
        pairs.append((args.core.upper(), args.chain))
    else:
        raw_lines = Path(args.cores_file).read_text().splitlines()
        for line_no, line in enumerate(raw_lines, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) < args.core_column:
                print(f"line {line_no}: only {len(fields)} column(s), need "
                      f"{args.core_column}; skipped", file=sys.stderr)
                continue
            core = fields[args.core_column - 1].upper()
            if not pairs and not looks_like_core(core):
                print(f"Skipping header line: {line}", file=sys.stderr)
                continue
            chain = args.chain
            if args.chain_column is not None:
                if len(fields) < args.chain_column:
                    print(f"line {line_no}: no column {args.chain_column}; skipped",
                          file=sys.stderr)
                    continue
                chain = fields[args.chain_column - 1]
            pairs.append((core, chain))
        if not pairs:
            ap.error(f"no cores found in {args.cores_file} at column {args.core_column}")

    seen, cores = set(), []
    for core, chain in pairs:
        if (core, chain) not in seen:
            seen.add((core, chain))
            cores.append((core, chain))

    sequences = None
    if args.fasta:
        try:
            sequences = read_fasta(args.fasta)
        except (OSError, ValueError) as exc:
            ap.error(f"could not read --fasta: {exc}")
        unknown = sorted({c for _, c in cores if c and c not in sequences})
        if unknown:
            ap.error(f"chain identifier(s) not in the FASTA: {', '.join(unknown)} "
                     f"(FASTA has: {', '.join(sequences)})")

    anchor_rows, fasta_records, map_rows, problems = [], [], [], []
    counter = {"n": 0}
    if sequences is None:
        wt_ids = {}
    elif args.id_style == "short":
        wt_ids = {c: f"WT{i}" for i, c in enumerate(sequences, 1)}
    else:
        wt_ids = {c: f"WT__{c}" for c in sequences}

    for core, core_chain in cores:
        try:
            rows, combined_core = process_core(core, anchors, allowed)
        except ValueError as exc:
            problems.append(str(exc))
            continue
        anchor_rows.extend(rows)

        print(f"\nCore {core}   (never introduced: {''.join(sorted(excluded)) or 'none'})")
        print(f"  {'pos':>4}  {'WT':^3} {'->':^3} {'B62':>4} {'self':>4}  "
              f"{'ties':<8}  mutant core")
        for r in rows:
            print(f"  {'P' + str(r['anchor']):>4}  {r['wt']:^3} {r['chosen']:^3} "
                  f"{r['blosum62']:>4} {r['self_score']:>4}  {r['ties']:<8}  "
                  f"{r['single_mutant_core']}")
        print(f"  combined ({len(anchors)} anchors): {combined_core}")

        if sequences is None:
            continue

        try:
            chain, start = locate_core(sequences, core, core_chain)
        except ValueError as exc:
            problems.append(str(exc))
            print(f"  !! no FASTA output: {exc}")
            continue

        parent = sequences[chain]
        print(f"  located in {chain} at residues {start + 1}-{start + 9}")

        def emit(mut_core, label, detail):
            seq = parent[:start] + mut_core + parent[start + 9:]
            descriptive = f"{chain}__{label}"
            if args.id_style == "short":
                counter["n"] += 1
                header = f"{args.id_prefix}{counter['n']:02d}"
            else:
                header = descriptive
            fasta_records.append((header, seq))
            map_rows.append({
                "record_id": header,
                "label": descriptive,
                "wt_record_id": wt_ids[chain],
                "chain": chain,
                "core_wt": core,
                "core_mutant": mut_core,
                "core_start_1based": start + 1,
                "substitutions": detail,
                "n_substitutions": detail.count(",") + 1,
                "sequence_length": len(seq),
            })

        if args.mutants in ("single", "both"):
            for r in rows:
                pos = start + r["anchor"]          # 1-based position in the chain
                label = f"{r['wt']}{pos}{r['chosen']}"
                emit(r["single_mutant_core"], label, label)

        if args.mutants in ("combined", "both"):
            subs = [f"{r['wt']}{start + r['anchor']}{r['chosen']}" for r in rows]
            emit(combined_core, "_".join(subs), ",".join(subs))

    # ----------------------------------------------------------------- output
    if sequences is not None:
        if not args.no_wt:
            wt_records = [(wt_ids[c], s) for c, s in sequences.items()]
            fasta_records = wt_records + fasta_records

        if not fasta_records:
            print("\nNo FASTA records produced.", file=sys.stderr)
        elif args.split:
            outdir = Path(args.outdir)
            ext = args.ext.lstrip(".")
            labels = {r["record_id"]: r["label"] for r in map_rows}
            for header, seq in fasta_records:
                stem = header
                if header in labels and labels[header] != header:
                    stem = f"{header}_{labels[header]}"
                write_fasta(outdir / f"{stem}.{ext}", [(header, seq)])
            print(f"\nWrote {len(fasta_records)} files to {outdir}/ (.{ext})")
        else:
            out = Path(args.out_fasta)
            write_fasta(out, fasta_records)
            print(f"\nWrote {len(fasta_records)} records to {out}")

        if args.mapfile and map_rows:
            with open(args.mapfile, "w", newline="") as fh:
                w = csv.DictWriter(fh, delimiter="\t", fieldnames=list(map_rows[0].keys()))
                w.writeheader()
                w.writerows(map_rows)
            print(f"Wrote mapping for {len(map_rows)} mutants to {args.mapfile}")

    if args.tsv and anchor_rows:
        with open(args.tsv, "w", newline="") as fh:
            w = csv.DictWriter(fh, delimiter="\t", fieldnames=list(anchor_rows[0].keys()))
            w.writeheader()
            w.writerows(anchor_rows)
        print(f"Wrote {len(anchor_rows)} anchor rows to {args.tsv}")

    if problems:
        print(f"\n{len(problems)} problem(s):", file=sys.stderr)
        for p in problems:
            print(f"  - {p}", file=sys.stderr)


if __name__ == "__main__":
    main()
