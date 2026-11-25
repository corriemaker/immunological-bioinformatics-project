#!/usr/bin/env python3
"""Compare NetMHCIIpan predictions for mutant chains against their wild type.

Reads one or more NetMHCIIpan-4.x output files containing BOTH the wild-type
records (identity 'WT__<chain>') and the mutant records (identity
'<chain>__<substitutions>'), plus the mapping TSV written by
worst_anchor_substitution.py --map.

For each mutant it answers the two questions that matter at this step:

  1. Did the targeted epitope go away?  For every allele that bound somewhere in
     the mutated window in the wild type, what is the best %Rank now?
  2. Did a new epitope appear?  A substitution can shift the binding register and
     create a binder one position over. Any allele that was a non-binder in the
     window and is a binder after mutation is reported as NEW.

Only 15-mers that actually contain a mutated residue can change, so the
comparison is restricted to that window and joined wild-type-to-mutant on
(position, allele) -- an exact correspondence, since a substitution preserves
both length and numbering.

Verdicts
--------
  CLEARED   every target allele pushed above threshold, no new binders
  PARTIAL   some target alleles cleared, others still bind
  FAILED    no target allele cleared
  NEW_ONLY  nothing was bound before, but the mutation created a binder
  NEW       cleared or partial, BUT a new allele now binds (register shift)
  NO_TARGET nothing bound in this window in the wild type and nothing does now

Example
-------
  python compare_mutant_predictions.py \
      --predictions netMHCIIpan_mutants_out.txt \
      --map mutants.tsv \
      --threshold 2.0 \
      --out step6_summary.tsv --details step6_by_allele.tsv
"""
import argparse
import csv
import re
import sys
from collections import defaultdict
from pathlib import Path

PEPTIDE_WINDOW = 15  # NetMHCIIpan default peptide length
SUB_RE = re.compile(r"^([A-Z])(\d+)([A-Z])$")


# --------------------------------------------------------------------------- #
# Parsing
# --------------------------------------------------------------------------- #

def parse_netmhciipan(paths):
    """Yield dict rows from NetMHCIIpan output file(s).

    Tolerates CRLF line endings, '# Allele:' banners, dashed rules, the column
    header, and the trailing summary block. Works with or without BA columns.
    """
    rows = []
    for path in paths:
        text = Path(path).read_text(errors="replace")
        for raw in text.splitlines():
            line = raw.replace("\r", "").rstrip()
            if not line or line.lstrip().startswith(("#", "-", "Pos")):
                continue
            f = line.split()
            if not f or not f[0].isdigit():
                continue
            if len(f) < 9:
                continue
            row = {
                "pos": int(f[0]),
                "mhc": f[1],
                "peptide": f[2],
                "core": f[4],
                "identity": f[6],
                "rank_el": _tofloat(f[8]),
            }
            # BA block present when there are >= 13 fields
            row["rank_ba"] = _tofloat(f[12]) if len(f) >= 13 else None
            rows.append(row)
    return rows


def _tofloat(text):
    try:
        return float(text)
    except (TypeError, ValueError):
        return None


def parse_substitutions(text):
    """'W148D,D151L' -> [(148, 'W', 'D'), (151, 'D', 'L')]"""
    out = []
    for token in str(text).replace(";", ",").split(","):
        token = token.strip()
        if not token:
            continue
        m = SUB_RE.match(token)
        if not m:
            raise ValueError(f"cannot parse substitution {token!r}")
        wt, pos, mut = m.group(1), int(m.group(2)), m.group(3)
        out.append((pos, wt, mut))
    if not out:
        raise ValueError("no substitutions found")
    return out


def resolve_identities(wanted_ids, available):
    """Map each wanted record id to an identity in the prediction output.

    NetMHCIIpan truncates long identities, so an exact match is tried first and
    only the leftovers fall back to prefix matching. A prefix candidate is only
    considered if it is not itself the exact identity of some other record, and
    it is rejected if it could stand for more than one record -- truncation that
    collapses two names into one is unrecoverable, and guessing would silently
    attribute predictions to the wrong mutant.

    Returns (resolved: dict, problems: list[str]).
    """
    resolved, problems = {}, []
    exact = {w for w in wanted_ids if w in available}
    for w in exact:
        resolved[w] = w

    # identities not already claimed by an exact match are truncation candidates
    free = set(available) - exact
    unmatched = [w for w in wanted_ids if w not in resolved]

    claims = defaultdict(list)
    for w in unmatched:
        for a in free:
            if a and w.startswith(a) and len(a) < len(w):
                claims[a].append(w)

    for w in unmatched:
        hits = [a for a in claims if w in claims[a]]
        if not hits:
            continue
        if len(hits) > 1:
            problems.append(
                f"{w}: several truncated identities could match "
                f"({', '.join(sorted(hits))}); rename the FASTA headers")
            continue
        a = hits[0]
        if len(claims[a]) > 1:
            problems.append(
                f"identity {a!r} is ambiguous -- it could be any of "
                f"{', '.join(sorted(claims[a]))}; NetMHCIIpan truncated the "
                "headers, so shorten them and re-run the prediction")
            continue
        resolved[w] = a
    return resolved, problems


# --------------------------------------------------------------------------- #

def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--predictions", nargs="+", required=True,
                    help="NetMHCIIpan output file(s) containing WT and mutant records")
    ap.add_argument("--map", dest="mapfile", required=True,
                    help="mutants.tsv written by worst_anchor_substitution.py --map")
    ap.add_argument("--threshold", type=float, default=2.0,
                    help="%%Rank at or below which a peptide counts as a binder "
                         "(default: 2.0, the NetMHCIIpan-4.1 weak-binder cutoff)")
    ap.add_argument("--metric", choices=["rank_el", "rank_ba"], default="rank_el",
                    help="Which %%Rank column to use (default: rank_el)")
    ap.add_argument("--peptide-length", type=int, default=PEPTIDE_WINDOW,
                    help="Peptide length used in the prediction (default: 15)")
    ap.add_argument("--out", help="Write the per-mutant summary to this TSV")
    ap.add_argument("--details", help="Write the per-mutant-per-allele table to this TSV")
    args = ap.parse_args()

    metric = args.metric
    thr = args.threshold
    plen = args.peptide_length

    rows = parse_netmhciipan(args.predictions)
    if not rows:
        sys.exit("No prediction rows parsed -- is this a NetMHCIIpan output file?")

    # index: identity -> (pos, mhc) -> row
    by_id = defaultdict(dict)
    for r in rows:
        if r[metric] is None:
            continue
        by_id[r["identity"]][(r["pos"], r["mhc"])] = r
    identities = set(by_id)
    print(f"Parsed {len(rows)} rows across {len(identities)} records "
          f"from {len(args.predictions)} file(s)")

    with open(args.mapfile, newline="") as fh:
        mapping = list(csv.DictReader(fh, delimiter="\t"))
    if not mapping:
        sys.exit(f"{args.mapfile} is empty")

    summary, details, warnings = [], [], []

    def wt_id_for(entry):
        """Prefer the explicit wt_record_id column; fall back to the old scheme."""
        return entry.get("wt_record_id") or f"WT__{entry['chain']}"

    wanted = [e["record_id"] for e in mapping]
    wanted += sorted({wt_id_for(e) for e in mapping})
    resolved, id_problems = resolve_identities(wanted, identities)
    warnings.extend(id_problems)

    for entry in mapping:
        rec_id = entry["record_id"]
        chain = entry["chain"]

        label = entry.get("label") or rec_id
        wt_id = wt_id_for(entry)
        mut_key = resolved.get(rec_id)
        wt_key = resolved.get(wt_id)
        if mut_key is None:
            warnings.append(f"{label}: no record {rec_id!r} in the prediction output")
            continue
        if wt_key is None:
            warnings.append(f"{label}: no record {wt_id!r} in the prediction output")
            continue

        try:
            subs = parse_substitutions(entry["substitutions"])
        except ValueError as exc:
            warnings.append(f"{rec_id}: {exc}")
            continue

        # 15-mers containing at least one mutated residue
        window = set()
        for pos, _, _ in subs:
            window.update(range(max(1, pos - plen + 1), pos + 1))

        wt_rows, mut_rows = by_id[wt_key], by_id[mut_key]

        best_wt, best_mut = {}, {}
        for (pos, mhc), r in wt_rows.items():
            if pos in window:
                cur = best_wt.get(mhc)
                if cur is None or r[metric] < cur[metric]:
                    best_wt[mhc] = r
        for (pos, mhc), r in mut_rows.items():
            if pos in window:
                cur = best_mut.get(mhc)
                if cur is None or r[metric] < cur[metric]:
                    best_mut[mhc] = r

        target_alleles = sorted(m for m, r in best_wt.items() if r[metric] <= thr)
        new_alleles = sorted(
            m for m, r in best_mut.items()
            if r[metric] <= thr and (m not in best_wt or best_wt[m][metric] > thr))

        cleared, retained = [], []
        for mhc in target_alleles:
            mr = best_mut.get(mhc)
            if mr is None or mr[metric] > thr:
                cleared.append(mhc)
            else:
                retained.append(mhc)

        if not target_alleles:
            verdict = "NEW_ONLY" if new_alleles else "NO_TARGET"
        elif retained and not cleared:
            verdict = "FAILED"
        elif retained:
            verdict = "PARTIAL"
        else:
            verdict = "CLEARED"
        if new_alleles and verdict in ("CLEARED", "PARTIAL", "FAILED"):
            verdict = "NEW"

        for mhc in sorted(set(target_alleles) | set(new_alleles)):
            wr, mr = best_wt.get(mhc), best_mut.get(mhc)
            if mhc in new_alleles and mhc not in target_alleles:
                status = "NEW"
            elif mhc in cleared:
                status = "cleared"
            else:
                status = "retained"
            details.append({
                "record_id": rec_id,
                "label": label,
                "chain": chain,
                "allele": mhc,
                "status": status,
                "wt_rank": f"{wr[metric]:.2f}" if wr else "",
                "mut_rank": f"{mr[metric]:.2f}" if mr else "",
                "wt_peptide": wr["peptide"] if wr else "",
                "mut_peptide": mr["peptide"] if mr else "",
                "wt_core": wr["core"] if wr else "",
                "mut_core": mr["core"] if mr else "",
            })

        wt_best = min((r[metric] for r in best_wt.values()), default=None)
        mut_best = min((r[metric] for r in best_mut.values()), default=None)
        summary.append({
            "record_id": rec_id,
            "label": label,
            "chain": chain,
            "substitutions": entry["substitutions"],
            "n_substitutions": entry.get("n_substitutions", ""),
            "core_wt": entry.get("core_wt", ""),
            "core_mutant": entry.get("core_mutant", ""),
            "verdict": verdict,
            "n_target_alleles": len(target_alleles),
            "n_cleared": len(cleared),
            "n_retained": len(retained),
            "n_new": len(new_alleles),
            "best_wt_rank": f"{wt_best:.2f}" if wt_best is not None else "",
            "best_mut_rank": f"{mut_best:.2f}" if mut_best is not None else "",
            "target_alleles": ",".join(target_alleles),
            "retained_alleles": ",".join(retained),
            "new_alleles": ",".join(new_alleles),
        })

    # ------------------------------------------------------------- reporting
    order = ["CLEARED", "PARTIAL", "NEW", "FAILED", "NEW_ONLY", "NO_TARGET"]
    counts = defaultdict(int)
    for s in summary:
        counts[s["verdict"]] += 1

    print(f"\nBinder threshold: %Rank_{'EL' if metric == 'rank_el' else 'BA'} <= {thr}")
    print(f"Mutants evaluated: {len(summary)}\n")
    for v in order:
        if counts[v]:
            print(f"  {v:<10} {counts[v]:>3}")

    interesting = [s for s in summary if s["verdict"] in ("CLEARED", "NEW", "PARTIAL")]
    if interesting:
        print(f"\n{'record':<40} {'verdict':<9} {'WT':>6} {'mut':>6}  retained / new")
        for s in sorted(interesting, key=lambda x: order.index(x["verdict"])):
            note = s["retained_alleles"]
            if s["new_alleles"]:
                note = (note + " / " if note else "") + "NEW:" + s["new_alleles"]
            print(f"  {s['label'][:38]:<38} {s['verdict']:<9} "
                  f"{s['best_wt_rank']:>6} {s['best_mut_rank']:>6}  {note}")

    if args.out and summary:
        with open(args.out, "w", newline="") as fh:
            w = csv.DictWriter(fh, delimiter="\t", fieldnames=list(summary[0].keys()))
            w.writeheader()
            w.writerows(summary)
        print(f"\nWrote summary for {len(summary)} mutants to {args.out}")

    if args.details and details:
        with open(args.details, "w", newline="") as fh:
            w = csv.DictWriter(fh, delimiter="\t", fieldnames=list(details[0].keys()))
            w.writeheader()
            w.writerows(details)
        print(f"Wrote {len(details)} allele rows to {args.details}")

    if warnings:
        print(f"\n{len(warnings)} warning(s):", file=sys.stderr)
        for w_ in warnings:
            print(f"  - {w_}", file=sys.stderr)


if __name__ == "__main__":
    main()
