#!/usr/bin/env python3
"""Ligand-derived PSSM substitution proposal for the infliximab project.

Drop-in CLI and output schema replacement for worst_anchor_substitution.py.
The seven ligand datasets from Archive 2.zip are embedded as position counts
(20 amino acids, nine positions), deduplicating identical assigned cores.
The 17 known core-to-HLA mappings come from the original NetMHCIIpan 4.0
EL rank <= 1% predictions in full_allele_list_netMHCIIpan_out.txt.

Default: inspect all nine positions and retain up to four candidate positions.
For each position, consider only changes that decrease the motif score for
EVERY available target HLA; among changes within one bit of the best mean
reduction, prefer higher BLOSUM62, then larger reduction, then alphabetically.
Choose positions by mean reduction. Native cysteine is protected by default.
For a core with no available target matrix, the legacy lowest-BLOSUM62
method is used on canonical positions (with a warning on stderr); --missing-policy
skip or error can instead prohibit this fallback. Missing targets are never
imputed from another allele. A partially covered core is scored only against
available matrices, with a warning. The original --anchors option still
restricts positions if explicitly specified.

All matrices use log2((count+1)/(N+20)/(1/20)). Uniform background is an
explicit simplifying assumption; scores are motif compatibility, NOT binding
energies. The one-bit tolerance is a heuristic, not a validated cutoff.
Use --motif-report for per-HLA score changes without modifying original TSVs.
Use --hla to override the embedded core mapping for a new core; optionally
supply additional log2-odds matrices with --motif-matrices (TSV).
Review protected structural residues and re-run NetMHCIIpan on full sequences
for all original 15 HLA molecules. No clinical or structural validation here.

Original usage examples (unchanged):
  python pssm_anchor_substitution.py --core SGGGLVQPG
  python pssm_anchor_substitution.py --cores-file cores.txt --fasta parent.fasta --out-fasta mutants.fasta --map mutants.tsv
  python pssm_anchor_substitution.py --cores-file cores.txt --fasta parent.fasta --split --outdir mutant_fastas --ext txt
"""
import argparse
import csv
import math
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


# Embedded data: unique aligned cores from the seven supplied ligand files.
MOTIF_COUNTS = {'DRB1_0301': [[55, 6, 8, 1, 10, 12, 17, 11, 4, 994, 1024, 7, 96, 208, 19, 10, 16, 29, 44, 630],
               [244, 250, 125, 103, 27, 149, 212, 83, 96, 275, 349, 230, 55, 104, 10, 153, 232, 24, 106, 374],
               [252, 172, 98, 60, 42, 144, 219, 117, 57, 231, 349, 190, 96, 85, 267, 221, 177, 15, 115, 294],
               [93, 7, 47, 2422, 14, 8, 190, 15, 17, 11, 27, 6, 10, 25, 14, 220, 37, 2, 8, 28],
               [177, 177, 194, 185, 19, 128, 278, 256, 60, 102, 209, 220, 51, 77, 402, 235, 198, 19, 75, 139],
               [215, 521, 112, 62, 6, 137, 164, 69, 71, 79, 65, 971, 13, 9, 133, 140, 209, 0, 3, 222],
               [321, 76, 157, 126, 36, 145, 322, 144, 27, 193, 474, 84, 62, 87, 105, 234, 211, 6, 97, 294],
               [204, 270, 112, 79, 17, 192, 216, 183, 127, 153, 435, 276, 35, 87, 217, 170, 143, 13, 102, 170],
               [181, 804, 37, 76, 44, 111, 163, 101, 51, 100, 219, 271, 116, 154, 13, 159, 115, 15, 278, 193]],
 'DRB1_0401': [[17, 7, 7, 7, 5, 6, 10, 14, 11, 661, 785, 18, 99, 1353, 26, 44, 16, 222, 1277, 322],
               [287, 516, 196, 105, 58, 262, 276, 112, 147, 349, 487, 374, 63, 126, 5, 378, 326, 60, 178, 602],
               [417, 293, 205, 110, 84, 210, 220, 206, 61, 390, 452, 284, 133, 144, 330, 479, 269, 49, 194, 377],
               [455, 19, 119, 1180, 34, 400, 893, 53, 34, 250, 168, 13, 208, 90, 9, 392, 327, 29, 21, 213],
               [381, 289, 169, 235, 41, 213, 344, 282, 126, 207, 269, 352, 71, 98, 562, 496, 324, 29, 91, 328],
               [235, 33, 619, 758, 8, 14, 43, 49, 16, 49, 46, 23, 2, 13, 97, 1224, 1225, 1, 6, 446],
               [562, 92, 468, 250, 33, 252, 491, 141, 49, 199, 510, 33, 74, 69, 486, 458, 231, 18, 60, 431],
               [498, 390, 149, 113, 29, 249, 359, 273, 138, 163, 516, 402, 68, 113, 431, 377, 218, 38, 110, 273],
               [918, 23, 115, 91, 70, 549, 225, 766, 96, 101, 103, 54, 76, 42, 86, 1199, 189, 2, 37, 165]],
 'DRB1_0701': [[64, 13, 36, 28, 1, 9, 93, 40, 23, 817, 1193, 11, 166, 1519, 48, 42, 62, 357, 1226, 582],
               [402, 567, 236, 198, 18, 445, 401, 204, 359, 335, 535, 610, 84, 202, 11, 515, 376, 93, 205, 534],
               [530, 357, 219, 150, 36, 392, 273, 276, 244, 325, 541, 487, 146, 359, 555, 471, 241, 74, 354, 300],
               [526, 8, 59, 97, 21, 61, 89, 64, 205, 369, 198, 29, 34, 215, 31, 1352, 1929, 46, 24, 973],
               [456, 538, 231, 238, 29, 347, 292, 353, 282, 283, 385, 721, 79, 131, 371, 547, 475, 23, 124, 425],
               [1307, 54, 296, 213, 21, 122, 177, 784, 25, 54, 12, 29, 23, 135, 389, 1831, 459, 6, 55, 338],
               [339, 336, 456, 120, 6, 432, 509, 281, 335, 397, 649, 413, 96, 233, 179, 332, 335, 84, 238, 560],
               [449, 547, 257, 173, 33, 256, 365, 304, 205, 264, 577, 566, 94, 127, 750, 513, 327, 26, 119, 378],
               [319, 60, 66, 114, 32, 97, 178, 150, 31, 1394, 1527, 45, 156, 360, 58, 131, 216, 36, 158, 1202]],
 'DRB1_1101': [[21, 23, 7, 8, 2, 18, 31, 39, 29, 982, 1355, 30, 218, 1897, 19, 23, 49, 359, 1600, 523],
               [435, 701, 250, 220, 14, 474, 583, 173, 247, 393, 736, 696, 73, 237, 22, 502, 455, 64, 193, 765],
               [797, 393, 222, 153, 82, 340, 500, 325, 202, 551, 701, 373, 175, 221, 431, 555, 303, 55, 325, 529],
               [650, 93, 220, 99, 88, 207, 409, 79, 39, 861, 1765, 232, 480, 248, 22, 368, 388, 43, 69, 873],
               [514, 563, 303, 253, 48, 407, 374, 464, 246, 239, 387, 829, 112, 193, 446, 643, 451, 71, 230, 460],
               [203, 2596, 168, 86, 7, 224, 167, 111, 332, 43, 42, 2402, 10, 28, 173, 229, 298, 6, 4, 104],
               [644, 548, 240, 197, 25, 294, 431, 612, 257, 428, 1085, 561, 145, 183, 356, 344, 262, 48, 209, 364],
               [524, 578, 195, 117, 38, 431, 431, 414, 284, 352, 899, 624, 134, 176, 491, 417, 371, 31, 234, 492],
               [1131, 112, 249, 269, 84, 997, 321, 1286, 224, 146, 200, 115, 90, 64, 184, 1140, 322, 14, 56, 229]],
 'DRB1_1501': [[31, 1, 6, 14, 2, 5, 6, 31, 11, 1282, 987, 3, 142, 217, 20, 3, 16, 5, 11, 700],
               [231, 390, 150, 89, 21, 183, 197, 91, 178, 201, 324, 338, 35, 116, 22, 233, 224, 37, 138, 295],
               [293, 152, 143, 69, 36, 179, 173, 138, 89, 337, 281, 137, 75, 141, 278, 277, 172, 72, 188, 263],
               [55, 60, 30, 42, 25, 49, 28, 48, 67, 229, 282, 31, 61, 743, 38, 59, 40, 192, 1235, 179],
               [186, 295, 189, 137, 20, 170, 170, 145, 149, 146, 189, 383, 36, 75, 445, 270, 183, 31, 92, 182],
               [266, 32, 766, 419, 19, 41, 63, 537, 23, 13, 21, 24, 29, 11, 73, 986, 124, 1, 12, 33],
               [193, 201, 203, 109, 9, 201, 194, 138, 226, 167, 189, 225, 46, 292, 277, 273, 107, 69, 262, 112],
               [198, 225, 127, 115, 5, 164, 173, 196, 115, 160, 441, 291, 70, 114, 386, 217, 168, 9, 113, 206],
               [473, 52, 65, 25, 82, 59, 107, 219, 32, 307, 787, 17, 60, 63, 58, 234, 287, 6, 19, 541]],
 'HLA-DQA10201-DQB10202': [[155, 12, 17, 68, 9, 105, 173, 54, 59, 308, 228, 18, 66, 162, 284, 145, 145, 73, 128, 408],
                           [254,
                            180,
                            117,
                            152,
                            1,
                            129,
                            298,
                            148,
                            119,
                            107,
                            137,
                            179,
                            22,
                            80,
                            20,
                            245,
                            152,
                            25,
                            98,
                            154],
                           [118, 18, 32, 156, 4, 68, 73, 64, 24, 23, 15, 53, 20, 61, 35, 819, 920, 22, 37, 55],
                           [210, 1, 63, 389, 64, 136, 290, 107, 28, 263, 312, 0, 43, 35, 7, 104, 171, 9, 6, 379],
                           [182, 144, 68, 215, 3, 123, 384, 122, 56, 111, 175, 162, 34, 89, 112, 207, 200, 10, 80, 140],
                           [332, 17, 40, 207, 31, 81, 571, 77, 25, 129, 148, 13, 20, 25, 240, 216, 150, 15, 10, 270],
                           [170, 60, 119, 436, 6, 148, 545, 120, 68, 99, 179, 31, 25, 43, 150, 141, 88, 19, 56, 114],
                           [199,
                            208,
                            108,
                            145,
                            4,
                            114,
                            316,
                            110,
                            87,
                            117,
                            245,
                            151,
                            32,
                            51,
                            202,
                            185,
                            110,
                            18,
                            63,
                            152],
                           [175, 56, 39, 183, 13, 109, 346, 78, 27, 268, 422, 50, 137, 167, 32, 113, 87, 56, 71, 188]],
 'HLA-DQA10501-DQB10201': [[33, 2, 2, 5, 0, 16, 34, 6, 11, 227, 106, 4, 30, 121, 143, 23, 28, 39, 96, 201],
                           [112, 105, 64, 63, 1, 51, 134, 34, 56, 62, 59, 47, 10, 34, 4, 109, 71, 10, 33, 68],
                           [132, 13, 31, 68, 12, 36, 39, 118, 27, 20, 25, 30, 12, 39, 35, 235, 173, 5, 37, 40],
                           [66, 5, 26, 238, 21, 44, 141, 41, 7, 124, 134, 1, 24, 15, 1, 37, 46, 8, 1, 147],
                           [66, 57, 52, 74, 0, 61, 176, 58, 36, 40, 57, 82, 15, 41, 43, 98, 64, 11, 31, 65],
                           [146, 13, 19, 97, 7, 37, 267, 33, 8, 36, 56, 10, 8, 20, 102, 59, 71, 0, 7, 131],
                           [73, 11, 54, 216, 1, 56, 274, 49, 24, 37, 93, 14, 5, 12, 52, 46, 26, 2, 24, 58],
                           [101, 77, 33, 85, 1, 43, 100, 54, 24, 71, 120, 53, 18, 32, 100, 53, 46, 14, 33, 69],
                           [61, 36, 21, 85, 3, 47, 110, 35, 6, 117, 205, 5, 69, 76, 6, 43, 46, 29, 54, 73]]}

CORE_HLAS = {'INTVESEDI': ['HLA-DQA10201-DQB10202', 'HLA-DQA10501-DQB10201'],
 'LIKYASESM': ['DRB1_1501'],
 'LQSGNSQES': ['HLA-DQA10103-DQB10301', 'HLA-DQA10401-DQB10301'],
 'QESVTEQDS': ['HLA-DQA10201-DQB10202'],
 'SGGGLVQPG': ['HLA-DQA10103-DQB10301', 'HLA-DQA10401-DQB10301'],
 'TSGGTAALG': ['HLA-DQA10401-DQB10301'],
 'VHTFPAVLQ': ['DRB1_1501'],
 'VYACEVTHQ': ['HLA-DQA10201-DQB10202', 'HLA-DQA10401-DQB10301', 'HLA-DQA10501-DQB10201'],
 'WKVDNALQS': ['DRB1_0301', 'DRB1_0401'],
 'WNSGALTSG': ['HLA-DQA10401-DQB10301'],
 'WVAEIRSKS': ['DRB1_1101'],
 'WYQQRTNGS': ['DRB1_0401'],
 'YACEVTHQG': ['DRB1_0401'],
 'YFPEPVTVS': ['DRB1_0401'],
 'YLQMTDLRT': ['DRB1_0401'],
 'YSLSSTLTL': ['DRB1_0701'],
 'YSLSSVVTV': ['DRB1_0701']}

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


def build_matrices():
    """Pseudocount of one per amino acid, uniform background 0.05."""
    return {hla: [{a: math.log2((n + 1) * 20 / (sum(row) + 20))
                   for a, n in zip(AA, row)} for row in counts]
            for hla, counts in MOTIF_COUNTS.items()}


def read_extra_matrices(path):
    """TSV: allele, position (1-9), and one column per amino acid (log2 odds)."""
    out = {}
    with open(path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if not set(["allele", "position", *AA]).issubset(reader.fieldnames or []):
            raise ValueError("matrix TSV needs allele, position, and all 20 amino acid columns")
        for row in reader:
            hla, pos = row["allele"].strip(), int(row["position"])
            if not hla or not 1 <= pos <= 9:
                raise ValueError("invalid allele or position in matrix TSV")
            matrix = out.setdefault(hla, {})
            if pos in matrix:
                raise ValueError(f"duplicate matrix row {hla} P{pos}")
            vals = {a: float(row[a]) for a in AA}
            if not all(math.isfinite(v) for v in vals.values()):
                raise ValueError(f"non-finite matrix score for {hla} P{pos}")
            matrix[pos] = vals
    if not out or any(set(m) != set(range(1, 10)) for m in out.values()):
        raise ValueError("each extra HLA must have exactly nine matrix rows")
    return {hla: [m[p] for p in range(1, 10)] for hla, m in out.items()}


def score(core, matrix):
    return sum(matrix[i][a] for i, a in enumerate(core))


def process_core(core, anchors, allowed, matrices, targets, tolerance=1.0,
                 max_positions=4, missing_policy="legacy", allow_native_c=False):
    """Return legacy-schema rows, combined core, method and target metadata."""
    core = core.strip().upper()
    if len(core) != 9 or set(core) - set(AA):
        raise ValueError(f"core must contain nine standard amino acids: {core!r}")
    supported = [a for a in targets if a in matrices]
    missing = [a for a in targets if a not in matrices]
    if not targets:
        print(f"WARNING {core}: no embedded HLA mapping; use --hla to specify targets.", file=sys.stderr)
    if missing:
        print(f"WARNING {core}: no motif matrix for {', '.join(missing)}; not scored.", file=sys.stderr)
    if not supported:
        if missing_policy != "legacy":
            raise ValueError(f"{core}: no available target motif matrix ({missing_policy})")
        print(f"WARNING {core}: using LEGACY BLOSUM62 fallback, NOT a motif prediction.", file=sys.stderr)
    method = "PSSM" if supported else "legacy_BLOSUM62"
    positions = anchors if anchors is not None else (range(1, 10) if supported else [1,4,6,9])
    candidates = []
    for p in positions:
        wt = core[p-1]
        if wt == "C" and not allow_native_c:
            continue
        replacements = [a for a in allowed if a != wt]
        if not replacements:
            continue
        if supported:
            options = []
            for a in replacements:
                deltas = [matrices[h][p-1][a] - matrices[h][p-1][wt] for h in supported]
                if all(d < -1e-12 for d in deltas):
                    options.append((sum(deltas)/len(deltas), BLOSUM62[wt,a], a))
            if not options:
                continue
            best_delta = min(o[0] for o in options)
            eligible = [o for o in options if o[0] <= best_delta + tolerance + 1e-12]
            delta, b62, chosen = min(eligible, key=lambda o: (-o[1], o[0], o[2]))
            picks = sorted(o[2] for o in eligible if o[1] == b62 and abs(o[0]-delta) < 1e-12)
        else:
            b62, picks = worst_substitutions(wt, allowed)
            chosen, delta = picks[0], 0.0
        row = {"core":core, "anchor":p, "wt":wt, "chosen":chosen,
               "blosum62":b62, "self_score":BLOSUM62[wt,wt], "ties":"".join(picks),
               "single_mutant_core":core[:p-1]+chosen+core[p:]}
        candidates.append((delta, p, row))
    if anchors is None and supported:
        candidates = sorted(candidates, key=lambda x: (x[0], x[1]))[:max_positions]
    rows = [r for _, _, r in sorted(candidates, key=lambda x:x[1])]
    if not rows:
        raise ValueError(f"{core}: no eligible substitutions; no mutants generated")
    combined = list(core)
    for r in rows:
        combined[r["anchor"]-1] = r["chosen"]
    return rows, "".join(combined), method, supported, missing


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

    ap.add_argument("--anchors", default=None,
                    help="Restrict 1-based positions; default: automatically assess all nine")
    ap.add_argument("--max-positions", type=int, default=4,
                    help="Maximum positions in automatic PSSM mode (default 4)")
    ap.add_argument("--motif-tolerance", type=float, default=1.0,
                    help="Allow this many bits above best mean delta before preferring BLOSUM62")
    ap.add_argument("--hla", help="Override HLA targets for all input cores: comma-separated exact NetMHCIIpan names")
    ap.add_argument("--motif-matrices", help="Optional extra TSV matrices: allele, position, 20 amino acid columns (log2 odds)")
    ap.add_argument("--missing-policy", choices=["legacy", "skip", "error"], default="legacy",
                    help="Action for cores with no supported HLA matrix (default legacy with warning)")
    ap.add_argument("--allow-native-c", action="store_true", help="Allow replacing existing cysteines (protected by default)")
    ap.add_argument("--motif-report", help="Optional additional TSV showing score changes per target HLA")
    ap.add_argument("--exclude", default="C",
                    help="Residues never INTRODUCED as replacements (default: C, to avoid "
                         "creating free thiols). Native C is protected unless --allow-native-c is set.  Use 'CP' to also block proline, '' for all 20.")

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

    anchors = None
    if args.anchors is not None:
        try:
            anchors = sorted(set(int(x) for x in args.anchors.split(",") if x.strip()))
        except ValueError:
            ap.error("--anchors must be comma-separated integers, e.g. 1,4,6,9")
        if not anchors or any(p < 1 or p > 9 for p in anchors):
            ap.error("--anchors values must be between 1 and 9")
    if not 1 <= args.max_positions <= 9:
        ap.error("--max-positions must be 1-9")
    if not math.isfinite(args.motif_tolerance) or args.motif_tolerance < 0:
        ap.error("--motif-tolerance must be finite and nonnegative")
    matrices = build_matrices()
    if args.motif_matrices:
        try:
            matrices.update(read_extra_matrices(args.motif_matrices))
        except (ValueError, OSError, TypeError, KeyError) as e:
            ap.error(str(e))
    override = sorted(set(x.strip() for x in args.hla.split(",") if x.strip())) if args.hla else None
    if args.hla and not override:
        ap.error("--hla must contain at least one allele")
    report = []

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
            targets = override if override is not None else CORE_HLAS.get(core, [])
            rows, combined_core, method, supported, missing = process_core(
                core, anchors, allowed, matrices, targets, args.motif_tolerance,
                args.max_positions, args.missing_policy, args.allow_native_c)
        except ValueError as exc:
            if args.missing_policy == "error" and not any(a in matrices for a in targets):
                ap.error(str(exc))
            problems.append(str(exc))
            continue
        anchor_rows.extend(rows)
        print(f"INFO {core}: {method}; HLA targets: {', '.join(targets) or 'unknown'}", file=sys.stderr)
        if args.motif_report:
            variants = [(r["single_mutant_core"], f"P{r['anchor']}") for r in rows]
            variants.append((combined_core, "combined"))
            for variant, position in variants:
                for hla in targets or [""]:
                    old = score(core, matrices[hla]) if hla in supported else None
                    new = score(variant, matrices[hla]) if hla in supported else None
                    report.append({"core":core,"mutant_core":variant,"position":position,
                                   "hla":hla,"method":method,"matrix_available":hla in supported,
                                   "wt_pssm":round(old,6) if old is not None else "",
                                   "mutant_pssm":round(new,6) if new is not None else "",
                                   "delta_pssm":round(new-old,6) if old is not None else ""})

        print(f"\nCore {core}   (never introduced: {''.join(sorted(excluded)) or 'none'})")
        print(f"  {'pos':>4}  {'WT':^3} {'->':^3} {'B62':>4} {'self':>4}  "
              f"{'ties':<8}  mutant core")
        for r in rows:
            print(f"  {'P' + str(r['anchor']):>4}  {r['wt']:^3} {r['chosen']:^3} "
                  f"{r['blosum62']:>4} {r['self_score']:>4}  {r['ties']:<8}  "
                  f"{r['single_mutant_core']}")
        print(f"  combined ({len(rows)} anchors): {combined_core}")

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

    if args.motif_report:
        with open(args.motif_report, "w", newline="") as fh:
            fields = ["core","mutant_core","position","hla","method","matrix_available","wt_pssm","mutant_pssm","delta_pssm"]
            writer = csv.DictWriter(fh, delimiter="\t", fieldnames=fields)
            writer.writeheader()
            writer.writerows(report)

    if problems:
        print(f"\n{len(problems)} problem(s):", file=sys.stderr)
        for p in problems:
            print(f"  - {p}", file=sys.stderr)


if __name__ == "__main__":
    main()
