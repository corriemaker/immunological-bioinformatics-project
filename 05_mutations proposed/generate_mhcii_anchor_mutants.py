#!/usr/bin/env python3
"""Generate full-chain single-substitution FASTA variants at P1/P4/P6/P9 of MHC-II cores.

Input binder rows: PEPTIDE HLA IDENTITY CORE (whitespace separated, no header).
FASTA sequence headers must have IDs matching the IDENTITY column, e.g. 4G3Y_1_Chain.
"""
import argparse
import csv
import re
from pathlib import Path

AA = "ACDEFGHIKLMNPQRSTVWY"
ANCHORS = (1, 4, 6, 9)


def read_fasta(path):
    sequences = {}
    name, parts = None, []
    for raw in Path(path).read_text().splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith('>'):
            if name is not None:
                sequences[name] = ''.join(parts).upper()
            name, parts = line[1:].split()[0], []
            if name in sequences:
                raise ValueError(f'Duplicate FASTA identifier: {name}')
        else:
            if name is None:
                raise ValueError('FASTA must start with a > header')
            parts.append(line)
    if name is not None:
        sequences[name] = ''.join(parts).upper()
    if not sequences:
        raise ValueError('Empty FASTA')
    return sequences


def unique_position(sequence, fragment):
    first = sequence.find(fragment)
    if first == -1 or sequence.find(fragment, first + 1) != -1:
        return None
    return first


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--fasta', required=True, help='FASTA with full heavy/light chain sequences')
    ap.add_argument('--binders', required=True, help='Four-column PopCover binder file')
    ap.add_argument('--out', default='anchor_mutants', help='Output filename prefix')
    ap.add_argument('--replacements', default=AA,
                    help='Allowed replacement amino acids (default: all 20; original is excluded)')
    args = ap.parse_args()
    replacements = ''.join(dict.fromkeys(args.replacements.upper().replace(',', '').replace(' ', '')))
    if not replacements or any(a not in AA for a in replacements):
        ap.error('--replacements must contain standard one-letter amino acid symbols')
    sequences = read_fasta(args.fasta)
    variants = {}  # (chain, zero_based_position, new_amino_acid) -> supporting prediction records
    skipped = []
    binders = 0
    for line_no, raw in enumerate(Path(args.binders).read_text().splitlines(), 1):
        if not raw.strip() or raw.lstrip().startswith('#'):
            continue
        fields = raw.split()
        if len(fields) != 4:
            skipped.append((line_no, 'Expected 4 columns'))
            continue
        peptide, hla, chain, core = fields
        peptide, core = peptide.upper(), core.upper()
        binders += 1
        if chain not in sequences:
            skipped.append((line_no, f'FASTA identifier {chain!r} missing'))
            continue
        if len(core) != 9 or len(peptide) < 9:
            skipped.append((line_no, 'Invalid peptide/core length'))
            continue
        core_offset = unique_position(peptide, core)
        peptide_start = unique_position(sequences[chain], peptide)
        if core_offset is None or peptide_start is None:
            skipped.append((line_no, 'Peptide or core absent/ambiguous; needs manual sequence mapping'))
            continue
        for anchor in ANCHORS:
            idx = peptide_start + core_offset + anchor - 1
            original = sequences[chain][idx]
            if original != core[anchor - 1]:
                skipped.append((line_no, 'Sequence/core mismatch'))
                continue
            for mutant in replacements:
                if mutant == original:
                    continue
                key = (chain, idx, mutant)
                variants.setdefault(key, []).append((hla, peptide, core, anchor))

    outpath = Path(args.out)
    outpath.parent.mkdir(parents=True, exist_ok=True)
    fasta_path = Path(str(outpath) + '.fasta')
    tsv_path = Path(str(outpath) + '.tsv')
    skips_path = Path(str(outpath) + '_skipped.tsv')
    with fasta_path.open('w') as fo, tsv_path.open('w', newline='') as co:
        writer = csv.writer(co, delimiter='\t')
        writer.writerow(['variant_id', 'chain', 'position_1based', 'original_AA',
                         'replacement_AA', 'HLA_alleles', 'anchor_positions',
                         'source_peptides'])
        for chain, seq in sequences.items():
            fo.write(f'>WT__{chain}\n{seq}\n')
        for (chain, idx, mutant), records in sorted(variants.items()):
            seq = sequences[chain]
            original = seq[idx]
            ident = f'{chain}__{original}{idx+1}{mutant}'
            changed = seq[:idx] + mutant + seq[idx+1:]
            fo.write(f'>{ident}\n{changed}\n')
            writer.writerow([ident, chain, idx+1, original, mutant,
                             ','.join(sorted({r[0] for r in records})),
                             ','.join(map(str, sorted({r[3] for r in records}))),
                             ','.join(sorted({r[1] for r in records}))])
    with skips_path.open('w', newline='') as fo:
        writer = csv.writer(fo, delimiter='\t')
        writer.writerow(['input_line', 'reason'])
        writer.writerows(skipped)
    print(f'Read {binders} binder rows, {len(sequences)} FASTA chains')
    print(f'Created {len(variants)} unique single-AA variants')
    print(f'FASTA: {fasta_path}; metadata: {tsv_path}')
    print(f'Skipped {len(skipped)} input rows/position checks; details: {skips_path}')
    if skipped:
        print('WARNING: Inspect skipped rows. Chain identities must match and map uniquely.')


if __name__ == '__main__':
    main()
