python worst_anchor_substitution.py \
    --cores-file "03_MHCII peptides predicted/selected_peptide_cores.txt" \
    --core-column 3 --chain-column 7 \
    --fasta "01_raw data/rcsb_pdb_4G3Y.fasta" \
    --out-fasta mutants.fasta --map mutants.tsv --tsv anchors.tsv
