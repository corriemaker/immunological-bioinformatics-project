python3 compare_mutant_predictions.py \
    --predictions netMHCIIpan_mutants_out.txt \
    --map mutants.tsv \
    --threshold 2.0 \
    --out step6_summary.tsv --details step6_by_allele.tsv

python3 compare_mutant_predictions.py --predictions netMHCIIpan_mutants_out.txt --map ../05_mutations\ proposed/mutants.tsv --threshold 2.0 --out step6_summary.tsv --details step6_by_allele.tsv