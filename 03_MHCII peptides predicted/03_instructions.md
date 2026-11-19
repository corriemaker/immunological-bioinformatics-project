## Step 1) Code to create popcover MHC-II input file
awk '$1 ~ /^[0-9]+$/ && $9 <= 1 {print $3, $2, substr($7,1,12), $5}' \
    full_allele_list_netMHCIIpan_out.txt > popcover_inp_classII

### Code to get HLA alleles and predicted number of binders
awk '$1 ~ /^[0-9]+$/ && $9 <= 1 {print $2}' \
    full_allele_list_netMHCIIpan_out.txt | sort | uniq -c

### Output:
      3 DRB1_1101
      4 DRB1_1501
      2 HLA-DQA10103-DQB10301
     11 HLA-DQA10201-DQB10202
     16 HLA-DQA10401-DQB10301
      9 HLA-DQA10501-DQB10201

## Step 2) run ./afnd_popcover_exact_match.py after changing the .py file to take the above output

### Code to fix the output format: <-- Run after ./afnd_popcover_exact_match.py
awk '{$3 = 2; print}' popcover_allele_frequencies.txt > popcover_allele_frequencies_fixed.txt

## Step 3) Go to PopCover and supply the following files:
- popcover_allele_frequencies_fixed.txt
- popcover_inp_classII
- 4G3Y_1_Chain and 4G3Y_2_Chain (in the "Include protein sequence information (optional)" menu) --- 4G3Y_3_Chain is not the protein but the target so we need to NOT include it 
Under Additional customization select "Use phenotypic frequencies for calculation"