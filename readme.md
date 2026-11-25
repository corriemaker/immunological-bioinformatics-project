# Protein Drug: infliximab
 
**Research Title:** Proposed De-Immunisation of Therapeutic Monoclonal Antibody Infliximab
 
## pipeline:
 
1. Retrieve the antibody sequences and structure
- Thera-SAbDab + PDB 4G3Y (may be a patterned seq so we need the hole Ab)
 
2. Identify antibody regions that should be protected
- PyMOL: CDRs and TNF-α-contacting residues
 
3. Predict baseline MHC-II-binding peptides
- HLA-II set selected based on both research articles and world coverage
- NetMHCIIpan-4.3 across selected HLA-II panel

4. Find optimal subset of HLA+peptide set 
- PopCover on predicted HLA+peptides from step 3
 
5. Propose targeted amino acid substitutions
- Mutate anchor points in each predicted peptide
 
6. Recalculate immunogenicity for mutant sequences
- Recalculate MHC-II-binding after mutation
- Compare wild type against mutants using the same HLA panel
 
7. Assess structural compatibility
- PyMOL + FoldX or Rosetta, with a small set of candidates
 
8. Rank de-immunisation candidates
- Balance predicted epitope reduction against predicted structural disruption (Use FoldX for ranking)
