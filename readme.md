# Protein Drug: infliximab
 
**Research Title:** Proposed De-Immunisation of Therapeutic Monoclonal Antibody Infliximab
 
## pipeline:
 
1. Retrieve the antibody sequences and structure
Thera-SAbDab + PDB 4G3Y (may be a patterned seq so we need the hole Ab)
 
2. Identify antibody regions that should be protected
PyMOL: CDRs and TNF-α-contacting residues
 
3. Predict baseline MHC-II-binding peptides
NetMHCIIpan-4.3 across a selected HLA-II panel
 
4. Narrow down HLA set to most prominent in population
PopCover on predicted peptides from step 3
 
5. Propose targeted amino acid substitutions
Mutate anchor points in each predicted peptide
 
6. Recalculate immunogenicity for mutant sequences
Compare wild type against mutants using the same HLA panel
 
7. Assess structural compatibility
PyMOL + FoldX or Rosetta, with a small set of candidates
 
8. Rank de-immunisation candidates
Balance predicted epitope reduction against predicted structural disruption (Use FoldX for ranking)