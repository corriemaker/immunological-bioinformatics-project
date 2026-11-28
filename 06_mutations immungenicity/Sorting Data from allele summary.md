Tested against your actual `step6_summary.tsv`. Three levels, depending on how strict you want to be.

**1. Clean wins — every target allele cleared, nothing new created:**

```bash
awk -F'\t' 'NR>1 && $8=="CLEARED" && $12==0' step6_summary.tsv | sort -t$'\t' -k9,9nr -k14,14nr
```

Sorts by `n_target_alleles` (how much the mutation bought you), then by the new rank. Gives 15 rows — top two are `Y192D` (5 alleles cleared, 0.08 → 2.77) and `W148D` (4 alleles, 0.02 → 3.34).

**2. Same, restricted to the variable domains** — since constant-domain hits are human self-sequence and won't make the final list:

```bash
awk -F'\t' 'NR>1 && $8=="CLEARED" && $12==0 {
  match($4,/[0-9]+/); p=substr($4,RSTART,RLENGTH)+0
  lim = ($3 ~ /1_Chain/) ? 107 : 120
  if (p<=lim) print
}' step6_summary.tsv | sort -t$'\t' -k9,9nr -k14,14nr
```

This is the one that matters — it cuts 15 down to **6**.

**3. With CDR flagging**, since two of those six sit on the paratope:

```bash
awk -F'\t' '
NR==1 { next }
$8=="CLEARED" && $12==0 {
  match($4,/[0-9]+/); p = substr($4,RSTART,RLENGTH)+0
  if ($3 ~ /1_Chain/) { lim=107; cdr=(p>=27&&p<=38)||(p>=50&&p<=52)||(p>=89&&p<=97) }
  else                { lim=120; cdr=(p>=27&&p<=33)||(p>=51&&p<=59)||(p>=99&&p<=109) }
  if (p<=lim) printf "%-38s %-12s %4d %-10s %4s %8s %8s\n", $2, $3, p, (cdr?"CDR":"framework"), $9, $13, $14
}' step6_pssm_summary.tsv | sort -k5,5nr -k7,7nr
```

```
label                                  chain         pos region      tgt   WTrank  MUTrank
4G3Y_2_Chain__V12D                     4G3Y_2_Chain   12 framework     2     0.06     7.67
4G3Y_2_Chain__D87M                     4G3Y_2_Chain   87 framework     2     0.26     6.80
4G3Y_2_Chain__T90W                     4G3Y_2_Chain   90 framework     2     0.26     5.68
4G3Y_1_Chain__Y50Q                     4G3Y_1_Chain   50 CDR           1     0.15     5.82
4G3Y_1_Chain__L47R                     4G3Y_1_Chain   47 framework     1     0.15     2.35
```

**Four framework candidates to carry into step 7** — that's your actual result. `V12D` is the standout: two alleles cleared, rank 0.12 → 11.87, VH framework 1, nowhere near the paratope.

From the by-allele file, if you want to see which restrictions each one removed:

```bash
awk -F'\t' 'NR==1 || $4=="cleared"' step6_by_allele.tsv | grep -E 'V12D|D87L|T90F|S7W'
```

**One caveat on the CDR ranges** — those are Kabat numbering hardcoded from the sequence, not from your structure. The authoritative list is the paratope in `02_`, which is still only a `.pse`. Once you export that to a residue list you can swap the hardcoded ranges for a join against the real file, which will also catch TNF-α contacts in the framework that Kabat ranges miss.
