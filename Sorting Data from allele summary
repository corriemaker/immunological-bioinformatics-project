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
  if ($3 ~ /1_Chain/) { lim=107; cdr=(p>=24&&p<=34)||(p>=50&&p<=56)||(p>=89&&p<=97) }
  else                { lim=120; cdr=(p>=31&&p<=35)||(p>=50&&p<=65)||(p>=95&&p<=109) }
  if (p<=lim) printf "%-38s %-12s %4d %-10s %4s %8s %8s\n", $2, $3, p, (cdr?"CDR":"framework"), $9, $13, $14
}' step6_summary.tsv | sort -k5,5nr -k7,7nr
```

```
label                                  chain         pos region      tgt   WTrank  MUTrank
4G3Y_2_Chain__V12D                     4G3Y_2_Chain   12 framework     2     0.12    11.87
4G3Y_2_Chain__S7W_G10I_V12D_G15I       4G3Y_2_Chain    7 framework     2     0.12     6.47
4G3Y_2_Chain__D87L                     4G3Y_2_Chain   87 framework     1     0.29     5.79
4G3Y_2_Chain__R52F                     4G3Y_2_Chain   52 CDR           1     0.33     4.63
4G3Y_2_Chain__T90F                     4G3Y_2_Chain   90 framework     1     0.29     3.64
4G3Y_1_Chain__S52W                     4G3Y_1_Chain   52 CDR           1     0.28     3.40
```

**Four framework candidates to carry into step 7** — that's your actual result. `V12D` is the standout: two alleles cleared, rank 0.12 → 11.87, VH framework 1, nowhere near the paratope.

From the by-allele file, if you want to see which restrictions each one removed:

```bash
awk -F'\t' 'NR==1 || $4=="cleared"' step6_by_allele.tsv | grep -E 'V12D|D87L|T90F|S7W'
```

**One caveat on the CDR ranges** — those are Kabat numbering hardcoded from the sequence, not from your structure. The authoritative list is the paratope in `02_`, which is still only a `.pse`. Once you export that to a residue list you can swap the hardcoded ranges for a join against the real file, which will also catch TNF-α contacts in the framework that Kabat ranges miss.
