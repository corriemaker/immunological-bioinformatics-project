 
```sh 
 awk '$1 ~ /^[0-9]+$/ && $9 <= 5 {print $2}' \
    669056_NetMHCIIpan.txt | sort | uniq -c
  41 DRB1_0701
  28 DRB1_1101
  20 DRB1_1501
  10 HLA-DQA10101-DQB10501
  35 HLA-DQA10201-DQB10202
  32 HLA-DQA10501-DQB10201
```