
# Commands

```bash
awk -v OFS="\t" '{ print $6, $7, $8, $4, $1, $5 }' HG00171_cens.bed > HG00171_cens_ordered.bed
```
