# Docs

### Make `rulegraph.svg`
```bash
snakemake -np \
--configfile test/config/config.yaml \
--rulegraph | \
sed 's/snakemake_dag {/snakemake_dag { ranksep=0.3;/g' | \
dot -Tsvg > docs/rulegraph.svg
```
