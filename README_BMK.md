
Run NucFlag benchmarks.
```bash
snakemake -j 200 --configfile config/config_hprc.yaml \
-p nucflag_bmk_only
```

To prevent downloading data if already present.
```bash
touch download.done
```
