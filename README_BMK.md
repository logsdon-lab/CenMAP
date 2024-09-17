
Run NucFlag benchmarks.
```bash
snakemake -j 200 --configfile config/config_hg002.yaml \
-p nucflag_bmk_only \
--workflow-profile workflow/profiles/lpc
```

If data is already available, to prevent redownloading, create a `download_data.done` in the main dir.
```bash
touch download_data.done
```
