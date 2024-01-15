# HGSVC3
WIP

### Usage
Local
```bash
# NOTE: dna-brnn must be installed locally
snakemake --use-conda -np --configfile config/config.yaml
```

Cluster via `singularity`
```bash
singularity pull docker://koisland/hgsvc3:latest
singlarity run --use-conda -np --configfile config/config.yaml
```

Make test files for workflow.
```bash
./test/make_test_files.sh
```

### TODO
* TODO: Move renaming to start of analysis. Rename `sampleid_hapinfo_ctg-#`
