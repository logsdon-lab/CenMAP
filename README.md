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

#### Download raw data with `globus-cli`
Login and check your globus home endpt direcory. Directories specified below are local to this.
```bash
globus login
globus endpoint local-id
# C:\Users\koshima
```

To download files.
```bash
./workflow/scripts/globus_download_data.sh projects/hgsvc3/data/asm projects/hgsvc3/data/raw_data
```

### TODO
* TODO: Move renaming to start of analysis. Rename `sampleid_hapinfo_ctg-#`
