# HGSVC3
WIP

![](docs/rulegraph.svg)

### Usage
Local
```bash
# NOTE: dna-brnn must be installed locally
snakemake --use-conda -np --configfile config/config.yaml
```

Cluster via `singularity`. In project dir.
```bash
singularity pull docker://koisland/hgsvc3:latest
singlarity run --use-conda -np --configfile config/config.yaml
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
* Remove old bedminmax scripts. Only leaving in until output discrepancies finished.
