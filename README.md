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
# Make fake ref
ref=data/chm13_t2t/nobackups/assemblies/chm13_v1.1_plus38Y_par_masked.fasta
mkdir -p $(dirname $ref) && touch $ref

# Make fake assemblies
awk 'NR > 1 { print $2 }' config/table.asm.tbl | \
    xargs -n1 -I [] bash -c 'mkdir -p $(dirname []) && touch []'
```

### TODO
* TODO: Move renaming to start of analysis. Rename `sampleid_hapinfo_ctg-#`
