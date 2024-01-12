# HGSVC3
WIP

### Usage
```bash
snakemake --use-conda -np --configfile config/config.yaml
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
