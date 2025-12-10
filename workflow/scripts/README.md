# Scripts

## `cenplot*.yaml`
Various configuration for `cenplot` scripts. Can be modified to change layout of scripts.

## `count_complete_cens.py`
Count complete centromeres. Assumes diploid # of chromosomes by default but can be altered if needed.

## `create_rm_bed.py`
Takes headerless `RepeatMasker` output and JSON file of regex patterns and colors (See `config/repeatmasker_colors.json`) and outputs a simplified `RepeatMasker` bedfile. If desired, also adjusts coordinates to absolute if name column follows format `chrom:st-end`.

## `filter_entropy_bed.py`
Determine if centromeres are valid based on shannon index and ALR/Alpha overlap. A complete centromere should have dip(s) in entropy overlapping ALR/Alpha that account for the majority (or more) of the length of ALR/Alpha annotated array. Also performs trimming if needed.

## `format_cenplot_yaml.py`
Take input files as a JSON string and a `cenplot` track config and generate a new configfile.

## `format_rm_sat_annot.py`
Take headerless `RepeatMasker` and JSON file of groups with patterns/colors (See `config/repeatmasker_sat_annot_colors.json`) and generates BED9 file.

## `make_env_all.sh`
Generate conda environment using conda-merge. Will need to be manually fixed to remove duplicate dependencies.

## `map_chroms.py`
Maps chromosomes from an assembly to reference alignment bed file produced by `rustybam stats`.

## `plot_complete_cen_counts.py`
Plots complete centromere counts. Takes input from [`count_complete_cens.py`](#count_complete_censpy).

## `plot_hor_length.py`
Plots HOR array length. Takes `censtats length` output.

## `plot_ideogram.py`
Creates ideogram from faidx lengths and centromere positions.

## `plot_multiple_cen.py`
Generates plot with all JSON string of input files. Each key-value pair should correspond to a path in the `cenplot` track config.

```json
{"cdr": ["results/8-cdr_finder/bed/sm_cdrs.bed"]}
```
```yaml
tracks:
  - position: relative
    type: label
    proportion: 0.005
    path: "cdr" # <- HERE
    options:
      color: black
      legend: false
      hide_x: true
```

This script produces a number of temporary files and a config as input data is partitioned by name. To keep these use `--keep-tempfiles`.

To sort output by chromosome, provide the `--sort_order` flag which takes a single file of chromosome names.

## `reformat_rm.py`
`RepeatMasker` produces awful, non-machine-readable output so we have to clean it up.
* Renaming sequences to their original name.
    * It has an [arbitrary limit of 50 characters for sequence names](https://github.com/Dfam-consortium/RepeatMasker/issues/12) so to ensure no input causes issues, we rename all sequences before running.
* Removing multi-row header.
* Converting delimiters from arbitrary # of spaces to tabs.

## `create_rm_overlay_bed.awk`
Creates a BED file of satellite regions for `NucFlag` to plot.

## `create_rm_nucflag_ignore_bed.py`
Creates a list of non-ALR/Alpha regions for `NucFlag` to ignore. If there are small repeats less that 10 kbp between ALR/Alpha (TEs), these are allowed to be called.
