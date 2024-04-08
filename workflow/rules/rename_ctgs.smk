# Module to rename oriented regions in an assembly.
#
# ### Params:
#   * bed_input_regions
#       * Bed file with regions. Filename is split by '\t' and ':'.
#       * Requires sm and ort wilcards
#   * fa_assembly
#       * Assembly with contig names to replace.
#   * output_dir
#   * log_dir
#   * samples
#   * orientation
#       * Default: ("fwd", "rev")
#   * bed_find_col
#       * Column of values to replace in assembly contig names.
#       * Default: 1
#   * bed_replace_w_joined_cols
#       * Columns of values to join into string to replace bed_find_col.
#       * Default: (7,4,1))
#       * Delimited by '_'

import os


SAMPLES = config["samples"]
OUTPUT_DIR = config.get("output_dir", "output")
LOG_DIR = config.get("logs_dir", "logs")
BED_FIND_COL = config.get("bed_find_col", 1)
BED_REPLACE_W_JOINED_COLS = config.get("bed_replace_w_joined_cols", (7, 4, 1))
# Put renamed bed in same dir.
BED_INPUT_DIR = os.path.dirname(str(config["bed_input_regions"]))


# 1. name: haplotype1-0000027
# 2. start: 96023560
# 3. end: 101450776
# 4. chr: chr7
# 5. orientation: +
# 6. start_end_diff: 5427216
# 7. sample: HG00171
# haplotype1-0000027    HG00171_chr7_haplotype1-0000027
rule create_renamed_bed_n_legend:
    input:
        regions=config["bed_input_regions"],
    output:
        legend=os.path.join(OUTPUT_DIR, "{sm}_legend.txt"),
        regions_renamed=os.path.join(BED_INPUT_DIR, "{sm}_renamed.bed"),
    log:
        os.path.join(LOG_DIR, "create_legend_{sm}.log"),
    params:
        legend_key=f"${BED_FIND_COL}",
        legend_val='"_"'.join(f"${col}" for col in BED_REPLACE_W_JOINED_COLS),
    conda:
        "../env/tools.yaml"
    shell:
        """
        {{ awk -v OFS="\\t" '{{print $0, "{wildcards.sm}"}}' {input.regions} | \
        awk -v OFS="\\t" '{{
            print {params.legend_key},{params.legend_val};
            print {params.legend_val}, $2, $3, $4, $5, $6, $7 >> "{output.regions_renamed}"
        }}' | \
        sort -k2,2 | uniq ;}} > {output.legend} 2> {log}
        """


# Before:
# >haplotype1-0000003:4-8430174
# After:
# >HG00171_chr16_haplotype1-0000003:4-8430174
rule rename_ctgs:
    input:
        legend=rules.create_renamed_bed_n_legend.output.legend,
        seq=config["fa_assembly"],
    output:
        os.path.join(
            OUTPUT_DIR,
            "{sm}_regions.renamed.fa",
        ),
    run:
        with (
            open(str(input.legend)) as legend_fh,
            open(str(input.seq)) as seq_fh,
            open(str(output), "wt") as out_seq_fh,
        ):
            legend = dict(line.strip().split() for line in legend_fh.readlines())

            for line in seq_fh.readlines():
                # >h1tg000001l#1-110442987
                # >haplotype1-0000001:56723941-58474753
                if line.startswith(">"):
                    seq_id, *coords = line.strip().strip(">").partition(":")
                    new_seq_id = legend.get(seq_id, seq_id)
                    out_seq_fh.write(f">{new_seq_id}{''.join(coords)}\n")
                    # The seq.
                else:
                    out_seq_fh.write(line)


rule index_renamed_ctgs:
    input:
        rules.rename_ctgs.output,
    output:
        os.path.join(
            OUTPUT_DIR,
            "{sm}_regions.renamed.fa.fai",
        ),
    conda:
        "../env/tools.yaml"
    log:
        os.path.join(LOG_DIR, "index_renamed_ctgs_{sm}.log"),
    shell:
        """
        samtools faidx {input} &> {log}
        """


rule rename_ctg_all:
    input:
        expand(
            rules.create_renamed_bed_n_legend.output,
            sm=SAMPLES,
        ),
        expand(
            rules.rename_ctgs.output,
            sm=SAMPLES,
        ),
        expand(rules.index_renamed_ctgs.output, sm=SAMPLES),
