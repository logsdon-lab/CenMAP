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
#       * Column of values to replace in assembly contig names. After filename split.
#       * Default: 1
#   * bed_replace_w_joined_cols
#       * Columns of values to join into string to replace bed_find_col. After filename split
#       * Default: (9,5,1))
#       * Delimited by '_'

import os


SAMPLES = config["samples"]
OUTPUT_DIR = config.get("output_dir", "output")
LOG_DIR = config.get("logs_dir", "logs")
BED_FIND_COL = config.get("bed_find_col", 1)
BED_REPLACE_W_JOINED_COLS = config.get("bed_replace_w_joined_cols", (9, 5, 1))
ORIENTATION = config.get("orientation", ("fwd", "rev"))
SED_CMD = config.get(
    "sed_cmd", "sed -e 's/> />/g' -e 's/\([0-9]\) \([0-9]\)/\\1:\\2/g'"
)


# 1. name: haplotype1-0000027
# 2. start: 96023560
# 3. end: 101450776
# 4. length: 158780978
# 5. chr: chr7
# 6. chr_coords: 58924390-64604808
# 7. orientation: +
# 8. start_end_diff: 5427216
# 9-12. results/cens/HG00171    1       centromeric     regions.fwd.bed
# haplotype1-0000027    HG00171_chr7_haplotype1-0000027
rule create_oriented_ctg_name_legend:
    input:
        regions=config["bed_input_regions"],
    output:
        os.path.join(OUTPUT_DIR, "{sm}.legend.{ort}.txt"),
    log:
        os.path.join(LOG_DIR, "create_{ort}_ctg_name_legend_{sm}.log"),
    params:
        # Replaced awk FILENAME with (vvv) because run from multiple dirs above.
        # Would include subdirs otherwise.
        # test/HG00171_1_centromeric_regions.rev.bed -> HG00171_1_centromeric_regions
        file_bname=lambda wc, input: os.path.basename(str(input)).split(".")[0],
        legend_key=f"${BED_FIND_COL}",
        legend_val='"_"'.join(f"${col}" for col in BED_REPLACE_W_JOINED_COLS),
    conda:
        "../env/tools.yaml"
    shell:
        """
        {{ awk -v OFS="\\t" '{{print $0, "{params.file_bname}"}}' {input.regions} | \
        sed -e 's/_/\\t/g' -e 's/:/\\t/g' | \
        awk -v OFS="\\t" '{{print {params.legend_key},{params.legend_val}}}' | \
        sort -k2,2 | uniq ;}} > {output} 2> {log}
        """


# :Before:
# >haplotype1-0000001:56723941-58474753
# :After:
# >
# haplotype1-0000001
# 56723941-58474753
rule split_oriented_assembly_fasta:
    input:
        config["fa_assembly"],
    output:
        os.path.join(OUTPUT_DIR, "{sm}.{ort}.txt"),
    conda:
        "../env/tools.yaml"
    log:
        os.path.join(LOG_DIR, "split_{ort}_assembly_fasta_{sm}.log"),
    shell:
        """
        sed -e 's/>/>\\n/g' -e 's/:/\\n/g' {input} > {output} 2> {log}
        """


# Before:
# >
# haplotype1-0000003
# 4-8430174
# After:
# >HG00171_chr16_haplotype1-0000003:4-8430174
rule rename_oriented_ctgs:
    input:
        legend=rules.create_oriented_ctg_name_legend.output,
        split_seq=rules.split_oriented_assembly_fasta.output,
    output:
        os.path.join(
            OUTPUT_DIR,
            "{sm}_regions.renamed.{ort}.fa",
        ),
    conda:
        "../env/tools.yaml"
    params:
        sed_cmd=SED_CMD,
    log:
        os.path.join(LOG_DIR, "rename_{ort}_ctgs_{sm}.log"),
    shell:
        # Construct associative array from input legend. (line 2)
        """
        {{ awk 'BEGIN{{FS=OFS="\\t"}} \
        NR==FNR {{legend[$1]=$2; next}} \
        {{print ($1 in legend ? legend[$1] : $1)}}' {input.legend} {input.split_seq} | \
        awk '{{printf "%s%s", (/>/ ? ors : OFS), $0; ors=ORS}} END{{print ":"}}' | \
        {params.sed_cmd};}} > {output} 2> {log}
        """


rule index_renamed_ctgs:
    input:
        rules.rename_oriented_ctgs.output,
    output:
        os.path.join(
            OUTPUT_DIR,
            "{sm}_regions.renamed.{ort}.fa.fai",
        ),
    conda:
        "../env/tools.yaml"
    log:
        os.path.join(LOG_DIR, "index_renamed_{ort}_ctgs_{sm}.log"),
    shell:
        """
        samtools faidx {input} &> {log}
        """


rule rename_ctg_all:
    input:
        expand(
            rules.create_oriented_ctg_name_legend.output,
            sm=SAMPLES,
            ort=ORIENTATION,
        ),
        expand(
            rules.split_oriented_assembly_fasta.output,
            sm=SAMPLES,
            ort=ORIENTATION,
        ),
        expand(
            rules.rename_oriented_ctgs.output,
            sm=SAMPLES,
            ort=ORIENTATION,
        ),
        expand(rules.index_renamed_ctgs.output, sm=SAMPLES, ort=ORIENTATION),
