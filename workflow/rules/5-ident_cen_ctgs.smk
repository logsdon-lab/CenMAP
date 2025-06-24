# Identify centromeric regions
# Extract centromeric regions


include: "utils.smk"
include: "common.smk"
include: "2-concat_asm.smk"
include: "4-align_asm_to_ref.smk"


IDENT_CEN_CTGS_OUTDIR = join(OUTPUT_DIR, "5-ident_cen_ctgs")
IDENT_CEN_CTGS_LOGDIR = join(LOG_DIR, "5-ident_cen_ctgs")
IDENT_CEN_CTGS_BMKDIR = join(BMK_DIR, "5-ident_cen_ctgs")


# Map each contig to a chromosome.
rule map_chroms:
    input:
        script=workflow.source_path("../scripts/map_chroms.py"),
        alns=expand(rules.asm_ref_aln_to_bed.output, ref=REF_NAME, sm="{sm}"),
    output:
        # (old_name, new_name)
        rename_key=join(
            IDENT_CEN_CTGS_OUTDIR,
            "bed",
            "{sm}_rename_key.tsv",
        ),
    params:
        allow_multi_chr_prop=config["ident_cen_ctgs"]["perc_multi_chrom"],
    conda:
        "../envs/py.yaml"
    log:
        join(IDENT_CEN_CTGS_LOGDIR, "map_chroms_{sm}.log"),
    shell:
        """
        python {input.script} -i {input.alns} --allow_multi_chr_prop {params.allow_multi_chr_prop} | \
        awk -v OFS="\\t" '{{
            old_name=$1;
            chrom=$2;
            rc=($3== "true") ? "rc-" : "";
            ctg_len=$4;
            new_name="{wildcards.sm}_"rc""chrom"_"old_name
            print old_name, new_name, ctg_len
        }}' > {output.rename_key} 2> {log}
        """


rule rename_reort_asm:
    input:
        fa=rules.concat_asm.output.fa,
        idx=rules.concat_asm.output.idx,
        rename_key=rules.map_chroms.output,
    output:
        fa=join(
            CONCAT_ASM_OUTDIR,
            "{sm}-asm-renamed-reort.fa",
        ),
        idx=join(
            CONCAT_ASM_OUTDIR,
            "{sm}-asm-renamed-reort.fa.fai",
        ),
    params:
        pattern=r"'^(\S+)\s*'",
        replacement=lambda wc: "'{kv}'",
    conda:
        "../envs/tools.yaml"
    log:
        join(IDENT_CEN_CTGS_LOGDIR, "fix_{sm}_asm_orientation.log"),
    shell:
        """
        # Get the reverse cens and reverse them.
        # Get all the non-reversed contigs.
        # Then replace the names.
        seqkit replace -p {params.pattern} -r {params.replacement} \
        -k {input.rename_key} \
        <(cat \
            <(seqtk subseq {input.fa} \
                <(awk '$2 ~ "rc-chr"' {input.rename_key} | cut -f1) | \
                seqtk seq -r) \
            <(seqtk subseq {input.fa} \
                <(grep -v -f <(awk '$2 ~ "rc-chr"' {input.rename_key} | cut -f1) {input.idx} | cut -f 1)) \
        ) \
        --keep-key > {output.fa} 2> {log}
        samtools faidx {output.fa} 2>> {log}
        """


# Reorient satellite bed.
# Don't wait for rename assembly as can do in parallel.
rule reorient_satellite_bed:
    input:
        bed=rules.slop_region_bed.output,
        rename_key=rules.map_chroms.output,
        idx=rules.concat_asm.output.idx,
    output:
        bed=join(
            IDENT_CEN_CTGS_OUTDIR,
            "bed",
            "{sm}.bed",
        ),
    conda:
        "../envs/tools.yaml"
    log:
        join(IDENT_CEN_CTGS_LOGDIR, "reorient_satellite_bed_{sm}.log"),
    shell:
        """
        {{ join <(sort -k1,1 {input.bed}) <(sort -k1,1 {input.rename_key}) | \
        join - <(cut -f 1,2 {input.idx}| sort -k1,1) | \
        awk -v OFS="\\t" '{{
            ctg_len=$6;
            if ($5 ~ "rc-chr") {{
                st=ctg_len-$3;
                end=ctg_len-$2;
            }} else {{
                st=$2;
                end=$3;
            }}
            print $5, st, end, end-st
        }}' ;}} > {output.bed} 2> {log}
        """


# Extract complete centromeric contigs.
rule extract_cens_regions:
    input:
        bed=rules.reorient_satellite_bed.output,
        fa=rules.rename_reort_asm.output.fa,
    output:
        seq=temp(
            join(
                IDENT_CEN_CTGS_OUTDIR,
                "seq",
                "{sm}_satellite_regions.fa",
            )
        ),
        idx=temp(
            join(
                IDENT_CEN_CTGS_OUTDIR,
                "seq",
                "{sm}_satellite_regions.fa.fai",
            )
        ),
    params:
        **params_shell_extract_and_index_fa,
    log:
        join(IDENT_CEN_CTGS_LOGDIR, "extract_regions_{sm}.log"),
    conda:
        "../envs/tools.yaml"
    shell:
        shell_extract_and_index_fa


rule ident_cen_ctgs_all:
    input:
        expand(
            rules.map_chroms.output,
            sm=SAMPLE_NAMES,
        ),
        expand(
            rules.rename_reort_asm.output,
            sm=SAMPLE_NAMES,
        ),
        expand(
            rules.extract_cens_regions.output,
            sm=SAMPLE_NAMES,
        ),
    default_target: True
