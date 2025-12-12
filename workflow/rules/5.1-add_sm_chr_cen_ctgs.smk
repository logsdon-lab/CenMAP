# Identify centromeric regions
# Extract centromeric regions


include: "utils.smk"
include: "common.smk"
include: "2-concat_asm.smk"
include: "3-srf.smk"
include: "4-align_asm_to_ref.smk"


# Map each contig to a chromosome.
rule create_rename_key:
    input:
        alns=expand(rules.asm_ref_aln_to_bed.output, ref=REF_NAME, sm="{sm}"),
    output:
        # (old_name, new_name, length)
        rename_key=join(
            IDENT_CEN_CTGS_OUTDIR,
            "bed",
            "{sm}_rename_key.tsv",
        ),
    params:
        script=workflow.source_path("../scripts/map_chroms.py"),
        allow_multi_chr_prop=config["ident_cen_ctgs"]["perc_multi_chrom"],
    conda:
        "../envs/py.yaml"
    log:
        join(IDENT_CEN_CTGS_LOGDIR, "create_rename_key_{sm}.log"),
    shell:
        """
        python {params.script} -i {input.alns} --allow_multi_chr_prop {params.allow_multi_chr_prop} | \
        awk -v OFS="\\t" '{{
            old_name=$1;
            chrom=$2;
            rc=($3== "true") ? "rc-" : "";
            ctg_len=$4;
            new_name="{wildcards.sm}_"rc""chrom"_"old_name
            print old_name, new_name, ctg_len
        }}' > {output.rename_key} 2> {log}
        """


rule create_final_asm:
    input:
        fa=rules.concat_asm.output.fa,
        idx=rules.concat_asm.output.idx,
        rename_key=rules.create_rename_key.output,
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
        join(IDENT_CEN_CTGS_LOGDIR, "create_{sm}_final_asm.log"),
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
rule create_final_satellite_bed:
    input:
        bed=rules.merge_slop_region_bed.output,
        rename_key=rules.create_rename_key.output,
    output:
        bed=join(
            IDENT_CEN_CTGS_OUTDIR,
            "bed",
            "{sm}.bed",
        ),
    conda:
        "../envs/tools.yaml"
    log:
        join(IDENT_CEN_CTGS_LOGDIR, "create_final_satellite_bed_{sm}.log"),
    shell:
        """
        {{ join <(sort -k1,1 {input.bed} | cut -f 1-4) <(sort -k1,1 {input.rename_key}) | \
        awk -v OFS="\\t" '{{
            ctg_len=$6;
            if ($5 ~ "rc-chr") {{
                st=ctg_len-$3;
                end=ctg_len-$2;
            }} else {{
                st=$2;
                end=$3;
            }}
            print $5, st, end
        }}' | \
        sort -k1,1 -k2,2n | \
        bedtools merge -i - ;}} > {output.bed} 2> {log}
        """


# Extract complete centromeric contigs.
rule extract_cens_regions:
    input:
        bed=rules.create_final_satellite_bed.output,
        fa=rules.create_final_asm.output.fa,
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
        bed=lambda wc, input: input.bed,
        # Filter by chromosome and covert to upper-case
        added_cmds=(
            f"| {cmd_filter_fa_chrom("seqkit seq --upper-case")}"
            if CHROMOSOMES
            else "| seqkit seq --upper-case"
        ),
    log:
        join(IDENT_CEN_CTGS_LOGDIR, "extract_regions_{sm}.log"),
    conda:
        "../envs/tools.yaml"
    shell:
        shell_extract_and_index_fa


rule ident_cen_ctgs_all:
    input:
        expand(
            rules.create_rename_key.output,
            sm=SAMPLE_NAMES,
        ),
        expand(
            rules.create_final_asm.output,
            sm=SAMPLE_NAMES,
        ),
    default_target: True
