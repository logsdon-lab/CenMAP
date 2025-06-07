# Identify centromeric regions
# Extract centromeric regions


include: "utils.smk"
include: "common.smk"
include: "1-concat_asm.smk"
include: "4-align_asm_to_ref.smk"


IDENT_CEN_CTGS_OUTDIR = join(OUTPUT_DIR, "5-ident_cen_ctgs")
IDENT_CEN_CTGS_LOGDIR = join(LOG_DIR, "5-ident_cen_ctgs")
IDENT_CEN_CTGS_BMKDIR = join(BMK_DIR, "5-ident_cen_ctgs")


# Convert rustybam stats bedfile by adjusting start and end positions.
rule format_hor_ref_aln_cen_contigs:
    input:
        aln_bed=expand(
            rules.asm_ref_aln_to_bed.output, ref=f"{REF_NAME}_cens", sm="{sm}"
        ),
    output:
        cen_regions=join(
            IDENT_CEN_CTGS_OUTDIR,
            "bed",
            "interm",
            "{sm}_cens.bed",
        ),
    conda:
        "../envs/tools.yaml"
    log:
        join(IDENT_CEN_CTGS_LOGDIR, "format_hor_ref_aln_cen_contigs_{sm}.log"),
    shell:
        """
        awk -v OFS="\\t" '{{
            if (NR == 1) {{
                print;
                next;
            }} 
            # Find starts/ends in contig name.
            match($1, ":(.+)-", ref_starts);
            # Remove coords from ctg name
            gsub(":.*-.*", "", $1)
            # Print columns.
            $2=$2+ref_starts[1]
            $3=$3+ref_starts[1]
            print
        }}' {input} > {output} 2> {log}
        """


# Map each centromeric contig to a chromosome.
rule map_collapse_cens:
    input:
        script=workflow.source_path("../scripts/map_cens.py"),
        regions=rules.format_hor_ref_aln_cen_contigs.output,
        fai=rules.concat_asm.output.idx,
    output:
        cens_key=join(
            IDENT_CEN_CTGS_OUTDIR,
            "bed",
            "interm",
            "{sm}_mapped_cens.bed",
        ),
        # old_name, new_name, coords, sample, chrom, is_reverse
        renamed_cens_key=join(
            IDENT_CEN_CTGS_OUTDIR,
            "bed",
            "interm",
            "{sm}_renamed_cens.tsv",
        ),
    params:
        # TODO: Should also affect humas-sd monomer lib creation.
        allow_multi_chr_prop=33.0,
    conda:
        "../envs/py.yaml"
    log:
        join(IDENT_CEN_CTGS_LOGDIR, "map_collapse_cens_{sm}.log"),
    shell:
        """
        python {input.script} -i {input.regions} --allow_multi_chr_prop {params.allow_multi_chr_prop} > {output.cens_key} 2> {log}
        awk -v OFS="\\t" '{{
            # 1-based index for seqtk
            st=$2 + 1;
            # Query length. seqtk caps at contig length.
            # Cannot use rb stats as only centromere length.
            end=($3 > $7) ? $7 : $3;
            coords=st"-"end;
            old_name=$1;
            new_name="{wildcards.sm}_"$4"_"$1;
            print old_name,new_name,coords,"{wildcards.sm}",$4,$6
        }}' <(join <(sort -k1,1 {output.cens_key}) <(cut -f1,2 {input.fai} | sort -k1,1)) > {output.renamed_cens_key} 2> {log}
        """


# Extract complete centromeric contigs.
rule extract_cens_regions:
    input:
        bed=rules.map_collapse_cens.output.cens_key,
        fa=rules.concat_asm.output.fa,
    output:
        seq=temp(
            join(
                IDENT_CEN_CTGS_OUTDIR,
                "seq",
                "interm",
                "{sm}_centromeric_regions.fa",
            )
        ),
        idx=temp(
            join(
                IDENT_CEN_CTGS_OUTDIR,
                "seq",
                "interm",
                "{sm}_centromeric_regions.fa.fai",
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
            rules.map_collapse_cens.output,
            sm=SAMPLE_NAMES,
        ),
    default_target: True
