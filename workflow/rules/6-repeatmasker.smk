
include: "common.smk"
include: "utils.smk"
include: "5-ident_cen_ctgs.smk"


RM_OUTDIR = join(OUTPUT_DIR, "6-repeatmasker")
RM_LOGDIR = join(LOG_DIR, "6-repeatmasker")
RM_BMKDIR = join(BMK_DIR, "6-repeatmasker")


rule setup_repeatmasker:
    output:
        chkpt=touch(join(RM_OUTDIR, "rm_setup.done")),
        seq=join(RM_OUTDIR, "rm_setup.fa"),
        rm_dir=directory(join(RM_OUTDIR, "rm_setup")),
    params:
        species=config["repeatmasker"]["species"],
        engine="rmblast",
    threads: 1
    log:
        join(RM_LOGDIR, "setup_repeatmasker.log"),
    conda:
        "../envs/tools.yaml"
    shell:
        """
        echo ">rm_setup" > {output.seq}
        echo "NNNNNNNNNNNNNNNNNNNNN" >> {output.seq}
        RepeatMasker \
            -engine {params.engine} \
            -species {params.species} \
            -dir {output.rm_dir} \
            -pa {threads} \
            {output.seq} &> {log}
        """


checkpoint split_cens_for_rm:
    input:
        fa=rules.extract_cens_regions.output.seq,
    output:
        directory(join(RM_OUTDIR, "seq", "interm", "{sm}")),
    log:
        join(RM_LOGDIR, "split_{sm}_cens_for_rm.log"),
    params:
        split_dir=lambda wc, output: output[0],
    conda:
        "../envs/tools.yaml"
    shell:
        # https://gist.github.com/astatham/621901
        """
        mkdir -p {params.split_dir}
        awk '{{
            if (substr($0, 1, 1)==">") {{
                ctg_name=substr($0,2)
                filename=("{params.split_dir}/" ctg_name ".fa")
            }}
            print $0 > filename
        }}' {input.fa} 2> {log}
        """


# RepeatMasker has a limit of 50 characters for sequence names.
rule rename_for_repeatmasker:
    input:
        fa=join(rules.split_cens_for_rm.output[0], "{fname}.fa"),
    output:
        original_fa_idx=temp(
            join(RM_OUTDIR, "seq", "{sm}_renamed", "{fname}_original.fa.fai"),
        ),
        renamed_fa=temp(
            join(
                RM_OUTDIR,
                "seq",
                "{sm}_renamed",
                "{fname}.fa",
            )
        ),
        renamed_fa_idx=temp(
            join(
                RM_OUTDIR,
                "seq",
                "{sm}_renamed",
                "{fname}.fa.fai",
            )
        ),
    params:
        prefix="seq",
    conda:
        "../envs/tools.yaml"
    log:
        join(RM_LOGDIR, "rename_for_repeatmasker_{sm}_{fname}.log"),
    shell:
        """
        samtools faidx {input.fa} -o {output.original_fa_idx} 2> {log}
        seqtk rename {input.fa} {params.prefix} > {output.renamed_fa} 2>> {log}
        if [ -s {output.renamed_fa} ]; then
            samtools faidx {output.renamed_fa} 2>> {log}
        else
            touch {output.renamed_fa_idx}
        fi
        """


rule run_repeatmasker:
    input:
        setup=rules.setup_repeatmasker.output,
        seq=rules.rename_for_repeatmasker.output.renamed_fa,
    output:
        temp(
            join(
                RM_OUTDIR,
                "repeats",
                "{sm}_renamed",
                "{fname}.fa.out",
            )
        ),
    threads: config["repeatmasker"]["threads"]
    params:
        output_dir=lambda wc, output: os.path.dirname(str(output)),
        species=config["repeatmasker"]["species"],
        engine="rmblast",
    conda:
        "../envs/tools.yaml"
    log:
        join(RM_LOGDIR, "repeatmasker_{sm}_{fname}.log"),
    benchmark:
        join(RM_BMKDIR, "repeatmasker_{sm}_{fname}.tsv")
    shell:
        """
        if [ -s {input.seq} ]; then
            RepeatMasker \
            -engine {params.engine} \
            -species {params.species} \
            -dir {params.output_dir} \
            -qq \
            -pa {threads} \
            {input.seq} &> {log}
        else
            touch {output}
        fi
        """


# Rename repeatmasker output to match the original sequence names.
rule reformat_repeatmasker_output:
    input:
        rm_out=rules.run_repeatmasker.output,
        original_fai=rules.rename_for_repeatmasker.output.original_fa_idx,
        renamed_fai=rules.rename_for_repeatmasker.output.renamed_fa_idx,
    output:
        join(
            RM_OUTDIR,
            "repeats",
            "{sm}",
            "{fname}.fa.out",
        ),
    params:
        script=workflow.source_path("../scripts/reformat_rm.py"),
    log:
        join(RM_LOGDIR, "reformat_repeatmasker_output_{sm}_{fname}.log"),
    conda:
        "../envs/py.yaml"
    shell:
        """
        python {params.script} -i {input.rm_out} -of {input.original_fai} -rf {input.renamed_fai} > {output} 2> {log}
        """


# Gather all RM output
def refmt_rm_output(wc):
    _ = checkpoints.split_cens_for_rm.get(**wc).output
    fa_glob_pattern = join(
        expand(rules.split_cens_for_rm.output, sm=wc.sm)[0], "{fname}.fa"
    )
    wcs = glob_wildcards(fa_glob_pattern)
    fnames = wcs.fname
    files = expand(rules.reformat_repeatmasker_output.output, sm=wc.sm, fname=fnames)
    if not files:
        raise RuntimeError(f"No {wc.sm} centromeres found.")
    return files


rule format_repeatmasker_output:
    input:
        rm_out=refmt_rm_output,
        # Force snakemake to not evaluate chkpt function until all dirs created.
        rm_fa_dirs=expand(rules.split_cens_for_rm.output, sm=SAMPLE_NAMES),
    output:
        join(
            RM_OUTDIR,
            "repeats",
            "all",
            "{sm}_cens.fa.out",
        ),
    conda:
        "../envs/tools.yaml"
    shell:
        """
        {{ awk -v OFS="\\t" '{{$1=$1; print}}' {input.rm_out} | cut -f 1-15 ;}} > {output}
        """


include: "6.1-plot_repeatmasker_chr.smk"
include: "6.1-plot_repeatmasker_sm.smk"


RM_OUTPUTS = []
if "chromosome" in config["repeatmasker"].get("partition_by", []):
    RM_OUTPUTS.extend(
        expand(
            rules.plot_og_rm_bed_by_chr.output,
            chr=CHROMOSOMES if CHROMOSOMES else ["all"],
        )
    )
if "sample" in config["repeatmasker"].get("partition_by", []):
    RM_OUTPUTS.extend(
        expand(
            rules.plot_og_rm_bed_by_sm.output,
            sm=SAMPLE_NAMES,
        )
    )


rule repeatmasker_all:
    input:
        RM_OUTPUTS,
    default_target: True
