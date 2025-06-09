
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
        species="human",
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
        # [old_name, new_name, coords, sm, chr, is_reversed]
        rename_key=rules.map_collapse_cens.output.renamed_cens_key,
    output:
        directory(join(RM_OUTDIR, "seq", "interm", "{sm}")),
    log:
        join(RM_LOGDIR, "split_{sm}_cens_for_humas_sd.log"),
    params:
        split_dir=lambda wc, output: output[0],
    conda:
        "../envs/tools.yaml"
    shell:
        # https://gist.github.com/astatham/621901
        """
        mkdir -p {params.split_dir}
        awk '{{
            # Read key values in first file.
            if (FNR == NR) {{
                # Add coords to name.
                kv[$1":"$3]=$2":"$3
                next;
            }}
            if (substr($0, 1, 1)==">") {{
                ctg_name=substr($0,2)
                new_ctg_name=kv[ctg_name]
                filename=("{params.split_dir}/" new_ctg_name ".fa")
                print ">"new_ctg_name > filename
            }} else {{
                print $0 > filename
            }}
        }}' {input.rename_key} {input.fa} 2> {log}
        """


# RepeatMasker has a limit of 50 characters for sequence names.
rule rename_for_repeatmasker:
    input:
        fa=join(rules.split_cens_for_rm.output[0], "{fname}.fa"),
    output:
        original_fa_idx=temp(
            join(RM_OUTDIR, "seq", "{sm}_renamed", "{fname}_original.fa.fai"),
        ),
        # Use different dir to avoid greedy wildcard recursively running rule.
        # TODO: rules.split_cens_for_rm.output[0]
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
        species="human",
        engine="rmblast",
    conda:
        "../envs/tools.yaml"
    log:
        join(RM_LOGDIR, "repeatmasker_{sm}_{fname}.log"),
    # Retry in case of .RepeatMaskerCache failure.
    retries: 2
    benchmark:
        join(RM_BMKDIR, "repeatmasker_{sm}_{fname}.tsv")
    shell:
        """
        if [ -s {input.seq} ]; then
            RepeatMasker \
            -engine {params.engine} \
            -species {params.species} \
            -dir {params.output_dir} \
            -pa {threads} \
            {input.seq} &> {log}
        else
            touch {output}
        fi
        """


# Rename repeatmasker output to match the original sequence names.
rule reformat_repeatmasker_output:
    input:
        script=workflow.source_path("../scripts/reformat_rm.py"),
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
    log:
        join(RM_LOGDIR, "reformat_repeatmasker_output_{sm}_{fname}.log"),
    conda:
        "../envs/py.yaml"
    shell:
        """
        python {input.script} -i {input.rm_out} -of {input.original_fai} -rf {input.renamed_fai} > {output} 2> {log}
        """


# Gather all RM output
def refmt_rm_output(wc):
    _ = checkpoints.split_cens_for_rm.get(**wc).output
    fa_glob_pattern = join(
        expand(rules.split_cens_for_rm.output, sm=wc.sm)[0], "{fname}.fa"
    )
    wcs = glob_wildcards(fa_glob_pattern)
    fnames = wcs.fname
    return expand(rules.reformat_repeatmasker_output.output, sm=wc.sm, fname=fnames)


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
        awk -v OFS="\\t" '{{$1=$1; print}}' {input.rm_out} > {output}
        """


# Merge repeatmasker and convert sep to tab.
# |1259|28.4|7.4|5.3|GM18989_chr1_hap1-0000003:9717731-15372230|8|560|(5653940)|+|Charlie2b|DNA/hAT-Charlie|120|683|(2099)|1|
rule merge_repeatmasker_output:
    input:
        ancient(expand(rules.format_repeatmasker_output.output, sm=SAMPLE_NAMES)),
    output:
        temp(
            join(
                RM_OUTDIR,
                "repeats",
                "all",
                "all_samples_cens.fa.out",
            )
        ),
    shell:
        """
        cat {input} > {output}
        """


# Format repeatmasker reference output.
rule merge_control_repeatmasker_output:
    input:
        # Contains header. Should be first.
        config["repeatmasker"]["ref_repeatmasker_output"],
    output:
        join(
            RM_OUTDIR,
            "repeats",
            "ref",
            "ref_ALR_regions.fa.out",
        ),
    shell:
        """
        cat {input} | awk -v OFS='\\t' '{{$1=$1; print}}' | cut -f 1-15 > {output}
        """


rule format_add_control_repeatmasker_output:
    input:
        ref_rm_output=rules.merge_control_repeatmasker_output.output,
        sample_rm_output=rules.merge_repeatmasker_output.output,
    output:
        temp(
            join(
                RM_OUTDIR,
                "repeats",
                "all",
                "all_samples_and_ref_cens.fa.out",
            )
        ),
    log:
        join(RM_LOGDIR, "format_add_control_repeatmasker_output.log"),
    conda:
        "../envs/tools.yaml"
    shell:
        """
        # Copy file and append reference repeatmasker output.
        cp {input.sample_rm_output} {output} 2> {log}
        {{ awk -v OFS="\\t" '{{$1=$1; print}}' {input.ref_rm_output} | grep "chr" ;}} >> {output} 2> {log}
        """


rule create_og_rm_bed:
    input:
        script=workflow.source_path("../scripts/create_rm_bed.py"),
        rm_out=rules.format_add_control_repeatmasker_output.output,
    output:
        rm_bed=join(
            RM_OUTDIR,
            "bed",
            "{chr}_og",
            "rm.bed",
        ),
    params:
        chr_rgx="{chr}[:_]",
        color_mapping=config["repeatmasker"]["repeat_colors"],
    log:
        join(RM_LOGDIR, "create_og_rm_bed_{chr}.log"),
    conda:
        "../envs/py.yaml"
    shell:
        shell_create_rm_bed


rule modify_og_rm_cenplot_tracks:
    input:
        script=workflow.source_path("../scripts/format_cenplot_toml.py"),
        plot_layout=workflow.source_path("../scripts/cenplot_repeatmasker_plot.toml"),
        infiles=rules.create_og_rm_bed.output,
    output:
        plot_layout=join(
            RM_OUTDIR,
            "plots",
            "{chr}_cens_og.yaml",
        ),
    conda:
        "../envs/py.yaml"
    log:
        join(RM_LOGDIR, "modify_og_rm_cenplot_tracks_{chr}.log"),
    params:
        indir=lambda wc, input: os.path.abspath(os.path.dirname(str(input.infiles[0]))),
    shell:
        """
        python {input.script} \
        -i {input.plot_layout} \
        -o {output.plot_layout} \
        -k indir={params.indir} &> {log}
        """


rule plot_og_rm_bed_by_chr:
    input:
        bed_files=[rules.create_og_rm_bed.output],
        script=workflow.source_path("../scripts/plot_multiple_cen.py"),
        plot_layout=rules.modify_og_rm_cenplot_tracks.output,
    output:
        plots=multiext(
            join(
                RM_OUTDIR,
                "plots",
                "{chr}_cens_og",
            ),
            ".pdf",
            ".png",
        ),
        plot_dir=directory(
            join(
                RM_OUTDIR,
                "plots",
                "{chr}_cens_og",
            )
        ),
    threads: 4
    log:
        join(RM_LOGDIR, "plot_og_rm_bed_{chr}.log"),
    conda:
        "../envs/py.yaml"
    shell:
        shell_plot_multiple_cen


rule repeatmasker_all:
    input:
        expand(rules.plot_og_rm_bed_by_chr.output, chr=CHROMOSOMES),
    default_target: True
