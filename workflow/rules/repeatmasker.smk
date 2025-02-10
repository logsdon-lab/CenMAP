
include: "common.smk"
include: "utils.smk"


checkpoint split_cens_for_rm:
    input:
        fa=os.path.join(
            config["ident_cen_ctgs"]["output_dir"],
            "seq",
            "interm",
            "{sm}_contigs.ALR.fa",
        ),
        # [old_name, new_name, coords, sm, chr, is_reversed]
        rename_key=os.path.join(
            config["ident_cen_ctgs"]["output_dir"],
            "bed",
            "interm",
            "{sm}_renamed_cens.tsv",
        ),
    output:
        directory(
            os.path.join(config["repeatmasker"]["output_dir"], "seq", "interm", "{sm}")
        ),
    log:
        "logs/humas_sd/split_{sm}_cens_for_humas_sd.log",
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
                kv[$1]=$2;
                next;
            }}
            if (substr($0, 1, 1)==">") {{
                ctg_name=substr($0,2)
                split(ctg_name, ctg_name_parts, ":")
                new_ctg_name=kv[ctg_name_parts[1]]":"ctg_name_parts[2]
                filename=("{params.split_dir}/" new_ctg_name ".fa")
                print ">"new_ctg_name > filename
            }} else {{
                print $0 > filename
            }}
        }}' {input.rename_key} {input.fa} 2> {log}
        """


# RepeatMasker has a limit of 50 characters for sequence names.
# While I was able to create a fork of RepeatMasker that allowed for longer sequence names, I still ran into issues.
# Bumping to v4.1.5 and also adding the increased limit fixed the issue but resulted in a 5-10x increase in runtime.
# So I downgraded the repeatmasker70 image to v4.1.0.
rule rename_for_repeatmasker:
    input:
        fa=os.path.join(rules.split_cens_for_rm.output[0], "{fname}.fa"),
    output:
        original_fa_idx=temp(
            os.path.join(rules.split_cens_for_rm.output[0], "{fname}.fa.fai"),
        ),
        # Use different dir to avoid greedy wildcard recursively running rule.
        # TODO: rules.split_cens_for_rm.output[0]
        renamed_fa=temp(
            os.path.join(
                config["repeatmasker"]["output_dir"],
                "seq",
                "{sm}_renamed",
                "{fname}.fa",
            )
        ),
        renamed_fa_idx=temp(
            os.path.join(
                config["repeatmasker"]["output_dir"],
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
        "logs/repeatmasker/rename_for_repeatmasker_{sm}_{fname}.log",
    shell:
        """
        samtools faidx {input.fa} 2> {log}
        seqtk rename {input.fa} {params.prefix} > {output.renamed_fa} 2>> {log}
        if [ -s {output.renamed_fa} ]; then
            samtools faidx {output.renamed_fa} 2>> {log}
        else
            touch {output.renamed_fa_idx}
        fi
        """


rule run_repeatmasker:
    input:
        seq=rules.rename_for_repeatmasker.output.renamed_fa,
    output:
        temp(
            os.path.join(
                config["repeatmasker"]["output_dir"],
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
        "logs/repeatmasker/repeatmasker_{sm}_{fname}.log",
    # Retry in case of .RepeatMaskerCache failure.
    retries: 2
    benchmark:
        "benchmarks/repeatmasker/repeatmasker_{sm}_{fname}.tsv"
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
        script="workflow/scripts/reformat_rm.py",
        rm_out=rules.run_repeatmasker.output,
        original_fai=rules.rename_for_repeatmasker.output.original_fa_idx,
        renamed_fai=rules.rename_for_repeatmasker.output.renamed_fa_idx,
    output:
        os.path.join(
            config["repeatmasker"]["output_dir"],
            "repeats",
            "{sm}",
            "{fname}.fa.out",
        ),
    log:
        "logs/repeatmasker/reformat_repeatmasker_output_{sm}_{fname}.log",
    conda:
        "../envs/py.yaml"
    shell:
        """
        python {input.script} -i {input.rm_out} -of {input.original_fai} -rf {input.renamed_fai} > {output} 2> {log}
        """


# Gather all RM output
def refmt_rm_output(wc):
    _ = checkpoints.split_cens_for_rm.get(**wc).output
    fa_glob_pattern = os.path.join(
        expand(rules.split_cens_for_rm.output, sm=wc.sm)[0], "{fname}.fa"
    )
    wcs = glob_wildcards(fa_glob_pattern)
    fnames = wcs.fname
    assert (
        len(fnames) != 0
    ), f"No fasta files found for repeatmasker in {fa_glob_pattern}"
    return expand(rules.reformat_repeatmasker_output.output, sm=wc.sm, fname=fnames)


rule format_repeatmasker_output:
    input:
        rm_out=refmt_rm_output,
        # Force snakemake to not evaluate chkpt function until all dirs created.
        rm_fa_dirs=expand(rules.split_cens_for_rm.output, sm=SAMPLE_NAMES),
    output:
        os.path.join(
            config["repeatmasker"]["output_dir"],
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
        expand(rules.format_repeatmasker_output.output, sm=SAMPLE_NAMES),
    output:
        temp(
            os.path.join(
                config["repeatmasker"]["output_dir"],
                "repeats",
                "all",
                "all_samples_cens.fa.out",
            )
        ),
    params:
        added_cmd="",
    shell:
        """
        cat {input} {params.added_cmd} > {output}
        """


# Format repeatmasker reference output.
# TODO: This should be merged.
use rule merge_repeatmasker_output as merge_control_repeatmasker_output with:
    input:
        # Contains header. Should be first.
        config["repeatmasker"]["ref_repeatmasker_chrY_output"],
        config["repeatmasker"]["ref_repeatmasker_output"],
    params:
        added_cmd="| awk -v OFS='\\t' '{{$1=$1; print}}' | cut -f 1-15",
    output:
        os.path.join(
            config["repeatmasker"]["output_dir"],
            "repeats",
            "ref",
            "ref_ALR_regions.fa.out",
        ),


# Convert reference relative coordinates to absolute.
rule convert_coords_control_repeatmasker_output:
    input:
        rules.merge_control_repeatmasker_output.output,
    output:
        os.path.join(
            config["repeatmasker"]["output_dir"],
            "repeats",
            "ref",
            "ref_ALR_regions.fa.abs.out",
        ),
    shell:
        """
        awk -v OFS="\\t" '{{
            match($5, ":(.+)-", starts);
            $6=$6+starts[1];
            $7=$7+starts[1];
            print
        }}' {input} > {output}
        """


rule format_add_control_repeatmasker_output:
    input:
        ref_rm_output=rules.convert_coords_control_repeatmasker_output.output,
        sample_rm_output=rules.merge_repeatmasker_output.output,
    output:
        temp(
            os.path.join(
                config["repeatmasker"]["output_dir"],
                "repeats",
                "all",
                "all_samples_and_ref_cens.fa.out",
            )
        ),
    log:
        "logs/repeatmasker/format_add_control_repeatmasker_output.log",
    conda:
        "../envs/tools.yaml"
    shell:
        """
        # Copy file and append reference repeatmasker output.
        cp {input.sample_rm_output} {output} 2> {log}
        {{ awk -v OFS="\\t" '{{$1=$1; print}}' {input.ref_rm_output} | grep "chr" ;}} >> {output} 2> {log}
        """


use rule create_rm_bed as create_og_rm_bed with:
    input:
        script="workflow/scripts/create_rm_bed.py",
        rm_out=rules.format_add_control_repeatmasker_output.output,
    output:
        rm_bed=os.path.join(
            config["repeatmasker"]["output_dir"],
            "bed",
            "{chr}_og",
            "rm.bed",
        ),
    params:
        chr_rgx="{chr}[:_]",
        color_mapping=config["repeatmasker"]["repeat_colors"],
    log:
        "logs/repeatmasker/create_og_rm_bed_{chr}.log",


rule modify_rm_cenplot_tracks:
    input:
        plot_layout="workflow/scripts/cenplot_repeatmasker_plot.toml",
        infile=rules.create_og_rm_bed.output,
    output:
        plot_layout=os.path.join(
            config["repeatmasker"]["output_dir"],
            "plots",
            "{chr}_cens_{typ}.yaml",
        ),
    params:
        indir=lambda wc, input: os.path.dirname(str(input.infile)),
    run:
        import tomllib, yaml

        with (
            open(input.plot_layout, "rb") as fh,
            open(output.plot_layout, "wt") as out_fh,
        ):
            settings = tomllib.load(fh)
            for trk in settings["tracks"]:
                if not "path" in trk:
                    continue
                trk["path"] = trk["path"].format(indir=params.indir)
            out_fh.write(yaml.dump(settings))


use rule plot_multiple_cen as plot_og_rm_bed_by_chr with:
    input:
        bed_files=[rules.create_og_rm_bed.output],
        script="workflow/scripts/plot_multiple_cen.py",
        plot_layout=expand(rules.modify_rm_cenplot_tracks.output, chr="{chr}", typ="og"),
    output:
        plots=multiext(
            os.path.join(
                config["repeatmasker"]["output_dir"],
                "plots",
                "{chr}_cens_og",
            ),
            ".pdf",
            ".png",
        ),
        plot_dir=directory(
            os.path.join(
                config["repeatmasker"]["output_dir"],
                "plots",
                "{chr}_cens_og",
            )
        ),
    log:
        "logs/repeatmasker/plot_og_rm_bed_{chr}.log",


rule repeatmasker_all:
    input:
        expand(rules.format_repeatmasker_output.output, sm=SAMPLE_NAMES),
        expand(rules.plot_og_rm_bed_by_chr.output, chr=CHROMOSOMES),
    default_target: True
