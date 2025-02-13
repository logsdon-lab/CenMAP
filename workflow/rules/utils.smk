rule extract_and_index_fa:
    input:
        fa="",
        bed="",
    output:
        seq="",
        idx="",
    log:
        "logs/extract_and_index_fa.log",
    params:
        bed=lambda wc, input: input.bed,
        added_cmds="",
    conda:
        "../envs/tools.yaml"
    shell:
        """
        seqtk subseq {input.fa} {params.bed} {params.added_cmds} > {output.seq} 2> {log}
        # Check if empty before attempting to index. Always create index file.
        if [ -s {output.seq} ]; then
            samtools faidx {output.seq} &> {log}
        else
            touch {output.idx}
        fi
        """


rule create_rm_bed:
    input:
        script="workflow/scripts/create_rm_bed.py",
        rm_out="",
    output:
        rm_bed="",
    params:
        chr_rgx="",
        color_mapping="",
    log:
        "logs/create_rm_bed.log",
    conda:
        "../envs/py.yaml"
    shell:
        """
        python {input.script} -i <(cat {input.rm_out}) -c {params.chr_rgx} -m {params.color_mapping} > {output.rm_bed} 2> {log}
        """


rule modify_cenplot_tracks:
    input:
        plot_layout="",
        infiles="",
    output:
        plot_layout="",
    run:
        import os, sys, tomllib, yaml

        # infile indicates where all other bedfiles are.
        indir = os.path.abspath(os.path.dirname(str(input.infiles[0])))
        with (
            open(input.plot_layout, "rb") as fh,
            open(output.plot_layout, "wt") as out_fh,
        ):
            settings = tomllib.load(fh)
            new_settings = {"settings": settings["settings"], "tracks": []}

            for trk in settings["tracks"]:
                path = trk.get("path")
                # Some tracks don't have paths.
                if not path:
                    new_settings["tracks"].append(trk)
                    continue
                try:
                    new_path = path.format(indir=indir, **dict(params.items()))
                except KeyError:
                    print(
                        f"Invalid format key in path {path} from cenplot track file, {input.plot_layout}.",
                        file=sys.stderr,
                    )
                    continue

                    # Skip if empty.
                if os.stat(new_path).st_size == 0:
                    continue

                new_trk = trk.copy()
                # Pass params from snakemake
                new_trk["path"] = path.format(indir=indir, **dict(params.items()))
                new_settings["tracks"].append(new_trk)

            out_fh.write(yaml.dump(new_settings))


rule plot_multiple_cen:
    input:
        bed_files=[],
        script="workflow/scripts/plot_multiple_cen.py",
        plot_layout="",
    output:
        plots="",
        plot_dir="",
    threads: 4
    log:
        "logs/plot_multiple_cen.log",
    conda:
        "../envs/py.yaml"
    shell:
        """
        # Then use custom script and cenplot.
        {{ python {input.script} \
        -t {input.plot_layout} \
        -d {output.plot_dir} \
        --share_xlim \
        -p {threads} \
        -c <(cut -f 1 {input.bed_files[0]} | sort | uniq) || true ;}} 2> {log}
        # Allow failure. Possible to have no correct cens.
        touch {output.plots}
        mkdir -p {output.plot_dir}
        """


rule wget:
    output:
        "",
    params:
        url="",
    resources:
        mem=1,
    log:
        "",
    shell:
        """
        wget --no-verbose {params.url} -O {output} 2> {log}
        """


rule reorient_bed:
    input:
        bed="",
        # By default:
        # chrom	st	end	new_chrom	new_st	new_end	chrom_len
        og_coords_key="",
    output:
        "",
    log:
        "",
    conda:
        "../envs/tools.yaml"
    params:
        legend_col_chrom_og="$1",
        legend_col_chrom_new='$4":"$5"-"$6',
        legend_col_chrom_len="$7",
        additional_cols="",
    shell:
        """
        nfs=$(awk 'NR < 2 {{ print NF }}' {input.bed})
        {{ join -1 1 -2 1 \
            <(sort -k1 {input.bed}) \
            <(awk -v OFS="\\t" '{{ print {params.legend_col_chrom_og}, {params.legend_col_chrom_new}, {params.legend_col_chrom_len} }}' {input.og_coords_key}) | \
        awk -v OG_NF="${{nfs}}" -v OFS="\\t" '{{
            new_name=$(OG_NF + 1)
            chrom_len=$(OG_NF + 2)
            is_rc=(new_name ~ "rc-");
            if (is_rc) {{
                print new_name, chrom_len-$3, chrom_len-$2 {params.additional_cols}
            }} else {{
                print new_name, $2, $3 {params.additional_cols}
            }}
        }}';}} > {output} 2> {log}
        """
