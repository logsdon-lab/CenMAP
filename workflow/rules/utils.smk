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


rule plot_rm_bed:
    input:
        rm_bed="",
    output:
        plots="",
        plot_dir="",
    params:
        indir=lambda wc, input: os.path.dirname(str(input.rm_bed)).replace("/", "\\/"),
        script="workflow/scripts/plot_repeatmasker.py",
        plot_layout="workflow/scripts/repeatmasker_plot.toml",
        tmp_plot_layout=lambda wc: f"/tmp/repeatmasker_plot_{wc.chr}.toml",
    threads: 4
    log:
        "logs/plot_rm_out.log",
    conda:
        "../envs/py.yaml"
    shell:
        """
        sed 's/{{indir}}/{params.indir}/g' {params.plot_layout} > {params.tmp_plot_layout}
        python {params.script} \
        -t {params.tmp_plot_layout} \
        -d {output.plot_dir} \
        --share_xlim \
        -p {threads} \
        -c <(cut -f 1 {input.rm_bed} | uniq) 2> {log}
        rm {params.tmp_plot_layout}
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
