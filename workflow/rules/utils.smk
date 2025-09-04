# Rule inheritance breaks in modules
# https://github.com/snakemake/snakemake/issues/1757

params_shell_extract_and_index_fa = {
    "bed": lambda wc, input: input.bed,
    "added_cmds": "",
}


shell_extract_and_index_fa = """
    seqtk subseq {input.fa} {params.bed} {params.added_cmds} > {output.seq} 2> {log}
    # Check if empty before attempting to index. Always create index file.
    if [ -s {output.seq} ]; then
        samtools faidx {output.seq} &> {log}
    else
        touch {output.idx}
    fi
"""


shell_create_rm_bed = """
    python {params.script} \
    -i <(awk -v OFS="\\t" '{{$1=$1; print}}' {input.rm_out} | cut -f1-15) \
    -c {params.chr_rgx} \
    -m {params.color_mapping} {params.to_abs} > {output.rm_bed} 2> {log}
"""


shell_plot_multiple_cen = """
    {{ python {params.script} \
    -i '{params.json_file_str}' \
    -t {params.plot_layout} \
    -o {params.output_prefix} \
    {params.options} || true ;}} 2> {log}
    touch {output.plots}
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
