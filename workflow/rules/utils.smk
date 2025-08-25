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
    python {params.script} -i <(cat {input.rm_out}) -c {params.chr_rgx} -m {params.color_mapping} {params.to_abs} > {output.rm_bed} 2> {log}
"""


shell_plot_multiple_cen = """
    # Then use custom script and cenplot.
    {{ python {params.script} \
    -t {input.plot_layout} \
    -d {output.plot_dir} \
    --share_xlim \
    -p {threads} \
    -c $(cut -f 1 {input.bed_files[0]} | sort | uniq) || true ;}} 2> {log}
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
