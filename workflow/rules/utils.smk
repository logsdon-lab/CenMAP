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
        added_cmds="",
    conda:
        "../env/tools.yaml"
    shell:
        """
        seqtk subseq {input.fa} {input.bed} {params.added_cmds} > {output.seq} 2> {log}
        # Check if empty before attempting to index. Always create index file.
        if [ -s {output.seq} ]; then
            samtools faidx {output.seq} &> {log}
        else
            touch {output.idx}
        fi
        """
