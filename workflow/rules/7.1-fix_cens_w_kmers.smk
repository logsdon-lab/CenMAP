
include: "5-ident_cen_ctgs.smk"


FIX_KMERS_OUTDIR = join(OUTPUT_DIR, "7.1-fix_cens_w_kmers")
FIX_KMERS_LOGDIR = join(LOG_DIR, "7.1-fix_cens_w_kmers")
FIX_KMERS_BMKDIR = join(BMK_DIR, "7.1-fix_cens_w_kmers")


rule calculate_entropy:
    input:
        fa=rules.extract_cens_regions.output.seq,
    output:
        # BED9
        entropy_dir=directory(
            join(
                FIX_KMERS_OUTDIR,
                "entropy",
                "{sm}",
            )
        ),
    params:
        window=100_000,
        kmer_size=100,
        omit_plot="--omit_plot",
    threads: config["ident_cen_ctgs"]["threads_srf"]
    log:
        join(FIX_KMERS_LOGDIR, "calculate_entropy_{sm}.log"),
    conda:
        "../envs/py.yaml"
    shell:
        """
        censtats kmer-entropy \
        -i {input.fa} \
        -w {params.window} \
        -k {params.kmer_size} \
        -o {output.entropy_dir} \
        -c {threads} \
        {params.omit_plot} 2> {log}
        """


# Filter valid cens based on entropy
rule filter_entropy_bed:
    input:
        entropy_dir=rules.calculate_entropy.output,
        # Renamed file
        putative_alr_bed=rules.make_srf_putative_alr_regions.output,
    output:
        # (chrom, st, end, old_chrom)
        bed=temp(
            join(
                FIX_KMERS_OUTDIR,
                "entropy",
                "interm",
                "{sm}.bed",
            )
        ),
    log:
        join(FIX_KMERS_LOGDIR, "filter_entropy_bed_{sm}.log"),
    params:
        script=workflow.source_path("../scripts/filter_entropy_bed_kmers.py"),
        bp_merge=config["ident_cen_ctgs"]["bp_merge"],
    conda:
        "../envs/py.yaml"
    shell:
        """
        python {params.script} \
        -i {input.entropy_dir}/*.bed \
        -b {input.putative_alr_bed} \
        -d {params.bp_merge} > {output} 2> {log}
        """


# Merge complete centromere coordinates. Prior to running NucFlag.
rule make_complete_cens_bed:
    input:
        # (sm_ctg, st, end, old_sm_ctg)
        bed=rules.filter_entropy_bed.output,
        # (ctg, sm_ctg, ctg_len)
        rename_key=rules.create_rename_key.output,
        idx=rules.create_final_asm.output.idx,
    output:
        # (sm_ctg, st, end, ctg, ctg_len)
        cen_bed=join(
            OUTPUT_DIR,
            "final",
            "bed",
            "{sm}_complete_cens.bed",
        ),
        # (old_sm_ctg, sm_ctg)
        rename_key=temp(
            join(
                OUTPUT_DIR,
                "final",
                "bed",
                "{sm}_complete_cens_rename_key.tsv",
            )
        ),
    conda:
        "../envs/tools.yaml"
    shell:
        """
        sort -k1,1 -k 2,2n {input.bed}  | \
        join -1 1 -2 2 - <(sort -k2,2 {input.rename_key}) | \
        awk -v FS=" " -v OFS="\\t" '{{
            # Adjust for seqtk 1-index
            $2=($2 == 0) ? 1 : $2
            new_name=$1":"$2+1"-"$3
            print $4, new_name >> "{output.rename_key}"
            print $1, $2, $3, $5, $6
        }}' > {output.cen_bed}
        """


rule fix_cens_w_kmers_all:
    input:
        expand(rules.make_complete_cens_bed.output, sm=SAMPLE_NAMES),
    default_target: True
