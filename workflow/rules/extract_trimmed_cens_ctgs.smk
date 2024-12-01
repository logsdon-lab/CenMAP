# Extract trimmed cens from dna-brnn.
include: "common.smk"


use rule extract_and_index_fa as extract_alr_region_sample_by_chr with:
    input:
        fa=os.path.join(config["concat_asm"]["output_dir"], "{sm}-asm-comb-dedup.fa"),
        bed=os.path.join(
            config["dna_brnn"]["output_dir"],
            "bed",
            "{sm}_contigs.ALR.bed",
        ),
    output:
        seq=temp(
            os.path.join(
                config["ident_cen_ctgs"]["output_dir"],
                "seq",
                "interm",
                "{sm}_contigs.ALR.fa",
            )
        ),
        idx=temp(
            os.path.join(
                config["ident_cen_ctgs"]["output_dir"],
                "seq",
                "interm",
                "{sm}_contigs.ALR.fa.fai",
            )
        ),
    log:
        "logs/extract_new_cens_ctgs/extract_alr_region_{sm}.log",
