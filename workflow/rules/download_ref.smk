include: "utils.smk"
include: "common.smk"


use rule wget as download_ref_asm with:
    output:
        REF_FA,
    log:
        f"logs/download_ref/get_asm_{REF_NAME}.log",
    params:
        url="https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz",
