include: "utils.smk"
include: "common.smk"


rule download_ref_asm:
    output:
        REF_FA,
    log:
        join(LOG_DIR, "1-download_ref", f"get_asm_{REF_NAME}.log"),
    params:
        url=REF_URL,
    shell:
        """
        wget --no-verbose {params.url} -O {output} 2> {log}
        """
