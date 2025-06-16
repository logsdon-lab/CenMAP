include: "utils.smk"
include: "common.smk"


use rule wget as download_ref_asm with:
    output:
        REF_FA,
    log:
        join(LOG_DIR, "0-download_ref", f"get_asm_{REF_NAME}.log"),
    params:
        url=REF_URL,
