include: "common.smk"
include: "utils.smk"
include: "5-ident_cen_ctgs.smk"


if RUN_REPEATMASKER:

    include: "7.1-fix_cens_w_repeatmasker.smk"

else:

    include: "7.1-fix_cens_w_kmers.smk"


rule finalize_cens_all:
    input:
        rules.fix_cens_w_repeatmasker_all.input
        if RUN_REPEATMASKER
        else rules.fix_cens_w_kmers_all.input,
    default_target: True
