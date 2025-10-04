# Identify and add sample and chromosome to centromeric regions or just add sample.


include: "common.smk"


if CHROMOSOMES:

    include: "5-add_sm_chr_cen_ctgs.smk"

else:

    include: "5-add_sm_cen_ctgs.smk"
