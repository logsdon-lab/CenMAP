
include: "common.smk"
include: "utils.smk"
include: "3-srf.smk"
include: "5-ident_cen_ctgs.smk"
include: "7-fix_cens_w_repeatmasker.smk"


HUMAS_ANNOT_OUTDIR = join(OUTPUT_DIR, "8-humas_annot")
HUMAS_ANNOT_LOGDIR = join(LOG_DIR, "8-humas_annot")
HUMAS_ANNOT_BMKDIR = join(BMK_DIR, "8-humas_annot")


rule extract_cens_for_humas_annot:
    input:
        fa=rules.rename_reort_asm.output.fa,
        bed=ancient(rules.make_complete_cens_bed.output.cen_bed),
    output:
        seq=temp(
            join(
                HUMAS_ANNOT_OUTDIR,
                "temp",
                "{sm}_cens.fa",
            )
        ),
        idx=temp(
            join(
                HUMAS_ANNOT_OUTDIR,
                "temp",
                "{sm}_cens.fa.fai",
            )
        ),
    conda:
        "../envs/tools.yaml"
    log:
        join(HUMAS_ANNOT_LOGDIR, "extract_cens_for_humas_annot_{sm}.log"),
    params:
        bed=lambda wc, input: input.bed,
        added_cmds="| seqkit seq --upper-case",
    shell:
        shell_extract_and_index_fa


checkpoint split_cens_for_humas_annot:
    input:
        fa=rules.extract_cens_for_humas_annot.output.seq,
    output:
        touch(
            join(
                HUMAS_ANNOT_OUTDIR,
                "split_cens_for_humas_annot_{sm}.done",
            )
        ),
    log:
        join(HUMAS_ANNOT_LOGDIR, "split_{sm}_cens_for_humas_annot.log"),
    params:
        split_dir=HUMAS_CENS_SPLIT_DIR,
    conda:
        "../envs/tools.yaml"
    shell:
        # https://gist.github.com/astatham/621901
        """
        mkdir -p {params.split_dir}
        awk '{{
            if (substr($0, 1, 1)==">") {{
                filename=("{params.split_dir}/" substr($0,2) ".fa")
            }}
            print $0 > filename
        }}' {input.fa} 2> {log}
        """


checkpoint split_srf_trf_monomers:
    input:
        beds=expand(rules.merge_slop_region_bed.output, sm=SAMPLE_NAMES),
    output:
        join(HUMAS_ANNOT_OUTDIR, "srf_monomers", "{fname}.fa"),
    conda:
        "../envs/tools.yaml"
    params:
        output_dir=lambda wc, output: dirname(output[0]),
        # Can't use sample wildcard so figure out from fname.
        sample=lambda wc: wc.fname.split("_")[0],
        contig=lambda wc: wc.fname.split("_", 2)[2].split(":", 1)[0],
    log:
        join(HUMAS_ANNOT_LOGDIR, "split_srf_trf_monomers_{fname}.log"),
    shell:
        """
        mkdir -p {params.output_dir}
        {{ awk -v OFS="\\t" '{{
            if ($1 != "{params.contig}") {{ next; }}
            split($4, monomers, ",");
            out_fa="{params.output_dir}/{wildcards.fname}.tmp_fa"
            for (i in monomers) {{
                print ">monomer" >> out_fa
                print monomers[i] >> out_fa
            }};
            print out_fa
        }}' $(find {input.beds} -name "{params.sample}.bed") | \
        xargs -I {{}} \
        bash -c 'seqkit rmdup -s {{}} 2>> {log} | seqkit rename 2>> {log} > {params.output_dir}/{wildcards.fname}.fa && rm -f {{}}' ;}} 2> {log}
        """


# Choose which annotation workflow.
added_as_cfg = {}
if config["humas_annot"]["mode"] == "sd":
    humas_module_smk = "Snakemake-HumAS-SD/workflow/Snakefile"
    humas_env = "Snakemake-HumAS-SD/workflow/envs/env.yaml"
elif config["humas_annot"]["mode"] == "srf-n-trf":
    humas_module_smk = "Snakemake-HumAS-SD/workflow/Snakefile"
    humas_env = "Snakemake-HumAS-SD/workflow/envs/env.yaml"
    added_as_cfg["run_stv"] = False
    added_as_cfg["monomer_dir"] = dirname(rules.split_srf_trf_monomers.output[0])
    # Remove hmm profile if added.
    config["humas_annot"].pop("hmm_profile")
else:
    humas_module_smk = "Snakemake-HumAS-HMMER/workflow/Snakefile"
    humas_env = "Snakemake-HumAS-HMMER/workflow/envs/env.yaml"
    added_as_cfg["mode"] = (
        config["humas_annot"]["mode"]
        if config["humas_annot"]["mode"] == "sf"
        else "hor"
    )


module HumAS_Annot:
    snakefile:
        humas_module_smk
    config:
        {
            **config["humas_annot"],
            "input_dir": HUMAS_CENS_SPLIT_DIR,
            "output_dir": HUMAS_ANNOT_OUTDIR,
            "logs_dir": HUMAS_ANNOT_LOGDIR,
            "benchmarks_dir": HUMAS_ANNOT_BMKDIR,
            **added_as_cfg,
        }


use rule * from HumAS_Annot as cens_*


rule format_filter_srf_trf_annot:
    input:
        rules.cens_convert_to_bed9.output
        if config["humas_annot"]["mode"] == "srf-n-trf"
        else [],
    output:
        join(
            HUMAS_ANNOT_OUTDIR,
            "{fname}",
            "stv_mon.bed",
        ),
    params:
        thr_ident=70.0,
    shell:
        """
        awk -v OFS="\\t" -v THR_IDENT={params.thr_ident} '$5 >= THR_IDENT' {input} > {output}
        """


rule format_monomer_sf_classes:
    input:
        sf_annot_colors=(
            config["plot_hor_stv"]["sf_annot_colors"]
            if config.get("plot_hor_stv")
            else []
        ),
        bed=lambda wc: expand(
            rules.cens_filter_reformat_hmm_tbl_to_bed_w_thr.output,
            fname=wc.fname,
            mode="sf",
        ),
    output:
        join(
            HUMAS_ANNOT_OUTDIR,
            "{fname}_sf_colored.bed",
        ),
    shell:
        """
        awk -v OFS="\\t" '{{
            if (FNR == NR) {{
                mon_sf[$2]=$1;
                mon_color[$2]=$3;
                next
            }};
            mon=$4;
            $4=mon_sf[mon];
            $9=mon_color[mon];
            $5=0;
            print
        }}' {input.sf_annot_colors} {input.bed} > {output}
        """


# https://stackoverflow.com/a/63040288
def humas_annot_sm_outputs(wc):
    _ = checkpoints.split_cens_for_humas_annot.get(**wc).output
    wcs = glob_wildcards(
        join(HUMAS_CENS_SPLIT_DIR, wc.sm + "_{chrom}_{ctg}.fa"),
    )
    fnames = [f"{wc.sm}_{chrom}_{ctg}" for chrom, ctg in zip(wcs.chrom, wcs.ctg)]
    chrs = wcs.chrom

    if config["humas_annot"]["mode"] == "srf-n-trf":
        _ = [
            checkpoints.split_srf_trf_monomers.get(fname=fname).output
            for fname in fnames
        ]
        return expand(rules.format_filter_srf_trf_annot.output, fname=fnames)
    elif config["humas_annot"]["mode"] == "sf":
        return expand(rules.format_monomer_sf_classes.output, fname=fnames)
    else:
        return expand(rules.cens_generate_stv.output, zip, fname=fnames, chr=chrs)


checkpoint run_humas_annot:
    input:
        expand(rules.split_cens_for_humas_annot.output, sm=SAMPLE_NAMES),
        humas_annot_sm_outputs,
    output:
        touch(join(HUMAS_ANNOT_OUTDIR, "humas_annot_{sm}.done")),


# https://stackoverflow.com/a/63040288
def humas_annot_chr_outputs(wc):
    _ = [checkpoints.run_humas_annot.get(sm=sm).output for sm in SAMPLE_NAMES]
    fastas = glob.glob(join(HUMAS_CENS_SPLIT_DIR, "*.fa"))
    sms, chroms, ctgs, coords = [], [], [], []
    for fasta in fastas:
        bname, _ = splitext(basename(fasta))
        fname, coord = bname.rsplit(":", 1)
        sm, chrom, ctg = fname.split("_", 2)
        sms.append(sm)
        chroms.append(chrom)
        ctgs.append(ctg)
        coords.append(coord)

    fnames = []
    # Sort by coords so if multiple chr, chr position in name (chr3-chr21) matches.
    sorted_wcs = sorted(
        zip(sms, chroms, ctgs, coords),
        key=lambda x: (x[1], x[2], x[3]),
        reverse=True,
    )

    # Store index of chrom per contig.
    # In cases of dicentric contigs (chr3-chr21). Don't want to include twice.
    ctg_counter = Counter()
    for sm, chrom, ctg, coord in sorted_wcs:
        chrom_names: list[str] = chrom.replace("rc-", "").split("-")
        if not wc.chr in chrom_names:
            continue
        ctg_id = (sm, chrom, ctg, coord)
        idx = ctg_counter[ctg_id]
        ctg_counter[ctg_id] += 1
        if wc.chr != chrom_names[idx]:
            continue
        fnames.append(f"{sm}_{chrom}_{ctg}:{coord}")

    outputs = config.get("plot_hor_stv", {}).get("ref_stv", [])
    if config["humas_annot"]["mode"] == "srf-n-trf":
        _ = [
            checkpoints.split_srf_trf_monomers.get(fname=fname).output
            for fname in fnames
        ]
        outputs.extend(expand(rules.format_filter_srf_trf_annot.output, fname=fnames))
    elif config["humas_annot"]["mode"] == "sf":
        outputs.extend(
            expand(rules.format_monomer_sf_classes.output, fname=fnames, mode="sf")
        )
    else:
        outputs.extend(expand(rules.cens_generate_stv.output, fname=fnames, chr=wc.chr))

    if not outputs:
        raise RuntimeError(f"No hor/sf outputs for {wc.chr}.")
    return outputs


# Force including conda so --containerize includes.
# Must be done since Snakemake won't know rule metadata until runtime.
rule _force_humas_annot_env_inclusion:
    output:
        plots=touch("conda_humas_annot.done"),
    conda:
        humas_env
    shell:
        "echo ''"


rule humas_annot_all:
    input:
        rules._force_humas_annot_env_inclusion.output if IS_CONTAINERIZE_CMD else [],
        expand(rules.run_humas_annot.output, sm=SAMPLE_NAMES),
    default_target: True
