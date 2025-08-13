
include: "common.smk"
include: "utils.smk"
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
                "seq",
                "interm",
                "{sm}_cens.fa",
            )
        ),
        idx=temp(
            join(
                HUMAS_ANNOT_OUTDIR,
                "seq",
                "interm",
                "{sm}_cens.fa.fai",
            )
        ),
    conda:
        "../envs/tools.yaml"
    log:
        join(HUMAS_ANNOT_LOGDIR, "extract_cens_for_humas_annot_{sm}.log"),
    params:
        # Seqtk outputs 1-based coords which causes name issues.
        bed=lambda wc, input: f"""<(awk -v OFS="\\t" '{{new_len=$2-1; print $1, (new_len < 0) ? 0 : new_len, $3}}' {input.bed})""",
        added_cmds="",
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
        expand(rules.merge_slop_region_bed.output, sm=SAMPLE_NAMES),
    output:
        join(HUMAS_ANNOT_OUTDIR, "srf_monomers", "{fname}.fa"),
    conda:
        "../envs/tools.yaml"
    params:
        output_dir=lambda wc, output: dirname(output[0]),
        # Can't use sample wildcard so figure out from fname.
        sample=lambda wc: wc.fname.split("_")[0],
        contig=lambda wc: wc.fname.split("_", 2)[2].split(":", 1)[0],
    shell:
        """
        mkdir -p {params.output_dir}
        awk -v OFS="\\t" '{{
            if ($1 != "{params.contig}") {{ next; }}
            split($4, monomers, ",");
            out_fa="{params.output_dir}/{wildcards.fname}.tmp_fa"
            for (i in monomers) {{
                print ">monomer" >> out_fa
                print monomers[i] >> out_fa
            }};
            print out_fa
        }}' $(find {input} -name "{params.sample}.bed") | \
        xargs -I {{}} \
        bash -c 'seqkit rmdup -s {{}} | seqkit rename > {params.output_dir}/{wildcards.fname}.fa && rm -f {{}}'
        """


added_as_cfg = {}
if (
    config["humas_annot"]["mode"] == "sd"
    or config["humas_annot"]["mode"] == "srf-n-trf"
):
    humas_module_smk = "Snakemake-HumAS-SD/workflow/Snakefile"
    humas_env = "Snakemake-HumAS-SD/workflow/envs/env.yaml"
    if config["humas_annot"]["mode"] == "srf-n-trf":
        added_as_cfg["run_stv"] = False
        added_as_cfg["monomer_dir"] = dirname(rules.split_srf_trf_monomers.output[0])
        # Remove hmm profile if added.
        config["humas_annot"].pop("hmm_profile")
else:
    humas_module_smk = "Snakemake-HumAS-HMMER/workflow/Snakefile"
    humas_env = "Snakemake-HumAS-HMMER/workflow/envs/env.yaml"


module HumAS_Annot:
    snakefile:
        humas_module_smk
    config:
        {
            **config["humas_annot"],
            # if hmmer, otherwise no effect.
            "mode": "hor",
            "input_dir": HUMAS_CENS_SPLIT_DIR,
            "output_dir": HUMAS_ANNOT_OUTDIR,
            "logs_dir": HUMAS_ANNOT_LOGDIR,
            "benchmarks_dir": HUMAS_ANNOT_BMKDIR,
            **added_as_cfg,
        }


use rule * from HumAS_Annot as cens_*


rule format_filter_srf_trf_annot:
    input:
        rules.cens_convert_to_bed9.output,
    output:
        join(
            HUMAS_ANNOT_OUTDIR,
            "{fname}",
            "stv_mon.bed",
        ),
    params:
        script=workflow.source_path("../scripts/merge_srf_mons.py"),
        thr_ident=70.0,
    shell:
        """
        python {params.script} -i {input} -t {params.thr_ident} > {output}
        """


# https://stackoverflow.com/a/63040288
def humas_annot_sm_outputs(wc):
    _ = checkpoints.split_cens_for_humas_annot.get(**wc).output
    wcs = glob_wildcards(
        join(HUMAS_CENS_SPLIT_DIR, wc.sm + "_{chrom}_{ctg}.fa"),
    )
    fnames = [f"{wc.sm}_{chrom}_{ctg}" for chrom, ctg in zip(wcs.chrom, wcs.ctg)]
    chrs = wcs.chrom
    _ = [checkpoints.split_srf_trf_monomers.get(fname=fname).output for fname in fnames]

    if config["humas_annot"]["mode"] == "srf-n-trf":
        return expand(rules.format_filter_srf_trf_annot.output, fname=fnames)
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
    wcs = glob_wildcards(
        join(HUMAS_CENS_SPLIT_DIR, "{sm}_{chr}_{ctg}:{coords}.fa"),
    )
    fnames = []
    # Sort by coords so if multiple chr, chr position in name (chr3-chr21) matches.
    sorted_wcs = sorted(
        zip(wcs.sm, wcs.chr, wcs.ctg, wcs.coords),
        key=lambda x: (x[1], x[2], x[3]),
        reverse=True,
    )

    # Store index of chrom per contig.
    # In cases of dicentric contigs.
    ctg_counter = Counter()
    for sm, chrom, ctg, coords in sorted_wcs:
        chroms: list[str] = chrom.replace("rc-", "").split("-")
        if not wc.chr in chroms:
            continue
        ctg_id = (sm, chrom, ctg)
        idx = ctg_counter[ctg_id]
        ctg_counter[ctg_id] += 1
        if wc.chr != chroms[idx]:
            continue
        fnames.append(f"{sm}_{chrom}_{ctg}:{coords}")

    return (
        [
            *config.get("plot_hor_stv", {}).get("ref_stv", []),
            *expand(rules.cens_generate_stv.output, fname=fnames, chr=wc.chr),
        ],
    )


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
