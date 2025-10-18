
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
        fa=rules.create_final_asm.output.fa,
        bed=ancient(rules.make_complete_cens_bed.output.cen_bed),
    output:
        seq=pipe(
            join(
                HUMAS_ANNOT_OUTDIR,
                "seq",
                "{sm}_cens.fa",
            )
        ),
        idx=pipe(
            join(
                HUMAS_ANNOT_OUTDIR,
                "seq",
                "{sm}_cens.fa.fai",
            )
        ),
    conda:
        "../envs/tools.yaml"
    log:
        join(HUMAS_ANNOT_LOGDIR, "extract_cens_for_humas_annot_{sm}.log"),
    params:
        bed=lambda wc, input: input.bed,
        added_cmds=f"| {cmd_filter_fa_chrom("seqkit seq --upper-case")}",
    shell:
        shell_extract_and_index_fa


checkpoint split_cens_for_humas_annot:
    input:
        fa=rules.extract_cens_for_humas_annot.output.seq,
        idx=rules.extract_cens_for_humas_annot.output.idx,
    output:
        fofn=join(
            HUMAS_ANNOT_OUTDIR,
            "split_cens_for_humas_annot_{sm}.fofn",
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
                print filename >> "{output.fofn}"
            }}
            print $0 > filename
        }}' {input.fa} 2> {log}
        """


# Choose which annotation workflow.
added_as_cfg = {}
if config["humas_annot"]["mode"] == "sd":
    humas_module_smk = "Snakemake-HumAS-SD/workflow/Snakefile"
    humas_env = "Snakemake-HumAS-SD/workflow/envs/env.yaml"
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


rule format_srf_trf_annot:
    input:
        srf_bed=(
            rules.srf_paf2bed.output
            if config["humas_annot"]["mode"] == "srf-n-trf"
            else []
        ),
        rename_key=rules.create_rename_key.output,
    output:
        # BED9:
        # new_name, st, end, srf_motif, srf_weighted_identity, ".", st, end, "0,0,0"
        join(
            HUMAS_ANNOT_OUTDIR,
            "{sm}",
            "stv_hor.bed",
        ),
    log:
        join(HUMAS_ANNOT_LOGDIR, "format_srf_trf_annot_{sm}.log"),
    conda:
        "../envs/tools.yaml"
    shell:
        """
        {{ cut -f1-5 {input.srf_bed} |
        sort -k1,1 | \
        join - <(sort -k1,1 {input.rename_key}) | \
        awk -v OFS="\\t" '{{
            new_name=$6;
            ctg_len=$7;
            if (new_name ~ "rc-chr") {{
                st=ctg_len-$3;
                end=ctg_len-$2;
            }} else {{
                st=$2;
                end=$3;
            }}
            # prefix#circ126-3075
            match($4, "-([0-9]+)$", motif_lengths);
            print new_name, st, end, "circ-"motif_lengths[1], $5, ".", st, end, "0,0,0"
        }}' ;}} > {output} 2> {log}
        """


rule filter_srf_trf_annot:
    input:
        rules.format_srf_trf_annot.output,
    output:
        join(
            HUMAS_ANNOT_OUTDIR,
            "{sm}",
            "{fname}.bed",
        ),
    log:
        join(HUMAS_ANNOT_LOGDIR, "filter_srf_trf_annot_{sm}_{fname}.log"),
    conda:
        "../envs/tools.yaml"
    params:
        coords_tsv=lambda wc: f"<(printf '{chrom_coord_to_tsv(wc.fname)}')",
    shell:
        """
        {{ bedtools intersect -a {input} -b {params.coords_tsv} | \
        awk -v OFS="\\t" '{{
            $1="{wildcards.fname}";
            # Make relative coordinates
            match($1, ":(.+)-", starts);
            $2=$2-starts[1];
            $3=$3-starts[1];
            $7=$7-starts[1];
            $8=$8-starts[1];
            print
        }}' ;}} > {output} 2> {log}
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
    conda:
        "../envs/tools.yaml"
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
    if CHROMOSOMES:
        wcs = glob_wildcards(
            join(HUMAS_CENS_SPLIT_DIR, wc.sm + "_{chrom}_{ctg}.fa"),
        )
        fnames = [f"{wc.sm}_{chrom}_{ctg}" for chrom, ctg in zip(wcs.chrom, wcs.ctg)]
        chrs = wcs.chrom
    else:
        wcs = glob_wildcards(join(HUMAS_CENS_SPLIT_DIR, wc.sm + "_{ctg}.fa"))
        fnames = [f"{wc.sm}_{ctg}" for ctg in wcs.ctg]
        chrs = []

    if config["humas_annot"]["mode"] == "srf-n-trf":
        return expand(rules.filter_srf_trf_annot.output, sm=wc.sm, fname=fnames)
    elif config["humas_annot"]["mode"] == "sf":
        return expand(rules.format_monomer_sf_classes.output, fname=fnames)
    else:
        if not chrs and config["humas_annot"]["mode"] == "sd":
            return []
        return expand(rules.cens_generate_stv.output, zip, fname=fnames, chr=chrs)


checkpoint run_humas_annot:
    input:
        (
            expand(rules.split_cens_for_humas_annot.output, sm=SAMPLE_NAMES)
            if config["humas_annot"]["mode"] != "srf-n-trf"
            else []
        ),
        humas_annot_sm_outputs,
    output:
        touch(join(HUMAS_ANNOT_OUTDIR, "humas_annot_{sm}.done")),


# https://stackoverflow.com/a/63040288
def humas_annot_chr_outputs(wc):
    _ = [checkpoints.run_humas_annot.get(sm=sm).output for sm in SAMPLE_NAMES]
    fastas = glob.glob(join(HUMAS_CENS_SPLIT_DIR, "*.fa"))
    fnames = get_valid_fnames(fastas, filter_chrom=wc.chr if wc.chr != "all" else None)

    outputs = config.get("plot_hor_stv", {}).get("ref_stv", [])
    if config["humas_annot"]["mode"] == "srf-n-trf":
        outputs.extend(
            expand(
                rules.filter_srf_trf_annot.output,
                zip,
                sm=[fname.sm for fname in fnames],
                fname=fnames,
            )
        )
    elif config["humas_annot"]["mode"] == "sf":
        outputs.extend(
            expand(rules.format_monomer_sf_classes.output, fname=fnames, mode="sf")
        )
    else:
        assert (
            wc.chr != "all" and config["humas_annot"]["mode"] == "sd"
        ), "Chromosomes needed for HumAS-SD."
        outputs.extend(expand(rules.cens_generate_stv.output, fname=fnames, chr=wc.chr))

    return outputs


# Force including conda so --containerize includes.
# Must be done since Snakemake won't know rule metadata until runtime.
rule _force_humas_sd_env_inclusion:
    output:
        plots=touch("conda_humas_sd.done"),
    conda:
        "Snakemake-HumAS-SD/workflow/envs/env.yaml"
    shell:
        "echo ''"


rule _force_humas_hmmer_env_inclusion:
    output:
        plots=touch("conda_humas_hmmer.done"),
    conda:
        "Snakemake-HumAS-HMMER/workflow/envs/env.yaml"
    shell:
        "echo ''"


rule humas_annot_all:
    input:
        rules._force_humas_sd_env_inclusion.output if IS_CONTAINERIZE_CMD else [],
        rules._force_humas_hmmer_env_inclusion.output if IS_CONTAINERIZE_CMD else [],
        expand(rules.run_humas_annot.output, sm=SAMPLE_NAMES),
    default_target: True
