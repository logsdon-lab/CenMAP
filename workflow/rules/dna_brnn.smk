# https://github.com/lh3/dna-nn/tree/master
rule run_dna_brnn:
    input:
        model=config["dna_brnn"]["model"],
        seqs=rules.extract_cens_oriented_regions.output.seq,
    output:
        repeat_regions=os.path.join(
            config["dna_brnn"]["output_dir"],
            "{sm}",
            "{sm}_centromeric_regions.renamed.{ort}.bed",
        ),
    threads: config["dna_brnn"]["threads"]
    resources:
        mem_mb=30_000,
    log:
        "logs/dna_brnn_{ort}_{sm}.log",
    benchmark:
        "benchmarks/dna_brnn_{ort}_{sm}.tsv"
    # No conda recipe. Use singularity if not installed locally.
    singularity:
        "docker://logsdonlab/dna-nn:latest"
    shell:
        """
        dna-brnn -t {threads} -Ai {input.model} {input.seqs} > {output} 2> {log}
        """


# Run dna-brnn only on the ref centromeres as baseline expectation of ALR repeats.
# (Per chr)
# /net/eichler/vol27/projects/hgsvc/nobackups/analysis/glennis/map_align_t2t_chm13/results/t2t_chm13_v2.0_cens/bed/centromeres/dna-brnn/chm13_cens.trimmed.bed
# dna-brnn on f"chm13.hor_arrays_masked.500kbp.fa"?
use rule run_dna_brnn as run_dna_brnn_ref_cens with:
    input:
        model=config["dna_brnn"]["model"],
        seqs=rules.extract_ref_hor_arrays.output,
    output:
        cens=os.path.join(
            config["dna_brnn"]["output_dir"], REF_NAME, f"{REF_NAME}_cens.trimmed.bed"
        ),
    log:
        f"logs/dna_brnn_{REF_NAME}_cens.log",
    benchmark:
        f"benchmarks/dna_brnn_{REF_NAME}_cens.tsv"


# TODO: Script different from notebook?
# grep "chr1:" chm13_cens.trimmed.bed | \
# sed 's/:/\t/g' | sed 's/-/\t/g' | \
# awk -v OFS="\t" '{print $1, $2+$4, $2+$5, $6, $5-$4}' | awk '$4==2' | awk '$5>1000' > chr1_tmp.fwd.bed
rule filter_dnabrnn_ref_cens_regions:
    input:
        repeats=(
            rules.run_dna_brnn_ref_cens.output
            if config["dna_brnn"].get("ref_alr_file") is None
            else config["dna_brnn"]["ref_alr_file"]
        ),
    output:
        temp(
            os.path.join(
                config["dna_brnn"]["output_dir"], REF_NAME, "{chr}_tmp.fwd.bed"
            )
        ),
    params:
        repeat_type_filter=2,
        repeat_len_thr=1000,
    log:
        "logs/filter_dnabrnn_ref_{chr}_cens_regions.log",
    conda:
        "../env/tools.yaml"
    shell:
        """
        {{ grep "{wildcards.chr}:" {input.repeats} | \
        sed -e 's/:/\\t/g' -e 's/-/\\t/g' | \
        awk -v OFS="\\t" '{{print $1, $2+$4, $2+$5, $6, $5-$4}}' | \
        awk '$4=={params.repeat_type_filter} && $5>{params.repeat_len_thr}';}} > {output} 2> {log}
        """


# ex. >HG00171_chr19_h1tg000004l#1-58112442
# ex. >HG00171_chr16_haplotype1-0000003:4-8430174
# HG00171_chr16_haplotype1|0000003|4|8430174|start|end|type
# (Per chr + sample)
# grep "chr2_" ${sample}.renamed.fwd.bed | \
# sed 's/:/\t/g' | \
# sed 's/-/\t/g' | \
# awk -v OFS="\t" '{print $1"-"$2, $3+$5, $3+$6, $7, $6-$5}' | \
# awk '$4==2' | \
# awk '$5>1000' >> chr2_tmp.fwd.bed
rule filter_dnabrnn_sample_cens_regions:
    input:
        repeats=rules.run_dna_brnn.output,
    output:
        tmp_alr_ctgs=temp(
            os.path.join(
                config["dna_brnn"]["output_dir"],
                "{sm}",
                "{chr}_{sm}_contigs.{ort}.ALR.bed",
            )
        ),
    params:
        repeat_type_filter=2,
        awk_repeat_len_thr_stmt=lambda wc: build_awk_cen_region_length_thr(str(wc.chr)),
    log:
        "logs/filter_dnabrnn_{ort}_{sm}_{chr}_cens_regions.log",
    conda:
        "../env/tools.yaml"
    shell:
        # TODO: This may fail. Replace with regex.
        """
        chr_repeats=$(grep "{wildcards.chr}_" {input.repeats} || true)
        if [ -z "${{chr_repeats}}" ]; then
            # Still make the file even if chr doesn't exist.
            touch {output.tmp_alr_ctgs}
        else
           {{ printf '%s\\n' "${{chr_repeats}}" | \
            sed -e 's/h1tg/h1-tg/g' -e 's/h2tg/h2-tg/g' | \
            sed -e 's/#/\\t/g' -e 's/:/\\t/g' -e 's/-/\\t/g' | \
            awk -v OFS="\\t" '{{
                is_hifiasm=$1 ~ "h1" || $1 ~ "h2"
                if ("{wildcards.ort}" == "rev") {{
                    if (is_hifiasm) {{
                        new_start=$6-$8
                        new_stop=$6-$7
                    }} else {{
                        new_start=$4-$6
                        new_stop=$4-$5
                    }}
                }} else {{
                    if (is_hifiasm) {{
                        new_start=$5+$7
                        new_stop=$5+$8
                    }} else {{
                        new_start=$3+$5
                        new_stop=$3+$6
                    }}
                }}
                if (is_hifiasm) {{
                    contig_name=$1""$2"#"$3"-"$4
                    print contig_name, new_start, new_stop, $9, $8-$7
                }} else {{
                    contig_name=$1"-"$2
                    print contig_name, new_start, new_stop, $7, $6-$5
                }}
            }}' | \
            awk '$4=={params.repeat_type_filter} && {params.awk_repeat_len_thr_stmt}';}} > {output.tmp_alr_ctgs} 2> {log}
        fi
        """


# Calculate the start and end *-terms
# ex. chr2
# awk -v OFS="\t" '{print $1, $2-(*467987), $3+(*522450), $3-$2}' | \
rule get_dnabrnn_ref_cens_pos:
    input:
        script="workflow/scripts/get_cen_pos.py",
        repeats=(
            rules.run_dna_brnn_ref_cens.output
            if config["dna_brnn"].get("ref_alr_file") is None
            else config["dna_brnn"]["ref_alr_file"]
        ),
    output:
        os.path.join(config["dna_brnn"]["output_dir"], REF_NAME, "{chr}_cens_data.json"),
    log:
        "logs/get_dnabrnn_ref_{chr}_cens_data.log",
    params:
        repeat_type_filter=2,
        repeat_len_thr=1000,
    conda:
        "../env/py.yaml"
    shell:
        # Read from stdin filtered ref alr repeats. Done separately from above because of formatting.
        """
        python {input.script} -o {output} -i &> {log} \
        <(grep "{wildcards.chr}:" {input.repeats} | \
        awk -v OFS="\\t" \
        '{{len=$3-$2; rp_typ=$4; if (rp_typ=={params.repeat_type_filter} && len>{params.repeat_len_thr}) print}}')
        """


def alr_region_threshold(wc):
    if wc.chr == "chrY":
        return 300_000
    else:
        return 1_000_000


# TODO: Annotate
# /net/eichler/vol28/home/glogsdon/utilities/bedminmax.py (modified bedminmax) \
# -i chr2_tmp.fwd.bed | \
# awk -v OFS="\t" '{print $1, $2, $3, $3-$2}' | \
# awk -v OFS="\t" '{print $1, $2-467987, $3+522450, $3-$2}' | \
# awk '$4>1000000' | \
# awk -v OFS="\t" '$2<0 {$2=0}1' > chr2_contigs.fwd.repeat.bed
rule aggregate_dnabrnn_alr_regions_by_chr:
    input:
        script="workflow/scripts/filter_cen_ctgs.py",
        cen_pos=rules.get_dnabrnn_ref_cens_pos.output,
        added_ref_cens=lambda wc: (
            rules.filter_dnabrnn_ref_cens_regions.output if wc.ort == "fwd" else []
        ),
        sample_cens=lambda wc: expand(
            rules.filter_dnabrnn_sample_cens_regions.output,
            sm=SAMPLE_NAMES,
            ort=[wc.ort],
            chr=[wc.chr],
        ),
    output:
        os.path.join(
            config["dna_brnn"]["output_dir"],
            "{chr}_contigs.{ort}.ALR.bed",
        ),
    params:
        repeat_len_thr=alr_region_threshold,
        input_cols=" ".join(
            [
                "chr",
                "start",
                "end",
                "repeat_type",
                "repeat_length",
            ]
        ),
        output_cols=" ".join(["chr", "start", "end", "repeat_type"]),
        grp_cols=" ".join(["chr", "repeat_type"]),
        sort_cols=" ".join(["chr", "start"]),
        # dna-brnn output may not contain repeats from a chr.
        allow_empty="--allow_empty",
    log:
        "logs/aggregate_dnabrnn_alr_regions_by_{chr}_{ort}.log",
    conda:
        "../env/py.yaml"
    # Aggregate and bedminmax all.
    # Select cols and calculate length.
    # Calculate vals.
    # Take only repeats greater than some value.
    # Take abs value.
    shell:
        """
        {{ python {input.script} bedminmax \
            -i {input.sample_cens} {input.added_ref_cens} \
            -ci {params.input_cols} \
            -co {params.output_cols} \
            -g {params.grp_cols} \
            -s {params.sort_cols} {params.allow_empty} | \
        awk -v OFS="\\t" '{{print $1, $2, $3, $3-$2}}' | \
        awk -v START_POS=$(jq .start {input.cen_pos}) \
            -v END_POS=$(jq .end {input.cen_pos}) \
            -v OFS="\\t" '{{print $1, $2-START_POS, $3+END_POS, $3-$2}}' | \
        awk '$4>{params.repeat_len_thr}' | \
        awk -v OFS="\\t" '$2<0 {{$2=0}}1';}} > {output} 2> {log}
        """


rule dna_brnn_all:
    input:
        expand(rules.run_dna_brnn.output, sm=SAMPLE_NAMES, ort=ORIENTATION),
        expand(
            rules.get_dnabrnn_ref_cens_pos.output,
            chr=CHROMOSOMES,
        ),
        expand(
            rules.aggregate_dnabrnn_alr_regions_by_chr.output,
            ort=ORIENTATION,
            chr=CHROMOSOMES,
        ),
