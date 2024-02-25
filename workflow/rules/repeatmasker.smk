
rule get_valid_regions_for_rm:
    input:
        bed=get_correct_assembly_regions,
    output:
        good_regions=os.path.join(
            config["repeatmasker"]["output_dir"],
            "{sm}_valid_ALR_regions.500kbp.bed",
        ),
    params:
        assembly_filter="good",
    shell:
        """
        grep "{params.assembly_filter}" {input.bed} > {output.good_regions}
        """


use rule extract_and_index_fa as extract_correct_alr_regions_rm with:
    input:
        fa=rules.concat_asm.output,
        bed=rules.get_valid_regions_for_rm.output,
    output:
        seq=os.path.join(
            config["repeatmasker"]["output_dir"], "{sm}_correct_ALR_regions.500kbp.fa"
        ),
        idx=os.path.join(
            config["repeatmasker"]["output_dir"],
            "{sm}_correct_ALR_regions.500kbp.fa.fai",
        ),
    params:
        added_cmds="",
    log:
        "logs/extract_alr_regions_repeatmasker_{sm}.log",


rule rc_correct_alr_regions_rm:
    input:
        fa=rules.extract_correct_alr_regions_rm.output.seq,
    output:
        rc_seq=os.path.join(
            config["repeatmasker"]["output_dir"],
            "{sm}_correct_ALR_regions.500kbp.rc.fa",
        ),
    conda:
        "../env/tools.yaml"
    log:
        "logs/rc_alr_regions_repeatmasker_{sm}.log",
    shell:
        "seqtk seq -r {input.fa} > {output.rc_seq} 2> {log}"


rule count_complete_cens:
    input:
        expand(rules.extract_correct_alr_regions_rm.output.idx, sm=SAMPLE_NAMES),
    output:
        all_idx=os.path.join(
            config["repeatmasker"]["output_dir"],
            "all_correct_ALR_regions.500kbp.fa.fai",
        ),
        cmp_cnts=os.path.join(
            config["repeatmasker"]["output_dir"],
            "complete_cen_counts.tsv",
        ),
    params:
        # dir to remove from output strings.
        strip_dir=config["repeatmasker"]["output_dir"].replace("/", "\/") + "\/",
        num_chrs=46,
    conda:
        "../env/tools.yaml"
    log:
        "logs/count_complete_cens.log",
    shell:
        """
        rm -f {output.all_idx}
        cat {input} > {output.all_idx}

        {{ wc -l {input} {output.all_idx} | \
        sed -e 's/{params.strip_dir}//g' -e 's/_/\\t/g' | \
        awk -v OFS="\\t" '{{print $2, $1, $1/{params.num_chrs}*100}}' | \
        head -n -1;}} > {output.cmp_cnts} 2> {log}
        """


# TODO: Doc uses modified version. May need to update ID length from 50 -> 70 so ask glennis.
# https://github.com/search?q=repo%3Armhubley%2FRepeatMasker%20%2050&type=code
rule run_repeatmasker:
    input:
        seq=rules.extract_correct_alr_regions_rm.output.seq,
    output:
        os.path.join(
            config["repeatmasker"]["output_dir"],
            "{sm}",
            "{sm}_correct_ALR_regions.500kbp.fa.out",
        ),
    conda:
        "../env/tools.yaml"
    threads: config["repeatmasker"]["threads"]
    params:
        output_dir=lambda wc, output: os.path.dirname(str(output)),
        species="human",
        engine="rmblast",
    log:
        "logs/repeatmasker_{sm}.log",
    benchmark:
        "benchmarks/repeatmasker_{sm}.tsv"
    shell:
        """
        RepeatMasker \
        -engine {params.engine} \
        -species {params.species} \
        -dir {params.output_dir} \
        -pa {threads} \
        {input.seq} &> {log}
        """


rule merge_legends_for_rm:
    input:
        legends=lambda wc: expand(
            rules.new_cens_create_oriented_ctg_name_legend.output,
            sm=[wc.sm],
            ort=ORIENTATION,
        ),
        idxs=rules.extract_correct_alr_regions_rm.output.idx,
    output:
        merged_legends=temp(
            os.path.join(
                config["repeatmasker"]["output_dir"], "{sm}_merged_legend.txt"
            )
        ),
        sorted_idxs=temp(
            os.path.join(config["repeatmasker"]["output_dir"], "{sm}_sorted.fa.fai")
        ),
        trimmed_fmted_legend=os.path.join(
            config["repeatmasker"]["output_dir"], "{sm}_fmt_merged_legend.txt"
        ),
    conda:
        "../env/tools.yaml"
    log:
        "logs/merge_legends_for_rm_{sm}.log",
    shell:
        """
        # Merge legend with full contig name w/o coordinates.
        cat {input.legends} | sort -k 1 | uniq > {output.merged_legends}

        # Get contigs that are good and their coordinates.
        {{ sed 's/:/\\t/g' {input.idxs} | \
        awk -v OFS="\\t" '{{print $1, $2}}' | \
        sort -k 1 | uniq;}} > {output.sorted_idxs} 2> {log}

        # Join on original contig name and get final contig name.
        # (sm)_(chr)_(ctg):(start)-(end)
        {{ join -1 1 -2 1 {output.sorted_idxs} {output.merged_legends} | \
        awk -v OFS="\\t" '{{ print $1, $3":"$2 }}';}} > {output.trimmed_fmted_legend} 2>> {log}
        """


# TODO: rules.format_repeatmasker_output can be joined here probably.
rule rename_contig_name_repeatmasker:
    input:
        rm_out=rules.run_repeatmasker.output,
        legend=rules.merge_legends_for_rm.output.trimmed_fmted_legend,
    output:
        renamed_out=os.path.join(
            config["repeatmasker"]["output_dir"],
            "{sm}",
            "{sm}_renamed_correct_ALR_regions.500kbp.fa.out",
        ),
    run:
        import csv

        with (
            open(str(input.legend), "rt") as legend_fh,
            open(str(input.rm_out), "rt") as rm_out_fh,
            open(str(output.renamed_out), "wt") as rm_renamed_out_fh,
        ):
            contig_name_legend = dict(
                tuple(line.strip().split("\t")) for line in legend_fh.readlines()
            )
            headers = [next(rm_out_fh) for _ in range(3)]
            rm_renamed_out_writer = csv.writer(rm_renamed_out_fh, delimiter="\t")
            # Write headers to file handle without csv writer to not screw up formatting.
            rm_renamed_out_fh.writelines(headers)

            for rm_row in rm_out_fh.readlines():
                rm_row = rm_row.strip().split()
                # From: {ctg}:{start}-{end}
                # To: {sm}_{chr}_{ctg}
                query_seq_name, coords = rm_row[4].split(":")
                new_query_seq_name = contig_name_legend[query_seq_name]
                rm_row[4] = new_query_seq_name
                rm_renamed_out_writer.writerow(rm_row)


# Run repeatmasker on reference t2t-chm13 as a control.
use rule run_repeatmasker as run_repeatmasker_ref with:
    input:
        seq=config["align_asm_to_ref"]["config"]["ref"][REF_NAME],
    output:
        os.path.join(
            config["repeatmasker"]["output_dir"],
            REF_NAME,
            f"{REF_NAME}_correct_ALR_regions.500kbp.fa.out",
        ),
    log:
        f"logs/repeatmasker_{REF_NAME}.log",
    benchmark:
        f"benchmarks/repeatmasker_{REF_NAME}.tsv"


# Reformat repeatmasker
# |1259|28.4|7.4|5.3|GM18989_chr1_hap1-0000003:9717731-15372230|8|560|(5653940)|+|Charlie2b|DNA/hAT-Charlie|120|683|(2099)|1|
rule format_repeatmasker_output:
    input:
        expand(
            rules.rename_contig_name_repeatmasker.output,
            sm=SAMPLE_NAMES,
        ),
    output:
        os.path.join(
            config["repeatmasker"]["output_dir"],
            "all_samples_correct_ALR_regions.500kbp.fa.out",
        ),
    log:
        "logs/format_repeatmasker_output_all_samples.log",
    conda:
        "../env/tools.yaml"
    shell:
        """
        {{ cat {input} | \
        sed -e 's/:/\\t/g' -e 's/-/\\t/g' | \
        awk -v OFS="\\t" '{{print $1, $2, $3, $4, $5"-"$6":"$7"-"$8, $9, $10, $11, $12, $13,
$14, $15, $16, $17, $18}}' | \
        tail -n +4;}} > {output} 2> {log}
        """


rule format_add_control_repeatmasker_output:
    input:
        ref_rm_output=(
            config["repeatmasker"]["ref_repeatmasker_output"]
            if config["repeatmasker"].get("ref_repeatmasker_output")
            else rules.run_repeatmasker_ref.output
        ),
        sample_rm_output=rules.format_repeatmasker_output.output,
    output:
        os.path.join(
            config["repeatmasker"]["output_dir"],
            "all_samples_and_ref_correct_ALR_regions.500kbp.fa.out",
        ),
    log:
        "logs/format_add_control_repeatmasker_output.log",
    conda:
        "../env/tools.yaml"
    shell:
        """
        # Copy file and append reference repeatmasker output.
        cp {input.sample_rm_output} {output} 2> {log}
        {{ grep "chr" {input.ref_rm_output} | \
        sed -e 's/chr/chm13_chr/g' -e 's/:/\\t/g' -e 's/-/\\t/g' | \
        awk -v OFS="\\t" \
        '{{print $1, $2, $3, $4, $5":"$6"-"$7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17}}';}} >> {output} 2> {log}
        """


rule format_add_censat_annot_repeatmasker_output:
    input:
        cen_sat_rm_output=config["repeatmasker"]["censat_annot_hor_output"],
        sample_ctrl_rm_output=rules.format_add_control_repeatmasker_output.output,
    output:
        os.path.join(
            config["repeatmasker"]["output_dir"],
            "all_correct_ALR_regions.500kbp.fa.out",
        ),
    params:
        chr_to_add="chrY",
    log:
        "logs/format_add_censat_annot_repeatmasker_output.log",
    conda:
        "../env/tools.yaml"
    shell:
        """
        cp {input.sample_ctrl_rm_output} {output} 2> {log}

        {{ grep "{params.chr_to_add}" {input.cen_sat_rm_output} | \
        sed -e 's/chr/chm13_chr/g' -e 's/:/\\t/g' -e 's/-/\\t/g' | \
        awk -v OFS="\\t" '{{print $1, $2, $3,
$4, $5":"$6"-"$7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17}}';}} >> {output} 2> {log}
        """


rule reverse_complete_repeatmasker_output:
    input:
        rules.format_add_censat_annot_repeatmasker_output.output,
    output:
        os.path.join(
            config["repeatmasker"]["output_dir"],
            "all_correct_ALR_regions.500kbp.rc.fa.out",
        ),
    log:
        "logs/reverse_complete_repeatmasker_output.log",
    conda:
        "../env/tools.yaml"
    shell:
        """
        {{ tac {input} | \
        sed -e 's/:/\\t/g' -e 's/-/\\t/g' | \
        awk -v OFS="\\t" '{{print $1, $2, $3, $4, $5"-"$6":"$7"-"$8, $8-$7-$10+1,
$8-$7-$9+1, $11, $12, $13, $14, $15, $16, $17, $18}}' | \
        sed 's/chr/rc_chr/g' | \
        grep -v "chm13";}} > {output} 2> {log}
        """


rule extract_rm_out_by_chr:
    input:
        rm_out=rules.format_add_censat_annot_repeatmasker_output.output,
    output:
        rm_out_by_chr=os.path.join(
            config["repeatmasker"]["output_dir"], "{chr}_cens.fa.out"
        ),
    log:
        "logs/extract_{chr}_cens_from_rm.log",
    conda:
        "../env/tools.yaml"
    shell:
        """
        grep "{wildcards.chr}_" {input.rm_out} > {output.rm_out_by_chr} 2> {log}
        grep "{wildcards.chr}:" {input.rm_out} >> {output.rm_out_by_chr} 2>> {log}
        """


rule check_cens_status:
    input:
        script="workflow/scripts/check_cens_status.py",
        rm_out=rules.extract_rm_out_by_chr.output.rm_out_by_chr,
        rm_ref=(
            config["repeatmasker"]["ref_repeatmasker_output"]
            if config["repeatmasker"].get("ref_repeatmasker_output")
            else rules.run_repeatmasker_ref.output
        ),
    output:
        cens_status=os.path.join(
            config["repeatmasker"]["output_dir"], "{chr}_cens_status.tsv"
        ),
    log:
        "logs/check_cens_status_{chr}.log",
    conda:
        "../env/py.yaml"
    shell:
        """
        python {input.script} \
        -i {input.rm_out} \
        -r {input.rm_ref} \
        -o {output} 2> {log}
        """


rule create_correct_oriented_cens:
    input:
        rm_chr_out=rules.extract_rm_out_by_chr.output.rm_out_by_chr,
        rm_rev_out=rules.reverse_complete_repeatmasker_output.output,
        cens_correction_list=rules.check_cens_status.output.cens_status,
    output:
        reoriented_rm_out=os.path.join(
            config["repeatmasker"]["output_dir"], "reoriented_{chr}_cens.fa.out"
        ),
    log:
        "logs/create_correct_oriented_{chr}_cens_list.log",
    conda:
        "../env/tools.yaml"
    shell:
        """
        contigs_to_reverse=$(awk '{{ if ($3=="rev") {{ print $1 }} }}' {input.cens_correction_list} | sort | uniq)
        joined_contigs_to_reverse=$(echo "${{contigs_to_reverse[*]}}" | tr '\\n' '|' | sed 's/.$//')
        joined_contigs_to_reverse_rc=$(echo "${{joined_contigs_to_reverse[@]}}" | sed 's/chr/rc_chr/g' )

        # If nothing to reverse, just copy file.
        if [ -z $joined_contigs_to_reverse ]; then
            cp {input.rm_chr_out} {output.reoriented_rm_out}
        else
            # non-matching so everything correctly oriented
            grep -vE "$joined_contigs_to_reverse" {input.rm_chr_out} > {output.reoriented_rm_out}
            grep -E "$joined_contigs_to_reverse_rc" {input.rm_rev_out} >> {output.reoriented_rm_out}
        fi
        """


rule fix_incorrect_merged_legend:
    input:
        cens_correction_list=expand(
            rules.check_cens_status.output.cens_status, chr=CHROMOSOMES
        ),
        merged_legend=rules.merge_legends_for_rm.output.trimmed_fmted_legend,
    output:
        corrected_legend=os.path.join(
            config["repeatmasker"]["output_dir"], "{sm}_corrected_merged_legend.txt"
        ),
    run:
        import csv

        cens_renamed = {}
        incomplete_cens = set()
        for file in input.cens_correction_list:
            with open(str(file)) as cens_list_fh:
                reader_cens_renamed = csv.reader(cens_list_fh, delimiter="\t")
                for old, new, ort, is_partial in reader_cens_renamed:
                    cens_renamed[old] = new
                    if is_partial == "true":
                        incomplete_cens.add(old)

        with (
            open(str(input.merged_legend)) as merged_legend_fh,
            open(str(output.corrected_legend), "wt") as out_merged_legend_fh,
        ):
            # Write legend.
            writer_legend = csv.writer(out_merged_legend_fh, delimiter="\t")
            for old, new in csv.reader(merged_legend_fh, delimiter="\t"):
                # Skip incomplete cens.
                if new in incomplete_cens:
                    continue
                new_contig_name = cens_renamed.get(new, new)
                writer_legend.writerow((old, new_contig_name))


rule fix_incorrect_mapped_cens:
    input:
        cens_correction_list=rules.check_cens_status.output.cens_status,
        reoriented_rm_out=rules.create_correct_oriented_cens.output,
    output:
        corrected_cens_list=os.path.join(
            config["repeatmasker"]["output_dir"], "corrected_{chr}_cens.list"
        ),
        corrected_rm_out=os.path.join(
            config["repeatmasker"]["output_dir"], "corrected_{chr}_cens.fa.out"
        ),
    run:
        import csv

        with (
            open(str(input.cens_correction_list)) as cens_list_fh,
            open(str(input.reoriented_rm_out)) as rm_fh,
            open(str(output.corrected_cens_list), "wt") as out_cens_list_fh,
            open(str(output.corrected_rm_out), "wt") as out_rm_fh,
        ):
            reader_cens_renamed = csv.reader(cens_list_fh, delimiter="\t")
            cens_renamed = {}
            incomplete_cens = set()
            for old, new, ort, is_partial in reader_cens_renamed:
                cens_renamed[old] = new
                if is_partial == "true":
                    incomplete_cens.add(old)

            writer_rm_out = csv.writer(out_rm_fh, delimiter="\t")
            reader_rm_out = csv.reader(rm_fh, delimiter="\t")
            new_contigs = set()
            for line in reader_rm_out:
                contig_name = line[4]
                # Skip incomplete cens.
                if contig_name in incomplete_cens:
                    continue

                new_contig_name = cens_renamed.get(contig_name, contig_name)
                line[4] = new_contig_name
                new_contigs.add(new_contig_name)
                writer_rm_out.writerow(line)

            for c in new_contigs:
                out_cens_list_fh.write(f"{c}\n")


rule plot_cens_from_rm_by_chr:
    input:
        script="workflow/scripts/repeatStructure_onlyRM.R",
        rm_out=rules.fix_incorrect_mapped_cens.output.corrected_rm_out,
    output:
        repeat_plot_by_chr=os.path.join(
            config["repeatmasker"]["output_dir"], "hgsvc3_{chr}_cens.additional.pdf"
        ),
    log:
        "logs/plot_{chr}_cens_from_rm.log",
    conda:
        "../env/r.yaml"
    shell:
        """
        Rscript {input.script} {input.rm_out} {output.repeat_plot_by_chr} 2> {log}
        """


use rule plot_cens_from_rm_by_chr as plot_og_cens_from_rm_by_chr with:
    input:
        script="workflow/scripts/repeatStructure_onlyRM.R",
        rm_out=rules.extract_rm_out_by_chr.output.rm_out_by_chr,
    output:
        repeat_plot_by_chr=os.path.join(
            config["repeatmasker"]["output_dir"], "hgsvc3_{chr}_cens.original.pdf"
        ),
    log:
        "logs/plot_{chr}_cens_from_rm_og.log",
