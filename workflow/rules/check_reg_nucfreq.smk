
rule make_bed_files_for_plot:
    input:
        # TODO: All fasta index files concatenated. Can we just get the individual section?
        "*ALR.fa.fai",
    output:
        "{sample}_ALR_regions.500kp.bed",
    conda:
        "../env/py.yaml"
    log:
        "logs/make_bed_files_for_plot.log",
    shell:
        """
        sorted_fname="{wildcards.sample}_ALR_regions_sorted.500kbp.bed"
        collapsed_fname="{wildcards.sample}_ALR_regions_collapsed.500kbp.bed"

        cat ../*ALR.fa.fai | \
        grep {wildcards.sample} | \
        sed 's/hap/haplotype/g' | \
        sed 's/un/unassigned/g' | \
        sed 's/_/\t/g' | \
        sed 's/:/\t/g' | \
        sed 's/-/\t/g' | \
        awk -v OFS="\t" '{print $3"-"$4, $5, $6, $2}' | \
        sort | \
        uniq > $sorted_fname 2> {log}

        # Collapse by group.
        ./bedminmax.py -i $sorted_fname -o $collapsed_fname 2> {log}

        # Format
        awk -v OFS="\t" '{print $1, $2, $3, $3-$2, $4}' $collapsed_fname > {output} 2> {log}
        """


rule gen_nucfreq_plot:
    input:
        script="workflow/scripts/NucPlot.py",
        bam_file="/net/eichler/vol27/projects/AlphaSatelliteMapping/nobackups/FindingAlphaSat/hgsvc3/map_align_hifi/results/{sample}_hifi2v1.4.bam",
        alr_regions="{sample}_ALR_regions.500kbp.bed",
    output:
        plot="plots/${sample}_hifi2v1.4.cens.png ",
    conda:
        "env/py.yaml"
    params:
        ylim=100,
        height=4,
    log:
        "logs/run_nucfreq_{sample}.log",
    shell:
        """
        python {input.script} \
        -y {params.ylim} \
        {input.bam_file} \
        {output} \
        --bed {input.alr_regions} \
        --height {params.height} &> {log}
        """


# Then review plots manually.
