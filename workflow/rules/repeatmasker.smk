# TODO: Only run only correctly assembled regions.
# need perl as well.
# TODO: Doc uses modified version. May need to update ID length from 50 -> 70 so ask glennis.
# https://github.com/search?q=repo%3Armhubley%2FRepeatMasker%20%2050&type=code
rule run_repeatmasker:
    input:
        alr_seqs="{sample}_correct_ALR_regions.500kbp.fa",
    output:
        directory(
            "/net/eichler/vol27/projects/AlphaSatelliteMapping/nobackups/FindingAlphaSat/hgsvc3/repeatmasker"
        ),
    conda:
        "env/env.yaml"
    threads: 40
    params:
        species="human",
    log:
        "logs/repeatmasker_{sample}.log",
    shell:
        """
        RepeatMasker \
        -species {param.species} \
        -dir {output} \
        -pa {threads} \
        {input.alr_seqs} &> {log}
        """
