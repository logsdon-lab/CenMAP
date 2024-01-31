import json


def get_correct_alr_regions(wc):
    """
    Read correct_alr_regions json file.
    {
        "HG00171": "HG00171_correct_ALR_regions.500kbp.fa",
        "HG00172": "HG00172_correct_ALR_regions.500kbp.fa"
    }
    """
    if correct_samples := config.get("correct_samples_map"):
        samples = correct_samples[wc.sm]
    else:
        with open(config["repeatmasker"]["correct_samples"]) as fh:
            all_samples = json.load(fh)
            config["correct_samples_map"] = all_samples
            samples = all_samples[wc.sm]

    return samples


# TODO: Prior to this step we checked regions with NucFreq to determine if there was a misassembly.
# * If it's good, do we take {sample}_ALR_regions.500kbp.bed and take subseq of ref?
# * {sample}_correct_ALR_regions.500kbp.fa
# * /net/eichler/vol27/projects/AlphaSatelliteMapping/nobackups/FindingAlphaSat/hgsvc3/repeatmasker
# TODO: Need perl as well.
# TODO: Doc uses modified version. May need to update ID length from 50 -> 70 so ask glennis.
# https://github.com/search?q=repo%3Armhubley%2FRepeatMasker%20%2050&type=code
rule run_repeatmasker:
    input:
        seq=get_correct_alr_regions,
    output:
        directory(os.path.join(config["repeatmasker"]["output_dir"], "{sm}")),
    conda:
        "../env/tools.yaml"
    threads: config["repeatmasker"]["num_threads"]
    params:
        species="human",
    log:
        "logs/repeatmasker_{sm}.log",
    benchmark:
        "benchmarks/repeatmasker_{sm}.tsv"
    shell:
        """
        RepeatMasker \
        -species {param.species} \
        -dir {output} \
        -pa {threads} \
        {input.seq} &> {log}
        """


# TODO: Grp by chr


# Run repeatmasker on reference t2t-chm13 as a control.
use rule run_repeatmasker as run_repeatmasker_ref with:
    input:
        ref=config["align_asm_to_ref"]["config"]["ref"][REF_NAME],
    output:
        directory(os.path.join(config["repeatmasker"]["output_dir"], REF_NAME)),
    log:
        f"logs/repeatmasker_{REF_NAME}.log",
    benchmark:
        f"benchmarks/repeatmasker_{REF_NAME}.tsv"
