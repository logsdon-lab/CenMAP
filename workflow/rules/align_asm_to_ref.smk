# Generate table
ASM_TABLE = f"/tmp/table_{hash(workflow.basedir)}.asm.tbl"
with open(ASM_TABLE, "wt") as tbl_fh:
    tbl_fh.write("sample\tasm\n")
    for sm in SAMPLE_NAMES:
        concat_asm_fa = os.path.join(
            config["concat_asm"]["output_dir"], f"{sm}-asm-comb-dedup.fa"
        )
        tbl_fh.write(f"{sm}\t{concat_asm_fa}\n")


ALN_CFG = {
    "ref": {f"{REF_NAME}_cens": rules.extract_ref_hor_arrays.output.seq},
    "tbl": ASM_TABLE,
    "aln_threads": config["align_asm_to_ref"]["threads"],
    "mm2_opts": "-x asm20 --secondary=no -s 25000 -K 15G",
    "second_aln": "no",
    ## THE FOLLOWING OPTIONS ARE ONLY USED FOR GENE CONVERSION ANALYSIS AND NOT ALIGNMENT AND CAN BE IGNORRED
    # subset the gene conversion analysis to just these regions on the reference
    # bed: /net/eichler/vol26/projects/chm13_t2t/nobackups/Assembly_analysis/SEDEF/chm13_v1.1_plus38Y.SDs.bed
    "break_paf": 10000,
    "bed": "config/SDs.and.lowid.bed",
    "window": 1000,
    "slide": 100,
    "min_aln_len": 1000000,
}


# Align assemblies to reference.
# Pull alignment workflow from github.
# See output files here. Most of everything gets dumped in subdirs from results/{ref}/{ftype}/{asm}.{ftype}:
# https://github.com/mrvollger/asm-to-reference-alignment/blob/main/workflow/rules/reference_alignment.smk
module align_asm_to_ref:
    snakefile:
        github(
            "koisland/asm-to-reference-alignment",
            path="workflow/Snakefile",
            branch="remove_sm_num_index",
        )
    config:
        ALN_CFG


use rule * from align_asm_to_ref exclude dipcall, ideogram as asm_ref_*


# Override rules from ^ because input getter functions run before concat_asm so cannot detect file.
use rule alignment from align_asm_to_ref as asm_ref_alignment with:
    input:
        ref=rules.extract_ref_hor_arrays.output.seq,
        # Take concat asm with all types included.
        query=rules.concat_asm.output.fa,
    resources:
        mem=config["align_asm_to_ref"].get("mem", 4),


use rule alignment2 from align_asm_to_ref as asm_ref_alignment2 with:
    input:
        ref_fasta=rules.extract_ref_hor_arrays.output.seq,
        query=rules.concat_asm.output.fa,
        # Weird. Blank if passing reference rule output from above. Use string instead.
        aln="temp/{ref}/{sm}.bam",


use rule all from align_asm_to_ref as align_asm_ref_all with:
    input:
        ra=rules.asm_ref_reference_alignment.input,
        SafFire=expand(
            rules.asm_ref_SafFire.output,
            sm=SAMPLE_NAMES,
            ref=[f"{REF_NAME}_cens"],
        ),
