def get_hifi_read_wildcards() -> dict[str, list[str]]:
    """
    Get hifi reads by sample automatically from hifi_reads_dir.
    Expects {hifi_reads_dir}/{sample}/*.bam
    """
    reads_dir = config["nuc_freq"]["hifi_reads_dir"]
    path_pattern = os.path.join(reads_dir, "{sm}", "{flowcell_id}.bam")
    reads_run_mdata_id = glob_wildcards(path_pattern)

    samples = defaultdict(list)
    for sm, flowcell_id in zip(reads_run_mdata_id.sm, reads_run_mdata_id.flowcell_id):
        samples[sm].append(flowcell_id)
    return samples


def load_samples_df() -> pd.DataFrame:
    # Load samples.
    # https://github.com/mrvollger/asm-to-reference-alignment/blob/fe0d4d5f1ccfd28bc1e3d926df83d64f18047e92/workflow/rules/reference_alignment.smk#L7-L12
    df = pd.read_csv(config["align_asm_to_ref"]["config"]["tbl"], sep="\t")
    df.asm = df.asm.map(os.path.abspath)
    df["asm"] = df.asm.str.split(",")
    df = df.explode("asm")
    df["num"] = df.groupby(level=0).cumcount() + 1
    if len(df["num"].unique()) > 1:
        new_index = df["sample"] + "_" + df["num"].astype(str)
    else:
        new_index = df["sample"].astype(str)

    df.set_index(new_index, inplace=True)
    return df


def build_awk_cen_region_length_thr(chr_name: str) -> str:
    region_thresholds: list[list[int | None, int | None]] = (
        config["dna_brnn"]
        .get("repeat_len_thr", {})
        .get(
            chr_name, [config["dna_brnn"].get("default_repeat_len_thr", [1_000, None])]
        )
    )
    stmts = []
    for thr_min, thr_max in region_thresholds:
        if thr_min and thr_max:
            stmts.append(f"$5>{thr_min} && $5<{thr_max}")
        elif thr_min:
            stmts.append(f"$5>{thr_min}")
        else:
            stmts.append(f"$5<{thr_max}")

    return f"({' || '.join(stmts)})"
