def get_ref_name() -> str:
    if ref_genome_name := config.get("ref_genome_name"):
        return ref_genome_name
    else:
        REFS = list(config["align_asm_to_ref"]["config"]["ref"].keys())
        assert len(REFS) == 1, "Only one reference genome expected."
        config["ref_genome_name"] = REFS[0]
        return REFS[0]


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
