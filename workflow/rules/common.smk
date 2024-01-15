def get_samples() -> pd.DataFrame:
    df = pd.read_csv(config["align_asm_to_ref"]["config"]["tbl"], sep="\t")
    df.asm = df.asm.map(os.path.abspath)
    df["asm"] = df.asm.str.split(",")
    df = df.explode("asm")
    df["num"] = df.groupby(level=0).cumcount() + 1
    df.set_index(df["sample"] + "_" + df["num"].astype(str), inplace=True)
    return df


def get_ref_name() -> str | None:
    if ref_genome_name := config.get("ref_genome_name"):
        return ref_genome_name
    else:
        REFS = list(config["align_asm_to_ref"]["config"]["ref"].keys())
        assert len(REFS) == 1, "Only one reference genome expected."
        config["ref_genome_name"] = REFS[0]
        return REFS[0]
