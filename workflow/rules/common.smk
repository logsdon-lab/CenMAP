def get_ref_name() -> str:
    if ref_genome_name := config.get("ref_genome_name"):
        return ref_genome_name
    else:
        REFS = list(config["align_asm_to_ref"]["config"]["ref"].keys())
        assert len(REFS) == 1, "Only one reference genome expected."
        config["ref_genome_name"] = REFS[0]
        return REFS[0]
