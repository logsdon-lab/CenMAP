
BEGIN {
    OFS="\t"
    ALR_COLOR="#522758"
} {
    name=$5; start=$6; end=$7; rType=$10; rClass=$11;

    # Find contig coordinates
    match(name, "^(.+):", ctg_name);

    # Split repeat class and replace specific repeat types.
    split(rClass, split_rClass, "/" );
    new_rClass=split_rClass[1];
    if (rClass == "Satellite/centr" || rClass == "Satellite") {
        new_rClass=rType
    }
    switch (new_rClass) {
        case "SAR":
            new_rClass="HSat1A";
            break;
        case "HSAT":
            new_rClass="HSat1B";
            break;
        case "HSATII":
            new_rClass="HSat2";
            break;
        case "(CATTC)n":
            new_rClass="HSat3";
            break;
        case "(GAATG)n":
            new_rClass="HSat3";
            break;
        default:
            break;
    }

    # Set action for NucFlag
    color=""
    if (new_rClass == "ALR/Alpha") {
        color=ALR_COLOR
    }
    print ctg_name[1], start, end, new_rClass, 0, ".", start, end, color
}
