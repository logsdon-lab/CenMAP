import pytest
import subprocess
import os
import pandas as pd


@pytest.mark.parametrize(
    "input_bed,expected_bed",
    [
        (
            "test/bedminmax/HG00171_cens.bed",
            "test/bedminmax/HG00171_centromeric_contigs.bed",
        )
    ],
)
def test_bedminmax(input_bed: str, expected_bed: str):
    test_fmt_tsv = "test/bedminmax/HG00171_cens_fmt_temp.tsv"
    test_tsv = "test/bedminmax/HG00171_cens_temp.tsv"

    cols = [6, 7, 8, 1, 5]
    df = pd.read_csv(input_bed, sep="\t", header=None, usecols=list(range(0, 9)))
    df = df[[c - 1 for c in cols]]
    df["length"] = df[2] - df[1]
    breakpoint()
    df.to_csv(test_fmt_tsv, sep="\t", header=None, index=False)

    subprocess.run(
        ["python3", "workflow/scripts/bedminmax.py", "-i", test_fmt_tsv, "-o", test_tsv]
    )

    df_output_bed = pd.read_csv(test_fmt_tsv, sep="\t", header=None)
    df_exp_bed = pd.read_csv(expected_bed, sep="\t", header=None)

    breakpoint()
    pd.testing.assert_frame_equal(df_output_bed, df_exp_bed)

    os.remove(test_tsv)
    os.remove(test_fmt_tsv)
