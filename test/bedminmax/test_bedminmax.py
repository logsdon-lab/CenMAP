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

    # awk -v OFS="\t" { print $6, $7, $8, $4, $1, $5 }
    cols = [c - 1 for c in (6, 7, 8, 4, 1, 5)]
    df = pd.read_csv(input_bed, sep="\t", header=None, usecols=cols)
    df = df.reindex(columns=cols)
    df.to_csv(test_fmt_tsv, sep="\t", header=None, index=False)

    subprocess.run(
        ["python3", "workflow/scripts/bedminmax.py", "-i", test_fmt_tsv, "-o", test_tsv]
    )

    # awk -v OFS="\\t" '{ print $1, $2, $3, $5, $6, $3-$2 }' | \
    # sort -k1,1
    df_output_bed = pd.read_csv(test_tsv, sep="\t", header=None)
    df_output_bed[3] = df_output_bed[2] - df_output_bed[1]
    df_output_bed = df_output_bed.reindex(columns=[0, 1, 2, 4, 5, 3])
    df_output_bed.columns = pd.Index(range(6))
    df_output_bed.sort_values(by=[1], inplace=True)
    df_output_bed.reset_index(drop=True, inplace=True)

    df_exp_bed = pd.read_csv(expected_bed, sep="\t", header=None)
    df_exp_bed.columns = pd.Index(range(6))
    df_exp_bed.sort_values(by=[1], inplace=True)
    df_exp_bed.reset_index(drop=True, inplace=True)

    pd.testing.assert_frame_equal(df_exp_bed, df_output_bed)

    os.remove(test_tsv)
    os.remove(test_fmt_tsv)
