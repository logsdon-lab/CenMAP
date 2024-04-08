import pytest
import subprocess
import os
import pandas as pd


def standardize_df(df: pd.DataFrame) -> pd.DataFrame:
    df.columns = pd.Index(range(len(df.columns)))
    df.sort_values(by=[0, 1], inplace=True)
    df.reset_index(drop=True, inplace=True)
    return df


@pytest.mark.parametrize(
    "input_bed,expected_bed",
    [
        (
            "test/bedminmax/HG00171_cens_ordered.bed",
            "test/bedminmax/HG00171_centromeric_contigs.bed",
        )
    ],
)
def test_bedminmax(input_bed: str, expected_bed: str):
    test_tsv = "test/bedminmax/HG00171_cens_temp.tsv"
    test_new_tsv = "test/bedminmax/HG00171_cens_temp_new.tsv"
    IO_COLS = ["chr", "start", "end", "length", "name", "orientation"]
    GRP_COLS = ["chr", "length", "name", "orientation"]
    subprocess.run(
        [
            "python",
            "test/bedminmax/bedminmax_cens.py",
            "-i",
            input_bed,
            "-o",
            test_tsv,
        ]
    )
    subprocess.run(
        [
            "python",
            "workflow/scripts/filter_cen_ctgs.py",
            "bedminmax",
            "-i",
            input_bed,
            "-o",
            test_new_tsv,
            "-ci",
            *IO_COLS,
            "-co",
            *IO_COLS,
            "-g",
            *GRP_COLS,
        ]
    )
    # awk -v OFS="\\t" '{ print $1, $2, $3, $5, $6, $3-$2 }' | \
    # sort -k1,1
    df_output_bed = pd.read_csv(test_tsv, sep="\t", header=None)
    df_output_bed[3] = df_output_bed[2] - df_output_bed[1]
    df_output_bed = df_output_bed.reindex(columns=[0, 1, 2, 4, 5, 3])
    df_output_bed = standardize_df(df_output_bed)

    df_output_bed_new = pd.read_csv(test_new_tsv, sep="\t", header=None)
    df_output_bed_new[3] = df_output_bed_new[2] - df_output_bed_new[1]
    df_output_bed_new = df_output_bed_new.reindex(columns=[0, 1, 2, 4, 5, 3])
    df_output_bed_new = standardize_df(df_output_bed_new)

    df_exp_bed = pd.read_csv(expected_bed, sep="\t", header=None)
    df_exp_bed = standardize_df(df_exp_bed)

    pd.testing.assert_frame_equal(df_exp_bed, df_output_bed)
    # New script matches old script.
    pd.testing.assert_frame_equal(df_output_bed_new, df_output_bed)

    os.remove(test_tsv)
    os.remove(test_new_tsv)
