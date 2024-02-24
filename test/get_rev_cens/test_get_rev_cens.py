import pytest
import subprocess


@pytest.mark.parametrize(
    "input_rm_out,expected_rc_list",
    [
        (
            "test/get_rev_cens/chr21_cens.fa.out",
            "test/get_rev_cens/correct_chr21_cens.tsv",
        ),
        (
            "test/get_rev_cens/chr22_cens.fa.out",
            "test/get_rev_cens/correct_chr22_cens.tsv",
        ),
    ],
)
def test_get_rev_cens(input_rm_out: str, expected_rc_list: str):
    process = subprocess.run(
        [
            "python",
            "workflow/scripts/get_reversed_cens.py",
            "-i",
            input_rm_out,
            "-r",
            "data/annotations/chm13_chm1_cens_v21.trimmed.fa.noheader.out",
        ],
        capture_output=True,
        check=True,
    )
    res = sorted(
        line.split("\t") for line in process.stdout.decode().split("\n") if line
    )
    with open(expected_rc_list, "rt") as exp_res_fh:
        exp_res = sorted(
            line.strip().split("\t") for line in exp_res_fh.readlines() if line
        )
        assert res == exp_res
