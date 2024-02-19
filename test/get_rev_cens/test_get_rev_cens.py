import pytest
import subprocess


@pytest.mark.parametrize(
    "input_rm_out,expected_rc_list",
    [
        (
            "test/get_rev_cens/chr21_cens.fa.out",
            "test/get_rev_cens/rc_chr21_cens.list",
        ),
        (
            "test/get_rev_cens/chr22_cens.fa.out",
            "test/get_rev_cens/rc_chr22_cens.list",
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
        ],
        capture_output=True,
        check=True,
    )
    res = process.stdout.decode().split()
    with open(expected_rc_list, "rt") as exp_res_fh:
        exp_res = [line.strip() for line in exp_res_fh.readlines()]

    assert res == exp_res
