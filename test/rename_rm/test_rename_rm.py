import pytest
import subprocess


@pytest.mark.parametrize(
    ["rm_out", "original_faidx", "renamed_faidx", "expected"],
    [
        (
            "test/rename_rm/HG01596_renamed.fa.out",
            "test/rename_rm/HG01596.fa.fai",
            "test/rename_rm/HG01596_renamed.fa.fai",
            "test/rename_rm/HG01596_expected.fa.out",
        )
    ],
)
def test_rename_rm(rm_out: str, original_faidx: str, renamed_faidx: str, expected: str):
    process = subprocess.run(
        [
            "python",
            "workflow/scripts/reformat_rm.py",
            "-i",
            rm_out,
            "-of",
            original_faidx,
            "-rf",
            renamed_faidx,
        ],
        capture_output=True,
        check=True,
    )
    res = [
        line.strip().split("\t") for line in process.stdout.decode().split("\n") if line
    ]
    with open(expected, "rt") as exp_res_fh:
        exp_res = [line.strip().split("\t") for line in exp_res_fh.readlines() if line]
        assert res == exp_res
