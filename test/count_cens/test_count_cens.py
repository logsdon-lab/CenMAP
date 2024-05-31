import pytest
import subprocess


@pytest.mark.parametrize(
    ["infile", "expected"],
    [
        (
            "test/count_cens/all_AS-HOR_lengths.tsv",
            "test/count_cens/expected_counts.tsv",
        )
    ],
)
def test_count_cens(infile: str, expected: str):
    process = subprocess.run(
        [
            "python",
            "workflow/scripts/count_complete_cens.py",
            "-i",
            infile,
        ],
        capture_output=True,
        check=True,
    )
    res = sorted(
        line.strip().split("\t") for line in process.stdout.decode().split("\n") if line
    )
    with open(expected, "rt") as exp_res_fh:
        exp_res = sorted(
            line.strip().split("\t") for line in exp_res_fh.readlines() if line
        )
        assert res == exp_res
