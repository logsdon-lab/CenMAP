import pytest
import subprocess


@pytest.mark.parametrize(
    "input,expected",
    [
        (
            "data/annotations/AS-HOR-vs-chm13_cens_v18.correctcoords.stv_row.all.bed",
            "test/calc_hor_length/chm13_hor_length.tsv",
        )
    ],
)
def test_calc_hor_length(input: str, expected: str):
    process = subprocess.run(
        ["python", "workflow/scripts/calculate_HOR_length.py", "-i", input],
        capture_output=True,
        check=True,
    )
    res = [line.split("\t") for line in process.stdout.decode().split("\n") if line]
    with open(expected, "rt") as exp_res_fh:
        exp_res = [line.strip().split("\t") for line in exp_res_fh.readlines() if line]
        assert res == exp_res
