import pytest
from test.helpers.integration import run_integration_test


@pytest.mark.parametrize(
    ["infile", "expected"],
    [
        (
            "test/format_rm_sat_annot/input/HG00731_correct_ALR_regions.fa.reformatted.out",
            "test/format_rm_sat_annot/expected/HG00731_correct_ALR_regions.fa.reformatted.bed",
        )
    ],
)
def test_format_rm_sat_annot(infile: str, expected: str):
    run_integration_test(
        "python",
        "workflow/scripts/format_rm_sat_annot.py",
        "-i",
        infile,
        expected_output=expected,
    )
