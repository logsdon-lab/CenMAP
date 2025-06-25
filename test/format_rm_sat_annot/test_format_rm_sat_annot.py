import os
import pytest
from test.helpers.integration import run_integration_test


WD = os.path.dirname(__file__)


@pytest.mark.parametrize(
    ["infile", "expected"],
    [
        (
            f"{WD}/input/HG00731_correct_ALR_regions.fa.reformatted.out",
            f"{WD}/expected/HG00731_correct_ALR_regions.fa.reformatted.bed",
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
