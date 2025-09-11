import os
import pytest
from test.helpers.integration import run_integration_test


WD = os.path.dirname(__file__)


@pytest.mark.parametrize(
    ["infile", "expected"],
    [
        (
            f"{WD}/input/{sm}.bed",
            f"{WD}/expected/{sm}.tsv",
        )
        for sm in ["HG008-N", "HG008-T"]
    ],
)
def test_map_cens(infile: str, expected: str):
    args = ["python", "workflow/scripts/map_chroms.py", "-i", infile]

    run_integration_test(*args, expected_output=expected)
