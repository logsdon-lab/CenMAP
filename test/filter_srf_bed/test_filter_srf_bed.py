import os
import pytest
from test.helpers.integration import run_integration_test


WD = os.path.dirname(__file__)


@pytest.mark.parametrize(
    ["input_srf", "input_trf", "expected"],
    [
        (
            f"{WD}/input/srf_{sm}.bed",
            f"{WD}/input/monomers_{sm}.tsv",
            f"{WD}/expected/{sm}.bed",
        )
        for sm in ["HG008-N", "HG008-T"]
    ],
)
def test_filter_srf_bed(input_srf: str, input_trf: str, expected: str):
    args = [
        "python",
        "workflow/scripts/filter_srf_bed.py",
        "-i",
        input_srf,
        "-m",
        input_trf,
    ]

    run_integration_test(*args, expected_output=expected)
