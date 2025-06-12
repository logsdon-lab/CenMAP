import pytest
from test.helpers.integration import run_integration_test


@pytest.mark.parametrize(
    ["infile", "expected"],
    [
        (
            f"test/map_cens/input/{sm}_cens.bed",
            f"test/map_cens/expected/{sm}.bed",
        )
        for sm in ["HG008-N", "HG008-T"]
    ],
)
def test_map_cens(infile: str, expected: str):
    args = ["python", "workflow/scripts/map_cens.py", "-i", infile]

    run_integration_test(*args, expected_output=expected)
