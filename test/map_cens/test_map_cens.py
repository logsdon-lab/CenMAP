import pytest
from test.helpers.integration import run_integration_test


@pytest.mark.parametrize(
    ["infile", "expected"],
    [
        (
            f"test/map_cens/input/{sm}_cens.bed",
            f"test/map_cens/expected/{sm}_expected_cens.bed",
        )
        # HG00512 - no p-arm
        for sm in ["HG00512", "NA19705", "NA20355"]
    ],
)
def test_map_cens(infile: str, expected: str):
    run_integration_test(
        "python", "workflow/scripts/map_cens.py", "-i", infile, expected_output=expected
    )
