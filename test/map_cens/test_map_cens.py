import pytest
from test.helpers.integration import run_integration_test


@pytest.mark.parametrize(
    ["infile", "expected", "allow_pq_mismatch"],
    [
        (
            f"test/map_cens/input/{sm}_cens.bed",
            f"test/map_cens/expected/{sm}_expected_cens.bed",
            allow_pq_mismatch,
        )
        # HG00512 - no p-arm
        # HG008 (non-acro translocation), NA19036 (acro misassembly) - two cens on one contig
        for sm, allow_pq_mismatch in zip(["HG00512", "NA19705", "NA20355", "HG008", "NA19036"], [False, False, False, True, True])
    ],
)
def test_map_cens(infile: str, expected: str, allow_pq_mismatch: bool):
    args = ["python", "workflow/scripts/map_cens.py", "-i", infile]
    if allow_pq_mismatch:
        args.append("--allow_pq_mismatch")

    run_integration_test(*args, expected_output=expected)
