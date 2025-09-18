import os
import pytest
from test.helpers.integration import run_integration_test


WD = os.path.dirname(__file__)


@pytest.mark.parametrize(
    ["input_entropy", "input_rm", "expected"],
    [
        (
            f"{WD}/input/{sm}_entropy.bed",
            f"{WD}/input/{sm}_rm.out",
            f"{WD}/expected/{sm}.bed",
        )
        # HG00358 - chr4 (ALR + HSAT1A)
        # HG02059 - chr20 ALR fragment (Zero entropy, no dip)
        # HG02818 - chr2 complete centromere
        # HG00513 - chr20 partial centromere
        for sm in ["HG00358", "HG02059", "HG02818", "HG00513", "empty"]
    ],
)
def test_filter_srf_bed(input_entropy: str, input_rm: str, expected: str):
    args = [
        "python",
        "workflow/scripts/filter_entropy_bed.py",
        "-i",
        input_entropy,
        "-r",
        input_rm,
    ]

    run_integration_test(*args, expected_output=expected)
