import os
import pytest
from test.helpers.integration import run_integration_test


WD = os.path.dirname(__file__)


@pytest.mark.parametrize(
    ["rm_out", "original_faidx", "renamed_faidx", "expected"],
    [
        (
            f"{WD}/input/HG01596_renamed.fa.out",
            f"{WD}/input/HG01596.fa.fai",
            f"{WD}/input/HG01596_renamed.fa.fai",
            f"{WD}/expected/HG01596_expected.fa.out",
        )
    ],
)
def test_rename_rm(rm_out: str, original_faidx: str, renamed_faidx: str, expected: str):
    run_integration_test(
        "python",
        "workflow/scripts/reformat_rm.py",
        "-i",
        rm_out,
        "-of",
        original_faidx,
        "-rf",
        renamed_faidx,
        expected_output=expected,
    )
