import pytest
from test.helpers.integration import run_integration_test


@pytest.mark.parametrize(
    ["rm_out", "original_faidx", "renamed_faidx", "expected"],
    [
        (
            "test/rename_rm/HG01596_renamed.fa.out",
            "test/rename_rm/HG01596.fa.fai",
            "test/rename_rm/HG01596_renamed.fa.fai",
            "test/rename_rm/HG01596_expected.fa.out",
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
