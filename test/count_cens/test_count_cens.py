import pytest
from test.helpers.integration import run_integration_test


@pytest.mark.parametrize(
    ["infile", "expected"],
    [
        (
            "test/count_cens/all_AS-HOR_lengths.tsv",
            "test/count_cens/expected_counts.tsv",
        )
    ],
)
def test_count_cens(infile: str, expected: str):
    run_integration_test(
        "python",
        "workflow/scripts/count_complete_cens.py",
        "-i",
        infile,
        expected_output=expected,
    )
