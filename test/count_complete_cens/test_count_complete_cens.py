import os
import pytest
from test.helpers.integration import run_integration_test


WD = os.path.dirname(__file__)


@pytest.mark.parametrize(
    ["infile", "expected"],
    [
        (
            f"{WD}/input/all_AS-HOR_lengths.tsv",
            f"{WD}/expected/expected_counts.tsv",
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
