import os
import pytest
from test.helpers.integration import run_integration_test


WD = os.path.dirname(__file__)


@pytest.mark.parametrize(
    ["infile", "expected", "add_ct"],
    [
        (
            f"{WD}/input/HG00731_correct_ALR_regions.fa.reformatted.out",
            f"{WD}/expected/HG00731_correct_ALR_regions.fa.reformatted{'.ct' if add_ct else ''}.bed",
            add_ct,
        )
        for add_ct in [True, False]
    ],
)
def test_format_rm_sat_annot(infile: str, expected: str, add_ct: bool):
    args = [
        "python",
        "workflow/scripts/format_rm_sat_annot.py",
        "-i",
        infile,
    ]
    if add_ct:
        args.append("--add_ct")

    run_integration_test(
        *args,
        expected_output=expected,
    )
