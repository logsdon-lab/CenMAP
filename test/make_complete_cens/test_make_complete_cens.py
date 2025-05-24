import pytest
from test.helpers.integration import run_integration_test


@pytest.mark.parametrize(
    ["status", "fai", "expected", "use_new_name"],
    [
        (
            f"test/make_complete_cens/input/{sm}.tsv",
            f"test/make_complete_cens/input/{sm}.fa.fai",
            f"test/make_complete_cens/expected/{sm}{"_new_name" if use_new_name else ""}.bed",
            use_new_name,
        )
        for sm in ["HG02583"]
        # See HG02583_chr21_HG02583#1#JBHIIU010000054.1
        for use_new_name in [True, False]
    ],
)
def test_make_complete_cens(status: str, fai: str, expected: str, use_new_name: bool):
    args = [
        "python", "workflow/scripts/make_complete_cens_bed.py", "-i", status, "-f", fai
    ]
    if use_new_name:
        args.append("--use_new_name")

    run_integration_test(*args, expected_output=expected)
