import pytest
from test.helpers.integration import run_integration_test


@pytest.mark.parametrize(
    ["infile", "thresholds", "chromosome", "expected"],
    [
        (
            "test/filter_dnabrnn_output/input/NA19129_chr16_hap2_forward.bed",
            "test/filter_dnabrnn_output/input/dnabrnn_thresholds.json",
            "chr16",
            "test/filter_dnabrnn_output/expected/NA19129_chr16_hap2_forward_filtered.bed",
        ),
        (
            "test/filter_dnabrnn_output/input/NA19317_chr16_hap2_reverse.bed",
            "test/filter_dnabrnn_output/input/dnabrnn_thresholds.json",
            "chr16",
            "test/filter_dnabrnn_output/expected/NA19317_chr16_hap2_reverse_filtered.bed",
        ),
        (
            "test/filter_dnabrnn_output/input/NA19317_chr3_hap2_forward.bed",
            "test/filter_dnabrnn_output/input/dnabrnn_thresholds.json",
            "chr3",
            "test/filter_dnabrnn_output/expected/NA19317_chr3_hap2_forward_filtered.bed",
        ),
        (
            "test/filter_dnabrnn_output/input/HG03248_chr16_hap1_reverse.bed",
            "test/filter_dnabrnn_output/input/dnabrnn_thresholds.json",
            "chr16",
            "test/filter_dnabrnn_output/expected/HG03248_chr16_hap1_reverse_filtered.bed",
        ),
    ],
)
def test_filter_dnabrnn_output(
    infile: str, thresholds: str, chromosome: str, expected: str
):
    run_integration_test(
        "python",
        "workflow/scripts/filter_dnabrnn_output.py",
        "-i",
        infile,
        "-t",
        thresholds,
        "-c",
        chromosome,
        expected_output=expected,
    )
