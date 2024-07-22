import subprocess


def run_integration_test(*args: str, expected_output: str) -> None:
    process = subprocess.run(
        [*args],
        capture_output=True,
        check=True,
    )
    res = sorted(
        line.strip().split("\t") for line in process.stdout.decode().split("\n") if line
    )
    with open(expected_output, "rt") as exp_res_fh:
        exp_res = sorted(
            line.strip().split("\t") for line in exp_res_fh.readlines() if line
        )
        assert res == exp_res
