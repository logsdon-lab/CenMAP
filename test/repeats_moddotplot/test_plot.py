# import subprocess
# import filecmp
# import pytest

# @pytest.mark.parametrize(
#     ["expected_png_path"],
#     [
#         ("/path/to/script1.R", "/path/to/expected1.png"),
#         ("/path/to/script2.R", "/path/to/expected2.png"),
#         ("/path/to/script3.R", "/path/to/expected3.png")
#     ]
# )
# def test_rscript_output(script_path, expected_png_path):
#     # Run the Rscript using subprocess
#     subprocess.run(['Rscript', script_path])

#     # Define the path to the actual PNG file
#     actual_png_path = script_path.replace(".R", ".png")

#     # Compare the expected and actual PNG files
#     assert filecmp.cmp(expected_png_path, actual_png_path)
