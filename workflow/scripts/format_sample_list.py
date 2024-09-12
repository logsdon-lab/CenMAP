import sys
import polars as pl


def main():
    sample_sheet = sys.argv[1]
    (
        pl.read_csv(sample_sheet, separator="\t")
        .with_columns(pl.col("hifi").str.split(","))
        .select("sample_id", "hifi")
        .explode("hifi")
        .with_columns(
            pl.col("hifi")
            .str.strip_chars("[]")
            .str.replace_all("'", "", literal=True)
            .str.replace_all(" ", "", literal=True)
        )
        .write_csv(sys.stdout, include_header=False, separator="\t")
    )


if __name__ == "__main__":
    raise SystemExit(main())
