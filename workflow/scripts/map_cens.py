import sys
import argparse
from typing import Any, Callable, Iterable
from collections import defaultdict, deque
import polars as pl

ALN_HEADER = [
    "reference_name",
    "reference_start",
    "reference_end",
    "reference_length",
    "strand",
    "query_name",
    "query_start",
    "query_end",
    "query_length",
    "perID_by_matches",
    "perID_by_events",
    "perID_by_all",
    "matches",
    "mismatches",
    "deletion_events",
    "insertion_events",
    "deletions",
    "insertions",
    "mon_reference_name",
    "mon_start",
    "mon_end",
    "arm",
]
ACRO_CHRS = {"chr21", "chr22", "chr13", "chr14", "chr15"}

Interval = tuple[int, int, Any]
                 
def merge_itvs(
    itvs: Iterable[Interval],
    dst: int = 1,
    fn_cmp: Callable[[Interval, Interval], bool] | None = None,
    fn_merge_itv: Callable[[Interval, Interval], Interval] | None = None,
) -> list[Interval]:
    if not fn_cmp:
        fn_cmp = lambda x, y: True
    if not fn_merge_itv:
        fn_merge_itv = lambda x, y: (x[0], y[1], None)

    final_itvs = []
    sorted_itvs = deque(sorted(itvs))
    while sorted_itvs:
        try:
            itv_1 = sorted_itvs.popleft()
        except IndexError:
            break
        try:
            itv_2 = sorted_itvs.popleft()
        except IndexError:
            final_itvs.append(itv_1)
            break
        dst_between = itv_2[0] - itv_1[1]
        passes_cmp = fn_cmp(itv_1, itv_2)
        if dst_between <= dst and passes_cmp:
            sorted_itvs.appendleft(fn_merge_itv(itv_1, itv_2))
        else:
            final_itvs.append(itv_1)
            sorted_itvs.appendleft(itv_2)

    return final_itvs

def main():
    ap = argparse.ArgumentParser("Map cens from alignments.")
    ap.add_argument(
        "-i",
        "--input",
        help="Input alignment bed intersected with chr p and q arms.",
        default=sys.stdin,
        type=argparse.FileType("rb"),
    )
    ap.add_argument("-o", "--output", default=sys.stdout, type=argparse.FileType("wt"))
    ap.add_argument("-m", "--merge", default=5_000_000, type=int, help="Merge pq arm alignments by n bps. ")
    ap.add_argument("-s", "--slop", default=500_000, type=int, help="Add n bps to ends. ")
    ap.add_argument("--allow_pq_mismatch", action="store_true", help="Allow non-matching chromosome pq arm pairs.")

    args = ap.parse_args()

    df = pl.read_csv(
        args.input, separator="\t", new_columns=ALN_HEADER, has_header=False
    )

    for qname, df_query in df.sort(by="query_name").group_by(["query_name"], maintain_order=True):
        qname = qname[0]
        # To pass, a given query sequence must have both p and q arms mapped.
        # The acrocentrics are exceptions and only require the q-arm.
        df_query = df_query.filter(
            pl.when(~pl.col("reference_name").is_in(ACRO_CHRS))
            .then(
                pl.all_horizontal(
                    (pl.col("arm") == "p-arm").any(), (pl.col("arm") == "q-arm").any()
                )
            )
            .otherwise(pl.all_horizontal((pl.col("arm") == "q-arm").any()))
        ).sort(by=["query_start"])
        if df_query.is_empty():
            print(f"Skipped {qname}.", file=sys.stderr)
            continue

        df_pqarms = df_query.filter(
            (pl.col("arm") == "p-arm") | (pl.col("arm") == "q-arm")
        ).select(
            "reference_name", "strand", "matches", "query_start", "query_end", "arm"
        )

        arms = defaultdict(dict)
        for prts, df_arm in df_pqarms.partition_by(["reference_name", "arm"], as_dict=True).items():
            ref, arm = prts
            itvs = merge_itvs(
                df_arm.select("query_start", "query_end", "strand").iter_rows(),
                dst=args.merge,
            )
            arms[arm][ref] = itvs
        
        for ref_chrom, itvs in arms["q-arm"].items():
            if len(itvs) > 1:
                print(f"Multiple {ref_chrom} q-arms after merging. {itvs}", file=sys.stderr)
            qarm_itv = itvs[0]

            if ref_chrom in ACRO_CHRS or args.allow_pq_mismatch:
                parm_itvs = [itv for pitvs in arms["p-arm"].values() for itv in pitvs]
            else:
                parm_itvs = arms["p-arm"][ref_chrom]
            try:
                closest_parm_itv = min(parm_itvs, key=lambda x: abs(x[0] - qarm_itv[1]))
            except ValueError:
                # No p-arm
                continue
            
            st = min(closest_parm_itv[0], qarm_itv[0])
            end = max(closest_parm_itv[1], qarm_itv[1])
            if qarm_itv[2] == "+" or closest_parm_itv[2] == "+":
                is_rc = "false"
            else:
                is_rc = "true"

            row = [
                qname, st - args.slop, end + args.slop, ref_chrom, end-st, is_rc
            ]
            print("\t".join(str(e) for e in row), file=args.output)


if __name__ == "__main__":
    main()
