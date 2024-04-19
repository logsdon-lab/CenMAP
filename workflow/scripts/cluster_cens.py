import os
import json
import argparse
import editdistance
import polars as pl
import numpy as np
import matplotlib.pyplot as plt

from collections import defaultdict
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
from scipy.cluster.hierarchy import linkage, dendrogram

COLUMNS = ["chr", "start", "stop", "mon", "_", "ort", "start2", "stop2", "color", "ctg"]


# https://gist.github.com/vals/150ec97a5b7db9c82ee9
def get_cluster_classes(den, label="ivl"):
    cluster_idxs = defaultdict(list)
    for c, pi in zip(den["color_list"], den["icoord"]):
        for leg in pi[1:3]:
            i = (leg - 5.0) / 10.0
            if abs(i - int(i)) < 1e-5:
                cluster_idxs[c].append(int(i))

    cluster_classes = {}
    for c, lst in cluster_idxs.items():
        i_lst = [den[label][i] for i in lst]
        cluster_classes[c] = i_lst

    return cluster_classes


def compute_dst(coords: tuple[int, int], words: list[str]) -> int:
    x, y = coords
    return editdistance.eval(words[x], words[y])


def main():
    ap = argparse.ArgumentParser(
        description="Cluster centromeres by their HOR monomers."
    )
    ap.add_argument(
        "-i",
        "--input",
        help=f"Input live HOR bed file. Expects columns {COLUMNS}",
        required=True,
    )
    ap.add_argument(
        "-p",
        "--output_plot",
        help="Output dendrogram plot.",
        default="cens_dendrogram.png",
    )
    ap.add_argument(
        "-d",
        "--cens_dir",
        help="Input directory with HOR monomer annotated centromere images.",
        default=None,
    )
    ap.add_argument(
        "--output_clusters",
        help="Output clusters list as JSON.",
        default="cens_cluster.json",
    )
    ap.add_argument(
        "--linkage_method",
        help="Linkage method to perform clustering.",
        default="average",
    )

    args = ap.parse_args()

    df = pl.read_csv(args.input, separator="\t", new_columns=COLUMNS)

    df_mons = df.group_by("ctg").agg(
        merged_mons=pl.when(pl.col("ctg").str.contains("rc-"))
        .then(pl.col("mon").reverse())
        .otherwise(pl.col("mon"))
    )
    upper_tri_coords = np.triu_indices(len(df_mons.get_column("merged_mons")), 1)

    computed_dsts = np.apply_along_axis(
        lambda c: compute_dst(c, df_mons.get_column("merged_mons").to_list()),
        0,
        upper_tri_coords,
    )
    clusters = linkage(computed_dsts, method=args.linkage_method, metric=None)

    den_info = dendrogram(
        clusters, labels=df_mons.get_column("ctg").to_list(), orientation="left"
    )
    plt_den = plt.gcf()

    plt_den.set_size_inches(
        30,
        # Scale height based on elements default to 8 in.
        max(8, 1.5 * len(df["ctg"].unique())),
    )
    ax_den = plt.gca()

    if args.cens_dir:
        for lbl in ax_den.get_yticklabels():
            ctg_name = lbl.get_text()
            ctg_img_path = os.path.join(args.cens_dir, f"{ctg_name}.png")
            _, y_pos = lbl.get_position()
            try:
                img = plt.imread(ctg_img_path)
                im = OffsetImage(img, zoom=0.1)
                im.image.axes = ax_den
                # TODO: This might not be correctly oriented.
                ab = AnnotationBbox(
                    im,
                    xy=(-800, y_pos - 3),
                    xycoords="data",
                    pad=0,
                    annotation_clip=False,
                )
                ax_den.add_artist(ab)
            except FileNotFoundError:
                continue

    plt.tight_layout()
    plt.savefig(args.plot, dpi=300)

    clusters = get_cluster_classes(den_info)
    with open(args.clusters, "wt") as fh:
        json.dump(clusters, fh)


if __name__ == "__main__":
    raise SystemExit(main())
