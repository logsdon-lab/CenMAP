#!/usr/bin/env Rscript
library(ggplot2)
library(dplyr)
require(scales)
library(RColorBrewer)
library(data.table)
library(cowplot)
library(stringr)
library(glue)
library(tidyr)
library(ggnewscale)
library(argparser, quietly = TRUE)


args <- commandArgs(trailingOnly = FALSE)
wd <- normalizePath(
  dirname(sub("^--file=", "", args[grep("^--file=", args)]))
)
if (length(wd) == 0) {
  wd <- "."
}

source(glue("{wd}/r_utils/stained_glass.R"))
source(glue("{wd}/r_utils/repeats.R"))

parser <- arg_parser(
  "Plot combined HOR + satellites + moddotplot for a single centromere."
)
parser <- add_argument(
  parser,
  "--bed",
  help = "bedfile with alignment information",
  type = "character"
)
parser <- add_argument(
  parser,
  "--hor",
  help = "bedfile with HOR information",
  type = "character"
)
parser <- add_argument(
  parser,
  "--sat",
  help = "bedfile with satellite repeat information",
  type = "character"
)
parser <- add_argument(
  parser,
  "--cdr",
  help = "bedfile with CDRs",
  default = NA
)
parser <- add_argument(
  parser,
  "--methyl",
  help = "bedfile with binned methylation percents.",
  default = NA
)
parser <- add_argument(
  parser,
  "--hor_ort",
  help = "bedfile with hor ort. Can be derived from --hor but this stopgap is here because there are no good interval tree libraries in R.",
)
parser <- add_argument(
  parser,
  "--outdir",
  default = "results",
  help = "Output directory.",
  type = "character"
)
parser <- add_argument(
  parser,
  "--mer_order",
  type = "character",
  help = "HOR monomer order", default = "small"
)
parser <- add_argument(
  parser,
  "--height",
  type = "numeric",
  help = "Plot height."
)
parser <- add_argument(
  parser,
  "--width",
  type = "numeric",
  help = "Plot width."
)
args <- parse_args(parser)
bed_seq_ident <- args$bed
bed_sat_annot <- args$sat
bed_hor_mon <- args$hor
bed_hor_mon_ort <- args$hor_ort
outdir <- args$outdir
mer_order <- args$mer_order
dir.create(outdir)

df_humas_hmmer_stv_out <- read_one_humas_hmmer_input(bed_hor_mon)
# First sort by val. This sorts the dataframe but NOT the factor levels.
# larger HORs on top
# smaller HORs on top
df_humas_hmmer_stv_out <- switch(args$mer_order,
  "large" = df_humas_hmmer_stv_out %>% arrange(mer),
  "small" = df_humas_hmmer_stv_out %>% arrange(-mer),
  stop(paste("Invalid mer reordering option:", args$mer_order))
)
df_hor_ort <- read_one_hor_mon_ort_input(bed_hor_mon_ort)
# Read CDR dataframe if provided.
if (!is.na(args$cdr)) {
  df_cdr <- read_one_cdr_input(args$cdr)
  height_adj_cdr <- 0.1
} else {
  df_cdr <- NA
  height_adj_methyl_binned <- 0
}

if (!is.na(args$methyl)) {
  df_methyl_binned <- read_one_methyl_bed_input(args$methyl)
  height_adj_methyl_binned <- 0.5
} else {
  df_methyl_binned <- NA
  height_adj_methyl_binned <- 0
}

df_rm_sat_out <- read_one_repeatmasker_sat_input(bed_sat_annot)
df_seq_ident <- read_bedpe(bed_seq_ident)
# Get original reference name.
original_rname <- df_seq_ident$or[[1]]
rname <- df_seq_ident$r[[1]]
plots <- make_cen_plot(
  # Trimmed reference name.
  rname,
  df_seq_ident,
  df_humas_hmmer_stv_out,
  df_rm_sat_out,
  df_cdr,
  df_hor_ort,
  df_methyl_binned
)
plot_methyl <- plots[["methyl"]]
plot_cen <- plots[["cen"]]
plot_cen_legend <- get_legend(
  plot_cen +
  theme(
    legend.justification="center",
    legend.box.just = "left",
  )
)
plot_cen <- plot_cen + theme(legend.position = "none")
plot_hist <- plots[["hist"]]

bottom_row <- cowplot::plot_grid(
  plotlist = list(NULL, plot_cen_legend, plot_hist, NULL),
  nrow = 1,
  rel_widths = c(0.1, 1, 1, 0.1),
  labels = NULL
)

title <- make_cen_plot_title(original_rname)

if (!typeof(plot_methyl) == "logical") {
  plot <- cowplot::plot_grid(
    title, plot_methyl, plot_cen, bottom_row,
    nrow=4,
    ncol=1,
    rel_heights = c(0.05, 0.3, 1.2, 0.8),
    labels = NULL
  )
} else {
  plot <- cowplot::plot_grid(
    title, plot_cen, bottom_row,
    nrow=3,
    ncol=1,
    rel_heights = c(0.05, 1.2, 0.8),
    labels = NULL
  )
}

if (is.na(args$width)) {
  final_width <- 9
} else {
  final_width <- args$width
}

if (is.na(args$height)) {
  final_height <- 9.1 + height_adj_methyl_binned + height_adj_cdr
} else {
  final_height <- args$height
}

ggsave(
  plot = plot,
  file = glue("{outdir}/{original_rname}_{mer_order}.tri.png"),
  height = final_height,
  width = final_width,
  dpi = 600
)
ggsave(
  plot = plot,
  file = glue("{outdir}/{original_rname}_{mer_order}.tri.pdf"),
  height = final_height,
  width = final_width
)
