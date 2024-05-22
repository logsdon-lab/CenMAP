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

df_rm_sat_out <- read_one_repeatmasker_sat_input(bed_sat_annot)
df_seq_ident <- read_bedpe(bed_seq_ident)
rname <- df_seq_ident$q[[1]]

plot <- make_cen_plot(
  rname,
  df_seq_ident,
  df_humas_hmmer_stv_out,
  df_rm_sat_out
)

# Adjust width and height based on number of HORs and contig length
if (is.na(args$width)) {
  num_uniq_mers <- df_humas_hmmer_stv_out %>% distinct(mer) %>% nrow()
  default_width <- 12
  width_adjustment <- ceiling(num_uniq_mers / 6) - 1
  final_width <- default_width + width_adjustment
} else {
  width_adjustment <- 0
  final_width <- args$width
}

if (is.na(args$height)) {
  contig_len_mb <- round(max(df_seq_ident$q_en) - min(df_seq_ident$q_st)) / 1000000
  default_height <- 5
  default_mb <- 4
  height_ctg_len_adjustment_factor <- round((contig_len_mb - default_mb) * 0.02)
  height_num_uniq_mers_adjustment_factor <- round(width_adjustment / 3)
  height_adjustment <- height_ctg_len_adjustment_factor + height_num_uniq_mers_adjustment_factor
  final_height <- default_height + height_adjustment
} else {
  final_height <- args$height
}

ggsave(
  plot = plot,
  file = glue("{outdir}/{rname}_{mer_order}.tri.png"),
  height = final_height,
  width = final_width,
  dpi = 600
)
ggsave(
  plot = plot,
  file = glue("{outdir}/{rname}_{mer_order}.tri.pdf"),
  height = final_height,
  width = final_width
)
