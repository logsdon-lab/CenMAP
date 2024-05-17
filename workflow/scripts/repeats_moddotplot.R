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
  help = "HOR monomer order", default = "small"
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

df_rm_sat_out <- read_repeatmasker_sat_input(bed_sat_annot)
df_seq_ident <- read_bedpe(bed_seq_ident)
rname <- df_seq_ident$q[[1]]

plot <- make_cen_plot(
  rname,
  df_seq_ident,
  df_humas_hmmer_stv_out,
  df_rm_sat_out
)

# TODO: Scale with the number of HORs and length of the centromere.
ggsave(
  plot = plot,
  file = glue("{outdir}/{rname}_{mer_order}.tri.png"),
  height = 5,
  width = 9,
  dpi = 600
)
ggsave(
  plot = plot,
  file = glue("{outdir}/{rname}_{mer_order}.tri.pdf"),
  height = 5,
  width = 9
)
