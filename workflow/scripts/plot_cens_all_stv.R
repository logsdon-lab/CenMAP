#!/usr/bin/env Rscript
library(argparser, quietly = TRUE)
library(ggplot2)
library(plyr)
library(ggnewscale)
require(gridExtra)
require(scales)
library(RColorBrewer)
library(reshape2)
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(ggpubr)
library(glue)


args <- commandArgs(trailingOnly = FALSE)
wd <- normalizePath(dirname(sub("^--file=", "", args[grep("^--file=", args)])))
source(glue("{wd}/r_utils/repeats.R"))


# Create a parser
p <- arg_parser("Plot combined SF + satellites")
p <- add_argument(p, "--input_rm_sat",
  help = "Input bed file made from annotated satellite RM output.",
  type = "character"
)
p <- add_argument(p, "--input_stv",
  help = "Input sample HumAS-HMMER formatted output.",
  type = "character"
)
p <- add_argument(p, "--input_stv_chm13",
  help = "Input CHM13 HumAS-HMMER formatted output",
  type = "character"
)
p <- add_argument(p, "--input_stv_chm1",
  help = "Input CHM1 HumAS-HMMER formatted output",
  type = "character"
)
p <- add_argument(p, "--chr",
  help = "Chromosome to plot. ex. chrX",
  type = "character"
)
p <- add_argument(p, "--output",
  help = "Output plot.",
  type = "character", default = NA
)
p <- add_argument(p, "--output_dir",
  help = "Output dir with each centromere split into a separate image.",
  type = "character", default = NA
)
p <- add_argument(p, "--hor_filter",
  help = "Filter for HORs that occur at least n times.",
  type = "numeric", default = 0
)
p <- add_argument(p, "--mer_order",
  help = "Reorder data to put 'large'-r or 'small'-r -mers on top.",
  type = "character", default = "large"
)

argv <- parse_args(p)

df_rm_sat_out <- read_multiple_repeatmasker_sat_input(argv$input_rm_sat)
# Filter out chr without annotations.
df_humas_hmmer_stv_out <- read_multiple_humas_hmmer_input(
  argv$input_stv,
  argv$input_stv_chm1,
  argv$input_stv_chm13,
  argv$chr,
  argv$hor_filter
) %>%
  filter(chr %in% df_rm_sat_out$chr)
df_rm_sat_out <- df_rm_sat_out %>% filter(chr %in% df_humas_hmmer_stv_out$chr)

# Set new minimum and standardize scales.
new_min <- min(df_rm_sat_out$start2)
df_rm_sat_out <- df_rm_sat_out %>%
  join(
    df_rm_sat_out %>%
      group_by(chr) %>%
      summarize(dst_diff = min(start2) - new_min),
    by = "chr"
  ) %>%
  group_by(chr) %>%
  arrange(start2) %>%
  mutate(
    start2 = start2 - dst_diff,
    stop2 = stop - dst_diff
  )

df_humas_hmmer_stv_out <- df_humas_hmmer_stv_out %>%
  join(
    df_rm_sat_out %>%
      group_by(chr) %>%
      summarize(dst_diff = min(start) - new_min),
    by = "chr"
  ) %>%
  group_by(chr) %>%
  arrange(start) %>%
  mutate(
    start = start - dst_diff,
    stop = stop - dst_diff
  ) %>%
  ungroup(chr)

# First sort by val. This sorts the dataframe but NOT the factor levels
# larger HORs on top
# smaller HORs on top
df_humas_hmmer_stv_out <- switch(argv$mer_order,
    "large" = df_humas_hmmer_stv_out %>% arrange(mer),
    "small" = df_humas_hmmer_stv_out %>% arrange(-mer),
    stop(paste("Invalid mer reordering option:", argv$mer_order))
  )

# Generate split plot.
if (!is.na(argv$output_dir)) {
  dir.create(argv$output_dir)

  for (ctg in unique(df_rm_sat_out$chr)) {
    plt_ctg <- plot_single_ctg(ctg, df_rm_sat_out, df_humas_hmmer_stv_out)

    ggsave(
      paste0(argv$output_dir, "/", ctg, ".png"),
      plot = plt_ctg,
      width = 12,
      height = 1
    )
    while (!is.null(dev.list())) dev.off()
  }
}

# Generate full plot.
plt <- plot_all_ctgs(df_rm_sat_out, df_humas_hmmer_stv_out)

# Scale height to fit number contigs.
height <- length(unique(df_humas_hmmer_stv_out$chr)) * 0.5
ggsave(
  argv$output,
  plot = plt,
  width = 14,
  height = height + 4,
  limitsize = FALSE
)

# And the legend.
if (!is.na(argv$output_dir)) {
  legend <- get_legend(plt)

  # Convert to a ggplot and print
  plt_legend <- as_ggplot(legend)

  ggsave(
    paste0(argv$output_dir, "/legend.png"),
    plot = plt_legend,
    width = 8,
    height = 4
  )
}
