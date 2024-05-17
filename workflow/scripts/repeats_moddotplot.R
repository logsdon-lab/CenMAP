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
require("argparse")


args <- commandArgs(trailingOnly = F)
scriptPath <- normalizePath(dirname(sub("^--file=", "", args[grep("^--file=", args)])))
if (length(scriptPath) == 0) {
  scriptPath <- "."
}
source(glue("workflow/scripts/r_utils/stained_glass.R"))
source(glue("workflow/scripts/r_utils/repeats.R"))

# source(glue("{scriptPath}/r_utils/stained_glass.R"))
# source(glue("{scriptPath}/r_utils/repeats.R"))

########################################################################################################
# EDIT THIS SECTION FOR YOUR INPUTS
# Rscript aln_plot_adjustedScale.R
# -b ./results/hgsvc_chr1.5000.5000.bed
# -p hgsvc_chr1.5000.5000.facet.equal
parser <- ArgumentParser()
parser$add_argument("-b", "--bed",  help="bedfile with alignment information")
parser$add_argument("-r", "--hor",  help="bedfile with HOR information")
parser$add_argument("-s", "--sat",  help="bedfile with satellite repeat information")
parser$add_argument("-o", "--output", default="results", help="Output plot.")
parser$add_argument("--mer_order",  help="HOR monomer order", default="small")
args <- parser$parse_args()

bed_seq_ident <- args$bed
# bed_seq_ident <- "HG00358_chr2_haplotype1-0000022:2-3544202.bed"
bed_sat_annot <- args$sat
# bed_sat_annot <- "HG00358_chr2_haplotype1-0000022:2-3544202_sat_annot.bed"
bed_hor_mon <- args$hor
# bed_hor_mon <- "AS-HOR-vs-HG00358_chr2_haplotype1-0000022:2-3544202_stv_row.bed"
outdir <- args$output
mer_order <- args$mer_order
dir.create(outdir)

df_humas_hmmer_stv_out <- read_one_humas_hmmer_input(bed_hor_mon)
df_humas_hmmer_stv_out <- switch(args$mer_order,
    # First sort by val. This sorts the dataframe but NOT the factor levels #larger HORs on top
    "large" = df_humas_hmmer_stv_out %>% arrange(mer),
    # First sort by val. This sorts the dataframe but NOT the factor levels #smaller HORs on top
    "small" = df_humas_hmmer_stv_out %>% arrange(-mer),
    stop(paste("Invalid mer reordering option:", args$mer_order))
  )

df_rm_sat_out <- read_repeatmasker_sat_input(bed_sat_annot)
df_seq_ident <- read_bedpe(bed_seq_ident)
rname <- df_seq_ident$q[[1]]

plot <- make_cen_plot(rname, df_seq_ident, df_humas_hmmer_stv_out, df_rm_sat_out)
ggsave(
  plot = plot,
  file = glue("{outdir}/{rname}_{mer_order}.tri.png"),
  height = 6,
  width = 10,
  dpi = 600
)
ggsave(
  plot = plot,
  file = glue("{outdir}/{rname}_{mer_order}.tri.pdf"),
  height = 6,
  width = 10
)
