#!/usr/bin/env Rscript
library(ggplot2)
library(dplyr)
require(scales)
library(RColorBrewer)
library(data.table)
library(cowplot)
library(glue)
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
args <- parser$parse_args()

# BED = args$bed
BED = "HG00358_chr2_haplotype1-0000022:2-3544202.bed"
# SAT = args$sat
SAT = "HG00358_chr2_haplotype1-0000022:2-3544202_sat_annot.bed"
# HOR = args$hor
HOR = "AS-HOR-vs-HG00358_chr2_haplotype1-0000022:2-3544202_stv_row.bed"
OUT = args$output
DPI=600
dir.create(OUT)

df = read_bedpe(BED)
Qs = df$q
N=length(Qs)
columns = ceiling(sqrt(N+1))
rows = ceiling( (N+1) / columns)
scale = 2/3

plots = make_plots(Qs)
N <- length(plots)
columns <- ceiling(sqrt(N))
rows <- ceiling((N) / columns)
p = cowplot::plot_grid(plotlist = plots, nrow=rows, ncol=columns, labels = "auto");
ggsave(glue("{OUT}/pdfs/{PRE}.tri.all.pdf"), plot=p, height = 6*rows*scale, width = 6*columns)
ggsave(glue("{OUT}/pngs/{PRE}.tri.all.png"), plot=p, height = 6*rows*scale, width = 6*columns, dpi=DPI)
