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
source(glue("{scriptPath}/r_utils/stained_glass.R"))

########################################################################################################
# EDIT THIS SECTION FOR YOUR INPUTS
# Rscript aln_plot_adjustedScale.R
# -b ./results/hgsvc_chr1.5000.5000.bed
# -p hgsvc_chr1.5000.5000.facet.equal
parser <- ArgumentParser()
parser$add_argument("-b", "--bed",  help="bedfile with alignment information")
parser$add_argument("-o", "--output", default="results", help="Output directory for the plots.")
parser$add_argument("-p", "--prefix", help="Prefix for the outputs")
args <- parser$parse_args()

PRE = args$prefix
GLOB = args$bed
OUT = args$output
OUT=glue("{OUT}/{PRE}_figures")
print(PRE)
print(GLOB)

DPI=600
#
# STOP EDITING
#
########################################################################################################
dir.create("results")
dir.create(OUT)
dir.create(glue("{OUT}/pdfs"))
dir.create(glue("{OUT}/pngs"))
all.files = Sys.glob(GLOB)
df = read_bedpe(all.files)
Qs = unique(df$q)
N=length(Qs)
columns = ceiling(sqrt(N+1))
rows = ceiling( (N+1) / columns)

scale = 2/3
plots = lapply(Qs, make_plots)
N <- length(plots)
columns <- ceiling(sqrt(N))
rows <- ceiling((N) / columns)
p = cowplot::plot_grid(plotlist = plots, nrow=rows, ncol=columns, labels = "auto");
ggsave(glue("{OUT}/pdfs/{PRE}.tri.all.pdf"), plot=p, height = 6*rows*scale, width = 6*columns)
ggsave(glue("{OUT}/pngs/{PRE}.tri.all.png"), plot=p, height = 6*rows*scale, width = 6*columns, dpi=DPI)

#
# big plot
#
facet_fig = cowplot::plot_grid(make_hist(df), make_dot(df), rel_heights = c(1,4), ncol=1)
ggsave(plot=facet_fig, file=glue("{OUT}/pdfs/{PRE}.facet.all.pdf"), height = 20, width = 16)
ggsave(plot=facet_fig, file=glue("{OUT}/pngs/{PRE}.facet.all.png"), height = 20, width = 16, dpi=DPI)
