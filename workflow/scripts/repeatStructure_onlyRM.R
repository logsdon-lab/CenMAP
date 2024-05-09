#!/usr/bin/env Rscript
# 231029 repeatStructure_onlyRM.R
library(argparser, quietly = TRUE)
library(ggplot2)
library(directlabels)
library(plyr)
require(gridExtra)
require(scales)
library(RColorBrewer)
library(reshape2)
library(data.table)
library(dplyr)
library(tidyr)

# Create a parser
p <- arg_parser("Plot ALR region from RepeatMasker output.")
p <- add_argument(p, "input", help = "Input repeatmasker.out file", type = "character")
p <- add_argument(p, "output", help = "Output plot file.", type = "character")
argv <- parse_args(p)

# setwd("~/Documents/Eichler Lab/Weekly plans/HGSVC3/RepeatMasker")
df <- fread(argv$input, select = c(1:15), stringsAsFactors = FALSE, fill = TRUE, quote = "", header = FALSE, skip = 2)
cols <- c("idx", "div", "deldiv", "insdiv", "contig", "start", "end", "left", "C", "type", "rClass", "right", "x", "y", "z")
colnames(df) <- cols

# filter by length
df2 <- df %>%
  group_by(contig) %>%
  filter(last(end) > 0000)

# trim length
# df2 = subset(df2, df2$start>400000)
# df2 = subset(df2, df2$start<7900000)

# rename type as rClass
{
  mask <- (df2$rClass == "Satellite/centr") | (df2$rClass == "Satellite")
  df2$rClass[mask] <- df2$type[mask]
  df2$rClass <- sub("/ERVK", "", df2$rClass)
  df2$rClass <- sub("/ERVL", "", df2$rClass)
  df2$rClass <- sub("/ERV1", "", df2$rClass)
  df2$rClass <- sub("/CR1", "", df2$rClass)
  df2$rClass <- sub("/L1", "", df2$rClass)
  df2$rClass <- sub("/L2", "", df2$rClass)
  df2$rClass <- sub("/RTE-X", "", df2$rClass)
  df2$rClass <- sub("/RTE-BovB", "", df2$rClass)
  df2$rClass <- sub("/Gypsy", "", df2$rClass)
  df2$rClass <- sub("-MaLR", "", df2$rClass)
  df2$rClass <- sub("/Alu", "", df2$rClass)
  df2$rClass <- sub("/Deu", "", df2$rClass)
  df2$rClass <- sub("/MIR", "", df2$rClass)
  df2$rClass <- sub("?", "", df2$rClass)
  df2$rClass <- sub("/hAT", "", df2$rClass)
  df2$rClass <- sub("/hAT-Blackjack", "", df2$rClass)
  df2$rClass <- sub("/hAT-Charlie", "", df2$rClass)
  df2$rClass <- sub("/MULE-MuDR", "", df2$rClass)
  df2$rClass <- sub("/PiggyBac", "", df2$rClass)
  df2$rClass <- sub("/TcMar-Mariner", "", df2$rClass)
  df2$rClass <- sub("/TcMar", "", df2$rClass)
  df2$rClass <- sub("/TcMar?", "", df2$rClass)
  df2$rClass <- sub("/hAT-Tip100", "", df2$rClass)
  df2$rClass <- sub("/TcMar-Tigger", "", df2$rClass)
  df2$rClass <- sub("/Dong-R4", "", df2$rClass)
  df2$rClass <- sub("/tRNA", "", df2$rClass)
  df2$rClass <- sub("DNA-Tc2", "DNA", df2$rClass)
  df2$rClass <- sub("DNA?", "DNA", df2$rClass)
  df2$rClass <- sub("DNA-Blackjack", "DNA", df2$rClass)
  df2$rClass <- sub("DNA-Charlie", "DNA", df2$rClass)
  df2$rClass <- sub("DNA-Tigger", "DNA", df2$rClass)
  df2$rClass <- sub("DNA-Tip100", "DNA", df2$rClass)
  df2$rClass <- sub("GSATX", "GSAT", df2$rClass)
  df2$rClass <- sub("LTR\\S", "LTR", df2$rClass)
  df2$type <- as.factor(df2$type)

  df2$start2 <- df2$start
  df2$end2 <- df2$end
  mask <- df2$C == "C"
  df2$start2[mask] <- df2$end[mask]
  df2$end2[mask] <- df2$start[mask]
}

# rename contigs for easy comparison of structure
df2$contig <- gsub("chm1_cen", "chr", df2$contig)

# rename certain repeat types
{
  df2$rClass <- sub("SAR", "HSat1A", df2$rClass)
  df2$rClass <- sub("HSAT", "HSat1B", df2$rClass)
  df2$rClass <- sub("HSATII", "HSat2", df2$rClass)
  df2$rClass <- sub("(CATTC)n", "HSat2", df2$rClass)
  df2$rClass <- sub("(GAATG)n", "HSat2", df2$rClass)
}

#### assign colors to rClassses

# centromere colors (original)
{
  myColors <- c(
    "#522758", "#5188c9",
    "#4ea836", "#FFFFFF",
    "#FFFFFF", "#FFFFFF",
    "#FFFFFF", "#FFFFFF",
    "#FFFFFF", "#FFFFFF",
    "#9370af", "#9370af",
    "#84bac6", "#84bac6",
    "#f2b80b", "#ad8c2a",
    "#E571AB", "#E571AB", "#E571AB",
    "#FFFFFF", "#FFFFFF",
    "#e1665e", "#FFFFFF",
    "#FFFFFF", "#FFFFFF",
    "#FFFFFF", "#FFFFFF",
    "#FFFFFF", "#FFFFFF",
    "#FFFFFF", "#508f52",
    "#47b1b5", "#47b1b5",
    "#47b1b5", "#FFFFFF",
    "#FFFFFF", "#FFFFFF",
    "#FFFFFF", "#FFFFFF",
    "#FFFFFF", "#FFFFFF", "#FFFFFF",
    "#FFFFFF", "#FFFFFF", "#FFFFFF",
    "#FFFFFF", "#FFFFFF", "#FFFFFF",
    "#FFFFFF", "#FFFFFF", "#FFFFFF",
    "#ad8c2a", "#136aaf"
  )
  names(myColors) <- levels(as.factor(c("ALR/Alpha", "BSR/Beta", "CER", "DNA", "DNA-Ac", "DNA-Tag1", "DNA-Tc1", "DNA?", "DNA/Merlin", "DNA/PIF-Harbinger", "GSAT", "GSATII", "HSat1A", "HSat1B", "HSat2", "HSat3", "LINE", "LINE-Tx1", "LINE/Penelope", "Low_complexity", "LSAU", "LTR", "RC/Helitron", "Retroposon/SVA", "rRNA", "Satellite/acro", "Satellite/subtelo", "SATR1", "SATR2", "scRNA", "Simple_repeat", "SINE", "SINE-RTE", "SINE/5S-Deu-L2", "snRNA", "srpRNA", "SST1", "tRNA")))
}

height <- 5
ggplot() +
  geom_segment(data = df2, aes(x = start2 / 1000, y = contig, xend = (end2 / 1000) + 1, yend = contig, color = rClass), size = height * 2) +
  # geom_segment(data=df2, aes(x=start2/1000, y=contig, xend=(end2/1000)+50, yend=contig), size=.05, arrow=arrow(angle=40, length=unit(0.02, "npc"))) +
  # geom_segment(data=df2, aes(x=86095438/1000, y=contig, xend=(86102489/1000), yend=contig), size=.1) +
  # geom_segment(data=df2, aes(x=86102270/1000, y=contig, xend=(86109326/1000), yend=contig), size=.1) +
  scale_color_manual(values = myColors) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 12)) +
  theme_classic() +
  theme(axis.line.y = element_line(colour = "white")) +
  theme(legend.position = "bottom") +
  theme(legend.key.width = unit(1, "cm")) +
  theme(legend.title = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 1))) +
  xlab("Position (kbp)")

if (!interactive()) pdf(NULL)
height <- length(unique(df$contig)) * 0.5
ggsave(argv$output, width = 14, height = height + 4, limitsize = FALSE)
