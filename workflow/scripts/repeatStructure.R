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


read_repeatmasker_sat_input <- function(input_file) {
  # read in BED file
  df <- fread(
    input_file,
    select = c(1:9),
    stringsAsFactors = TRUE,
    fill = TRUE,
    sep = "\t",
    quote = "",
    header = FALSE
  )
  cols <- c("chr", "start", "stop", "region", "value", "strand", "start2", "stop2", "rgb")
  colnames(df) <- cols

  # reorder rows so that live arrays are plotted on top
  df$region <- factor(df$region, levels = c("ct", "asat", "bsat", "gsat", "hsat1A", "hsat1B", "hsat2", "hsat3"), ordered = T)

  return(df)
}

get_humas_hmmer_sf_annot_colors <- function() {
  myColors <- c(
    "#A8275C", "#9AC78A", "#A53D63", "#3997C6", "#29A3CE", "#5EB2A7", "#38A49B", "#45B4CE", "#A53D63", "#AA1B63", "#3F66A0",
    "#D66C54", "#BFDD97", "#C0D875", "#E5E57A", "#B75361", "#F9E193", "#C6625D", "#E5D1A1", "#A1B5E5", "#9F68A5", "#81B25B",
    "#F4DC78", "#7EC0B3", "#A8275C", "#8CC49F", "#893F89", "#6565AA"
  )
  names(myColors) <- levels(as.factor(c(
    "1", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19",
    "2", "20", "21", "22", "24", "26", "3", "30", "32", "34", "35",
    "4", "5", "6", "7", "8", "9"
  )))
  return(myColors)
}

get_rm_sat_annot_colors <- function() {
  myColors <- c(
    "#58245B",
    "#3A3A3A",
    "#DDDDDD", # lighter gray
    # "#C2C2C2", #darker gray
    "#3A3A3A",
    "#3A3A3A", "#3A3A3A",
    "#3A3A3A", "#3A3A3A"
  )
  names(myColors) <- levels(as.factor(c(
    "asat",
    "bsat",
    "ct",
    "gsat",
    "hsat1A", "hsat1B",
    "hsat2", "hsat3"
  )))
  return(myColors)
}


read_humas_hmmer_input <- function(
    input_chr,
    input_chm1,
    input_chm13,
    chr_name,
    mer_order = "large",
    hor_filter = 0) {
  cols_to_take <- seq(5)
  monomer_len <- 170
  cols <- c("chr", "start", "stop", "hor", "strand")

  chm13 <- fread(input_chm13,
    sep = "\t",
    stringsAsFactors = TRUE,
    fill = TRUE, quote = "",
    header = FALSE, select = cols_to_take
  )
  chm1 <- fread(input_chm1,
    sep = "\t",
    stringsAsFactors = TRUE,
    fill = TRUE, quote = "",
    header = FALSE, select = cols_to_take
  )
  hgsvc3 <- fread(input_chr,
    sep = "\t",
    stringsAsFactors = TRUE,
    fill = TRUE, quote = "",
    header = FALSE, select = cols_to_take
  )

  colnames(chm13) <- cols
  colnames(chm1) <- cols
  colnames(hgsvc3) <- cols

  # select the chr
  chm13_select <- chm13 %>%
    filter(str_detect(chr, paste0(chr_name, "$")))
  chm1_select <- chm1 %>%
    filter(str_detect(chr, paste0(chr_name, "$")))

  # combine the CHM1 and HGSVC3 centromeres
  chm13_select$chr <- gsub("chr", "chm13_chr", chm13_select$chr)
  chm1_select$chr <- gsub("chr", "chm1_chr", chm1_select$chr)
  chm13_chm1_hgsvc3 <- rbind(chm13_select, chm1_select, hgsvc3)

  # determine distance between start and stop
  chm13_chm1_hgsvc3$length <- chm13_chm1_hgsvc3$stop - chm13_chm1_hgsvc3$start

  # calculate monomer size and round
  chm13_chm1_hgsvc3$mer <- as.numeric(round(chm13_chm1_hgsvc3$length / monomer_len))

  # filter monomers
  chm13_chm1_hgsvc3 <- switch(chr_name,
    "chr10" = subset(chm13_chm1_hgsvc3, as.numeric(mer) >= 5),
    "chr20" = subset(chm13_chm1_hgsvc3, as.numeric(mer) >= 5),
    "chrY" = subset(chm13_chm1_hgsvc3, as.numeric(mer) >= 30),
    "chr17" = subset(chm13_chm1_hgsvc3, as.numeric(mer) >= 4),
    chm13_chm1_hgsvc3
  )

  # filter for HORs that occur at least 20 times (10 times per haplotype)
  df_mer <- chm13_chm1_hgsvc3 %>%
    group_by(mer) %>%
    filter(n() > hor_filter) #

  # Add contig length to start and stop so aligns with repeatmasker annotations.
  df_final <- df_mer %>%
    separate_wider_delim(chr, delim = ":", names = c("chr", "range"), too_few = "align_start") %>%
    separate_wider_delim(range, delim = "-", names = c("ctg_start", "ctg_stop"), too_few = "align_start") %>%
    mutate(
      ctg_start = replace_na(as.numeric(ctg_start), 0),
      ctg_stop = replace_na(as.numeric(ctg_stop), 0)
    ) %>%
    mutate(
      start = ctg_start + start,
      stop = ctg_start + stop
    ) %>%
    mutate(new_chr = str_extract(chr, "([\\w_-]*?):", group = 1))

  # reorder data to put larger or smaller -mers on top
  df_final <- switch(mer_order,
    # First sort by val. This sorts the dataframe but NOT the factor levels #larger HORs on top
    "large" = df_final %>% arrange(mer),
    # First sort by val. This sorts the dataframe but NOT the factor levels #smaller HORs on top
    "small" = df_final %>% arrange(-mer),
    stop(paste("Invalid mer reordering option:", mer_order))
  )

  return(df_final)
}

# Create a parser
p <- arg_parser("Plot combined SF + satellites")
p <- add_argument(p, "--input_rm_sat", help = "Input bed file made from annotated satellite RM output.", type = "character")

p <- add_argument(p, "--input_sf",
  help = "Input HGSVC3 HumAS-HMMER formatted output.",
  type = "character"
)
p <- add_argument(p, "--input_sf_chm13",
  help = "Input CHM13 HumAS-HMMER formatted output",
  type = "character"
)
p <- add_argument(p, "--input_sf_chm1",
  help = "Input CHM1 HumAS-HMMER formatted output",
  type = "character"
)
p <- add_argument(p, "--chr",
  help = "Chromosome to plot. ex. chrX",
  type = "character"
)
p <- add_argument(p, "--output",
  help = "Output plot file. Defaults to {chr}_hgsvc3_{order_hor}erontop.png.",
  type = "character", default = NA
)
p <- add_argument(p, "--hor_filter",
  help = "Filter for HORs that occur at least 20 times (10 times per haplotype)",
  type = "numeric", default = 0
)
p <- add_argument(p, "--mer_order",
  help = "Reorder data to put 'large'-r or 'small'-r -mers on top.",
  type = "character", default = "large"
)

argv <- parse_args(p)

df_rm_sat_out <- read_repeatmasker_sat_input(argv$input_rm_sat)
df_humas_hmmer_sf_out <- read_humas_hmmer_input(argv$input_sf, argv$input_sf_chm1, argv$input_sf_chm13, argv$chr, argv$mer_order, argv$hor_filter) %>%
  filter(str_detect(chr, "^chm", negate = TRUE))

# Set new minimum and standardize scales.
new_min <- min(df_rm_sat_out$start2)
df_rm_sat_out <- df_rm_sat_out %>%
  join(
    df_rm_sat_out %>% group_by(chr) %>% summarize(dst_diff = min(start2) - new_min),
    by = "chr"
  ) %>%
  filter(!dst_diff == 0) %>%
  group_by(chr) %>%
  arrange(start2) %>%
  mutate(
    start2 = start2 - dst_diff,
    stop2 = stop - dst_diff
  )

df_humas_hmmer_sf_out <- df_humas_hmmer_sf_out %>%
  join(
    df_rm_sat_out %>% group_by(chr) %>% summarize(dst_diff = min(start) - new_min),
    by = "chr"
  ) %>%
  filter(!dst_diff == 0) %>%
  group_by(chr) %>%
  arrange(start) %>%
  mutate(
    start = start - dst_diff,
    stop = stop - dst_diff
  )

height <- 5
ggplot(data = df_rm_sat_out[order(df_rm_sat_out$region), ]) +
  geom_segment(aes(x = start2, y = chr, xend = stop2 + 1000, yend = chr, color = region), alpha = 1, size = height * 2) +
  scale_color_manual(values = get_rm_sat_annot_colors()) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 12)) +
  # New colorscale.
  new_scale_color() +
  geom_segment(data = df_humas_hmmer_sf_out, aes(x = start, xend = stop + 2000, y = chr, yend = chr, color = as.factor(mer)), linewidth = 9) +
  scale_color_manual(values = get_humas_hmmer_sf_annot_colors()) +
  theme_classic() +
  theme(axis.line.y = element_line(colour = "white")) +
  theme(legend.position = "bottom") +
  theme(legend.key.width = unit(1, "cm")) +
  theme(legend.title = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 6))) +
  xlab("Position (kbp)")


if (!interactive()) pdf(NULL)
height <- length(unique(df_humas_hmmer_sf_out$chr)) * 0.5
ggsave(argv$output, width = 10, height = height + 4)
