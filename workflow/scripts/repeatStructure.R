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
  # Filter duplicated chm13 rows.
  df <- df %>%
    mutate(chr=str_replace(chr, "cen", "chr")) %>%
    filter(!str_detect(chr, "^chr[0-9XY]+$")) %>%
    # Correct for version and different naming of chr. ex. chm1_cen1v8 -> chm1_chr1
    mutate(chr=str_remove(chr, "v\\d+"))

  # reorder rows so that live arrays are plotted on top
  df$region <- factor(df$region, levels = c("ct", "asat", "bsat", "gsat", "hsat1A", "hsat1B", "hsat2", "hsat3"), ordered = T)

  return(df)
}

get_humas_hmmer_stv_annot_colors <- function() {
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
  samples <- fread(input_chr,
    sep = "\t",
    stringsAsFactors = TRUE,
    fill = TRUE, quote = "",
    header = FALSE, select = cols_to_take
  )

  colnames(chm13) <- cols
  colnames(chm1) <- cols
  colnames(samples) <- cols

  # combine the CHM1 and samples centromeres
  chm13$chr <- gsub("chr", "chm13_chr", chm13$chr)
  chm1$chr <- gsub("chr", "chm1_chr", chm1$chr)

  # Requires coordinates in name. ex. chr?:1-2
  chm13_chm1_samples <- rbind(chm13, chm1, samples)

  # determine distance between start and stop
  chm13_chm1_samples$length <- chm13_chm1_samples$stop - chm13_chm1_samples$start

  # calculate monomer size and round
  chm13_chm1_samples$mer <- as.numeric(round(chm13_chm1_samples$length / monomer_len))

  # filter monomers
  chm13_chm1_samples <- switch(chr_name,
    "chr10" = subset(chm13_chm1_samples, as.numeric(mer) >= 5),
    "chr20" = subset(chm13_chm1_samples, as.numeric(mer) >= 5),
    "chrY" = subset(chm13_chm1_samples, as.numeric(mer) >= 30),
    "chr17" = subset(chm13_chm1_samples, as.numeric(mer) >= 4),
    chm13_chm1_samples
  )

  # filter for HORs that occur at least 20 times (10 times per haplotype)
  df_stv <- chm13_chm1_samples %>%
    group_by(mer) %>%
    filter(n() > hor_filter)

  # Fix orientation.
  df_rc_stv <- df_stv %>%
    filter(str_detect(chr, "rc")) %>%
    mutate(
      ctg_start=as.integer(replace_na(str_extract(chr, ":(\\d+)-", 1), 0)),
      ctg_stop=as.integer(replace_na(str_extract(chr, "-(\\d+)$", 1), 0))
    ) %>%
    mutate(
      new_start=ctg_start+abs(ctg_stop-stop),
      new_stop=ctg_start+abs(ctg_stop-start)
    ) %>%
    mutate(start=new_start, stop=new_stop) %>%
    select(chr, start, stop, hor, strand, length, mer)

  df_stv <- df_stv %>%
    filter(!str_detect(chr, "rc"))
  df_both <- rbind(df_rc_stv, df_stv)

  return(df_both)
}

# Create a parser
p <- arg_parser("Plot combined SF + satellites")
p <- add_argument(p, "--input_rm_sat", help = "Input bed file made from annotated satellite RM output.", type = "character")

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
  help = "Filter for HORs that occur at least 20 times (10 times per haplotype)",
  type = "numeric", default = 0
)
p <- add_argument(p, "--mer_order",
  help = "Reorder data to put 'large'-r or 'small'-r -mers on top.",
  type = "character", default = "large"
)

argv <- parse_args(p)
# test_chr <- "chr9"
# {
#   argv$input_rm_sat <- paste0("results/repeatmasker_sat_annot/repeats/all_cens_", test_chr, ".annotation.fa.out")
#   argv$input_stv <- paste0("results/hor_stv/bed/", test_chr, "_AS-HOR_stv_row.all.bed")
#   argv$input_stv_chm13 <- "data/annotations/AS-HOR-vs-chm13_cens_v18.correctcoords.stv_row.all2.bed"
#   argv$input_stv_chm1 <- "data/annotations/AS-HOR-vs-chm1_cens_v21.stv_row.all2.bed"
#   argv$output <- "test.png"
#   argv$chr <- test_chr
#   argv$mer_order <- "large"
# }
# {
#   input_file <- argv$input_rm_sat
#   input_chr <- argv$input_stv
#   input_chm1 <- argv$input_stv_chm1
#   input_chm13 <- argv$input_stv_chm13
#   chr_name <- test_chr
# }

df_rm_sat_out <- read_repeatmasker_sat_input(argv$input_rm_sat)
# Filter out chr without annotations.
df_humas_hmmer_stv_out <- read_humas_hmmer_input(argv$input_stv, argv$input_stv_chm1, argv$input_stv_chm13, argv$chr, argv$hor_filter) %>%
  filter(chr %in% df_rm_sat_out$chr)
df_rm_sat_out <- df_rm_sat_out %>% filter(chr %in% df_humas_hmmer_stv_out$chr)

# Set new minimum and standardize scales.
new_min <- min(df_rm_sat_out$start2)
df_rm_sat_out <- df_rm_sat_out %>%
  join(
    df_rm_sat_out %>% group_by(chr) %>% summarize(dst_diff = min(start2) - new_min),
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
    df_rm_sat_out %>% group_by(chr) %>% summarize(dst_diff = min(start) - new_min),
    by = "chr"
  ) %>%
  group_by(chr) %>%
  arrange(start) %>%
  mutate(
    start = start - dst_diff,
    stop = stop - dst_diff
  ) %>%
  ungroup(chr)

df_humas_hmmer_stv_out <- switch(argv$mer_order,
    # First sort by val. This sorts the dataframe but NOT the factor levels #larger HORs on top
    "large" = df_humas_hmmer_stv_out %>% arrange(mer),
    # First sort by val. This sorts the dataframe but NOT the factor levels #smaller HORs on top
    "small" = df_humas_hmmer_stv_out %>% arrange(-mer),
    stop(paste("Invalid mer reordering option:", argv$mer_order))
  )

height <- 5

# Generate split plot.
if (!is.null(argv$output_dir)) {
  dir.create(argv$output_dir)

  for (ctg in unique(df_rm_sat_out$chr)) {
    ggplot(data = df_rm_sat_out[order(df_rm_sat_out$region), ] %>% filter(chr == ctg)) +
      geom_segment(
        aes(x = start2, y = chr, xend = stop2 + 1000, yend = chr, color = region),
        alpha = 1,
        size = height * 2
      ) +
      scale_color_manual(values = get_rm_sat_annot_colors()) +
      scale_x_continuous(breaks = scales::pretty_breaks(n = 12)) +
      # New colorscale.
      new_scale_color() +
      geom_segment(
        data = df_humas_hmmer_stv_out %>% filter(chr == ctg),
        aes(x = start, xend = stop + 2000, y = chr, yend = chr, color = as.factor(mer)),
        size = height * 2
      ) +
      scale_color_manual(values = get_humas_hmmer_stv_annot_colors()) +
      theme_classic() +
      # Remove everything.
      theme(
        axis.line.x = element_line(colour = "white"),
        axis.line.y = element_line(colour = "white"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none"
      ) +
      guides(color = guide_legend(override.aes = list(size = 6))) +
      xlab("Position (kbp)")

    if (!interactive()) pdf(NULL)
    ggsave(paste0(argv$output_dir, ctg, ".png"), width = 12, height = 1)
  }
}

# Generate full plot.
plt <- ggplot(data = df_rm_sat_out[order(df_rm_sat_out$region), ]) +
  geom_segment(
    aes(x = start2, y = chr, xend = stop2 + 1000, yend = chr, color = region),
    alpha = 1,
    size = height * 2
  ) +
  scale_color_manual(values = get_rm_sat_annot_colors()) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 12)) +
  # New colorscale.
  new_scale_color() +
  geom_segment(
    data = df_humas_hmmer_stv_out,
    aes(x = start, xend = stop + 2000, y = chr, yend = chr, color = as.factor(mer)),
    size = height * 2
  ) +
  scale_color_manual(values = get_humas_hmmer_stv_annot_colors()) +
  theme_classic() +
  theme(axis.line.y = element_line(colour = "white")) +
  theme(legend.position = "bottom") +
  theme(legend.key.width = unit(1, "cm")) +
  theme(legend.title = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 6))) +
  xlab("Position (kbp)")


if (!interactive()) pdf(NULL)
# Scale height to fit number contigs.
height <- length(unique(df_humas_hmmer_stv_out$chr)) * 0.5
ggsave(argv$output, plot = plt, width = 14, height = height + 4)

# And the legend.
if (!is.null(argv$output_dir)) {
  legend <- get_legend(plt)

  # Convert to a ggplot and print
  plt_legend <- as_ggplot(legend)

  ggsave(paste0(argv$output_dir, "/legend.png"), plot = plt_legend, width = 8, height = 4)
}
