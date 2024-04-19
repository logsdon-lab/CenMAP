library(ggplot2)
library(dplyr)
library(data.table)
library(stringr)
library(RColorBrewer)
library(argparser)
library(ggnewscale)
library(ggbeeswarm)


p <- arg_parser("Plot cumulative centromere HOR array lengths.")
p <- add_argument(
  p, "--input",
  help = "Input centromere HOR array lengths.", type = "character"
)
p <- add_argument(
  p, "--input_chm1",
  help = "Input chm1 centromere HOR array lengths.", type = "character"
)
p <- add_argument(
  p, "--input_chm13",
  help = "Input chm13 centromere HOR array lengths.", type = "character"
)
p <- add_argument(
  p, "--output",
  help = "Output plot file.", type = "character"
)
argv <- parse_args(p)

df_chm1_cen_lengths <- fread(
  argv$input_chm1,
  sep = "\t",
  col.names = c("chr", "start", "stop", "len")
)
df_chm13_cen_lengths <- fread(
  argv$input_chm13,
  sep = "\t",
  col.names = c("chr", "start", "stop", "len")
)
df_samples_cen_lengths <- fread(
  argv$input,
  sep = "\t",
  col.names = c("sm", "start", "stop", "len")
)

df_chm1_cen_lengths <- df_chm1_cen_lengths %>%
  mutate(Source = "CHM1") %>%
  group_by(chr, Source) %>%
  summarise(
    start = min(start),
    stop = max(stop),
    len = sum(len)
  )

df_chm13_cen_lengths <- df_chm13_cen_lengths %>%
  mutate(Source = "CHM13") %>%
  group_by(chr, Source) %>%
  summarise(
    start = min(start),
    stop = max(stop),
    len = sum(len),
  )

df_samples_adj_cen_lengths <- df_samples_cen_lengths %>%
  group_by(sm) %>%
  summarise(
    start = min(start),
    stop = max(stop),
    len = sum(len)
  ) %>%
  mutate(
    sm = str_replace(sm, "rc-chr", "chr"),
    Source = "Samples"
  ) %>%
  mutate(chr = str_extract(sm, "chr(\\d{1,2}|X|Y)")) %>%
  select(chr, start, stop, len, Source)

df_all_lengths <- rbind(
  df_chm1_cen_lengths,
  df_chm13_cen_lengths,
  df_samples_adj_cen_lengths
)
mean_length <- mean(df_all_lengths$len)

# Reorder chrs.
df_all_lengths$chr <- factor(
  df_all_lengths$chr,
  levels = c(paste0("chr", seq(22)), "chrX", "chrY")
)

chr_colors <- c(
  "#403E80",
  "#2C477D",
  "#2C477D",
  "#346B9B",
  "#3C80AA",
  "#4587A2",
  "#589D96",
  "#73ACA4",
  "#87BAAF",
  "#94BE9F",
  "#9FC38C",
  "#A1C27C",
  "#A4C165",
  "#C2C969",
  "#B5A957",
  "#DFD06C",
  "#F2D46C",
  "#E3C765",
  "#E5BA61",
  "#D89E56",
  "#C8824A",
  "#BB6B3E",
  "#AB5C40",
  "#9B4C41",
  "#8C3C42",
  "#731F37"
)

ggplot(df_all_lengths, aes(factor(chr), len)) +
  geom_violin(
    aes(factor(chr), len, fill = chr),
    scale = "width",
    alpha = 0.3,
    trim = FALSE,
    size = 0.8
  ) +
  scale_fill_manual(name = "Chromosome", values = chr_colors) +
  new_scale_fill() +
  geom_quasirandom(aes(factor(chr), len, color = chr)) +
  scale_color_manual(values = chr_colors) +
  # Remove legend for dotplot fill.
  guides(color = "none") +
  scale_fill_manual(values = chr_colors) +
  new_scale_fill() +
  # Add chm1 and chm13 dots with separate colorscale.
  geom_dotplot(
    data = df_chm1_cen_lengths,
    aes(factor(chr), len, fill = Source),
    binaxis = "y",
    stackdir = "center",
    dotsize = 0.4,
  ) +
  geom_dotplot(
    data = df_chm13_cen_lengths,
    aes(factor(chr), len, fill = Source),
    binaxis = "y",
    stackdir = "center",
    dotsize = 0.4,
  ) +
  scale_fill_manual(name = "Source", values = c("red", "black")) +
  # Add mean hline and label at end.
  geom_hline(yintercept = mean(df_all_lengths$len), linetype = 2) +
  annotate(
    "text",
    x = 23.5,
    y = mean_length + 200000,
    label = paste(round(mean_length / 1e6, 1), "Mbp")
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  ylab("Cumulative length of a-satellite HOR array(s) (Mbp)") +
  scale_y_continuous(
    breaks = seq(0, 10e6, 500000),
    labels = seq(0, 10, 0.5)
  )

ggsave(argv$output, width = 14, height = 6)
