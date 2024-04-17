library(ggplot2)
library(dplyr)
library(data.table)
library(stringr)
library(RColorBrewer)
library(argparser)


p <- arg_parser("Plot cumulative centromere HOR array lengths.")
p <- add_argument(p, "--input", help = "Input centromere HOR array lengths.", type = "character")
p <- add_argument(p, "--input_chm1", help = "Input chm1 centromere HOR array lengths.", type = "character")
p <- add_argument(p, "--input_chm13", help = "Input chm13 centromere HOR array lengths.", type = "character")
p <- add_argument(p, "--output", help = "Output plot file.", type = "character")
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
# At this point, "rc_chr" replaced with "rc-chr" and '_' replaced with '\t'
df_samples_cen_lengths <- fread(
  argv$input,
  sep = "\t",
  col.names = c("sm", "chr", "ctg", "start", "stop", "len")
)

df_chm1_cen_lengths <- df_chm1_cen_lengths %>%
  mutate(Source = "CHM1")

df_chm13_cen_lengths <- df_chm13_cen_lengths %>%
  mutate(Source = "CHM13")

df_samples_adj_cen_lengths <- df_samples_cen_lengths %>%
  mutate(
    chr = str_replace(chr, "rc-chr", "chr"),
    hap = str_extract(ctg, "haplotype\\d", group = 0),
    Source = "Samples"
  )


df_all_lengths <- rbind(
  df_chm1_cen_lengths,
  df_chm13_cen_lengths,
  df_samples_adj_cen_lengths %>%
    select(chr, start, stop, len, Source)
)
# Reorder chrs.
df_all_lengths$chr <- factor(
  df_all_lengths$chr,
  levels = c(paste0("chr", seq(22)), "chrX", "chrY")
)

# TODO: Group by chr and source and calculate cumlen.

ggplot(df_all_lengths, aes(factor(chr), len, fill = Source)) +
  geom_violin(
    fill = "grey",
    scale = "width",
  ) +
  ## TODO: Get colorscale.
  # +scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))
  geom_dotplot(
    binaxis = "y",
    stackdir = "center",
    dotsize = 0.2,
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank()
  ) +
  ylab("Centromere Length") +
  scale_y_continuous(
    breaks = seq(0, 6e6, 1e6),
    labels = formatC(seq(0, 6e6, 1e6), format = "d", big.mark = ",")
  )
