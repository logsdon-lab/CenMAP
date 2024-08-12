# 230816 repeatStructure.R
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
library(argparser, quietly = TRUE)

p <- arg_parser("Plot satellite RepeatMasker output.")
p <- add_argument(p, "input", help = "Input bed file made from annotated satellite RM output.", type = "character")
p <- add_argument(p, "output", help = "Output plot file.", type = "character")
argv <- parse_args(p)

# read in BED file
df <- fread(argv$input, select = c(1:9), stringsAsFactors = TRUE, fill = TRUE, quote = "", header = FALSE)
cols <- c("chr", "start", "stop", "region", "value", "strand", "start2", "stop2", "rgb")
colnames(df) <- cols

# adjust start2 and stop2
df <- df %>%
  group_by(chr) %>%
  mutate(start2 = start - min(start)) %>%
  mutate(stop2 = stop - min(start))

# reorder rows so that live arrays are plotted on top
df$region <- factor(df$region, levels = c("ct", "asat", "bsat", "gsat", "hsat1A", "hsat1B", "hsat2", "hsat3"), ordered = T)

# centromere colors
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


height <- 5
ggplot(data = df[order(df$region), ]) +
  geom_segment(aes(x = start2, y = chr, xend = stop2 + 1000, yend = chr, color = region), alpha = 1, size = height * 2) +
  scale_color_manual(values = myColors) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 12)) +
  theme_classic() +
  theme(axis.line.y = element_line(colour = "white")) +
  theme(legend.position = "bottom") +
  theme(legend.key.width = unit(1, "cm")) +
  theme(legend.title = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 6))) +
  xlab("Position (kbp)")


if (!interactive()) pdf(NULL)
ggsave(argv$output, width = 10, height = 22)
