#!/usr/bin/env Rscript
# 230303_StVHORorganization_simplified_2datasets_basedOnStartStop_newcolors.R

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
library(stringr)
library(argparser, quietly = TRUE)


p <- arg_parser("Plot STV from formatted HumAS-HMMER output.")
p <- add_argument(p, "--input",
                  help="Input HGSVC3 HumAS-HMMER formatted output.",
                  type="character")
p <- add_argument(p, "--input_chm13",
                  help="Input CHM13 HumAS-HMMER formatted output",
                  type="character")
p <- add_argument(p, "--input_chm1",
                  help="Input CHM1 HumAS-HMMER formatted output",
                  type="character")
p <- add_argument(p, "--chr",
                  help="Chromosome to plot. ex. chrX",
                  type="character")
p <- add_argument(p, "--output",
                  help="Output plot file. Defaults to {chr}_hgsvc3_{order_hor}erontop.png.",
                  type="character", default=NA)
p <- add_argument(p, "--output_tbl",
                  help="Output table file. Defaults to {chr}_length.tsv.",
                  type="character", default=NA)
p <- add_argument(p, "--hor_filter",
                   help="Filter for HORs that occur at least 20 times (10 times per haplotype)",
                   type="numeric", default=0)
p <- add_argument(p, "--order_hor",
                  help="Reorder data to put 'large'-r or 'small'-r -mers on top.",
                  type="character", default='large')
p <- add_argument(p,"--plot_width",
                  help="Plot width", default=10)
p <- add_argument(p, "--plot_height",
                  help="Plot height. If NA, scales to the number of haplotypes.", default=NA)
args <- parse_args(p)


if (is.na(args$chr)) {
  stop("Chromosome (--chr) must be provided.")
}

# Create default output names.
args$output_tbl <- ifelse(
  !is.null(args$output_tbl),
  paste0(args$chr, "_length.tsv"),
  args$output_tbl
)
args$output <- ifelse(
  !is.null(args$output),
  paste0(args$chr, "_hgsvc3_",
         args$order_hor,
         switch(args$order_hor, "small"="er", "large"="r"),
         "ontop.png"),
  args$output
)

cols_to_take <- seq(5)
monomer_len <- 170
cols =  c("chr", "start", "stop", "hor", "strand")

chm13 <- fread(args$input_chm13,
               sep="\t",
               stringsAsFactors=TRUE,
               fill=TRUE, quote="",
               header=FALSE, select=cols_to_take)
chm1 <- fread(args$input_chm1,
              sep="\t",
              stringsAsFactors=TRUE,
              fill=TRUE, quote="",
              header=FALSE, select=cols_to_take)
hgsvc3 <- fread(args$input,
                sep="\t",
                stringsAsFactors=TRUE,
                fill=TRUE, quote="",
                header=FALSE, select=cols_to_take)

colnames(chm13) <- cols
colnames(chm1) <- cols
colnames(hgsvc3) <- cols

#select the chr
chm13_select <- chm13 %>%
  filter(str_detect(chr, paste0(args$chr, "$")))
chm1_select <- chm1 %>%
  filter(str_detect(chr, paste0(args$chr, "$")))

#combine the CHM1 and HGSVC3 centromeres
chm13_select$chr <- gsub("chr", "chm13_chr", chm13_select$chr)
chm1_select$chr <- gsub("chr", "chm1_chr", chm1_select$chr)
chm13_chm1_hgsvc3 <- rbind(chm13_select, chm1_select, hgsvc3)

#determine distance between start and stop
chm13_chm1_hgsvc3$length <- chm13_chm1_hgsvc3$stop-chm13_chm1_hgsvc3$start

#calculate monomer size and round
chm13_chm1_hgsvc3$mer <- as.numeric(round(chm13_chm1_hgsvc3$length/monomer_len))

#filter monomers
chm13_chm1_hgsvc3 <- switch (
  args$chr,
  "chr10" = subset(chm13_chm1_hgsvc3, as.numeric(mer) >= 5),
  "chr20" = subset(chm13_chm1_hgsvc3, as.numeric(mer) >= 5),
  "chrY" = subset(chm13_chm1_hgsvc3, as.numeric(mer) >= 30),
  "chr17" = subset(chm13_chm1_hgsvc3, as.numeric(mer) >= 4),
  chm13_chm1_hgsvc3)

#determine length of HOR array
df_length <- chm13_chm1_hgsvc3 %>% group_by(chr)

if (args$chr == "chr8") {
  df_length <- df_length %>% slice(3:n()) %>% filter(row_number() <= n()-1)
}

df_length <- df_length %>%
  summarise(HOR_array_length = max(stop) - min(start),
            start=min(start),
            stop=max(stop))

#write bed output
write.table(df_length, file=args$output_tbl, quote=FALSE, sep='\t', row.names = FALSE)

# filter for HORs that occur at least 20 times (10 times per haplotype)
df_mer <- chm13_chm1_hgsvc3 %>% group_by(mer) %>% filter(n()>args$hor_filter) #

#change the start coordinate of each
df_final = df_mer %>%
  group_by(chr) %>%
  mutate(stop = stop - min(start)) %>%
  mutate(start = start - min(start))

# plot (new colors to emphasize novel HORs)
myColors <- c("#A8275C", "#9AC78A", "#A53D63", "#3997C6", "#29A3CE", "#5EB2A7", "#38A49B", "#45B4CE", "#A53D63", "#AA1B63", "#3F66A0",
              "#D66C54", "#BFDD97", "#C0D875", "#E5E57A", "#B75361", "#F9E193","#C6625D", "#E5D1A1", "#A1B5E5", "#9F68A5", "#81B25B",
              "#F4DC78", "#7EC0B3", "#A8275C", "#8CC49F", "#893F89", "#6565AA")
names(myColors) <- levels(as.factor(c("1", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19",
                                      "2", "20", "21", "22", "24", "26", "3", "30", "32", "34", "35",
                                      "4", "5", "6", "7", "8", "9")))

# reorder data to put larger or smaller -mers on top
df_final <- switch(
  args$order_hor,
  # First sort by val. This sorts the dataframe but NOT the factor levels #larger HORs on top
  "large" = df_final %>% arrange(mer),
  # First sort by val. This sorts the dataframe but NOT the factor levels #smaller HORs on top
  "small" =  df_final %>% arrange(-mer),
  stop(paste("Invalid mer reordering option:", args$order_hor)))

# Just the summary track
ggplot(df_final) +
  #geom_segment(aes(x = start, xend = stop+1000, y = mer, yend = mer, color = mer), size=7) +
  geom_segment(aes(x = start, xend = stop+2000, y = chr, yend = chr, color = as.factor(mer)), linewidth = 9) +
  #geom_segment(aes(x = 0, xend = 1000000, y = "scale", yend = "scale", color = "black"), size = 9) +
  scale_color_manual(values=myColors) +
  theme_classic() +
  xlab("Position") + ylab("HOR size") +
  guides(color = guide_legend(override.aes = list(size = 7)))
#facet_wrap(~chr, ncol=2, scales="free_y")

# save plot
if (is.na(args$plot_height)) {
  # 0.4 inches per chr
  args$plot_height <- round(length(levels(df_final$chr)) * 0.4)
}

# Prevent ggplot from creating Rplots.pdf when calling ggsave.
if(!interactive()) pdf(NULL)
ggsave(args$output, width=args$plot_width, height=args$plot_height, limitsize = FALSE)
