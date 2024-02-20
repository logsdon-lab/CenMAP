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


p <- arg_parser("Plot STV")
p <- add_argument(p, "--input", help="Input HGSVC3 file", type="character")
p <- add_argument(p, "--input_chm13", help="Input CHM13 HumAS-HMMER formatted output", type="character")
p <- add_argument(p, "--input_chm1", help="Input CHM1 HumAS-HMMER formatted output", type="character")
p <- add_argument(p, "--output", help="Output plot file.", type="character")
argv <- parse_args(p)


rm(list=ls(all=TRUE)) # clear variables
# setwd("~/Documents/Eichler Lab/Weekly plans/Ape_cen_assembly/CHM1/HORorganization/StV plots/CHM13_v1.1") # set to folder of your files
chm13 = fread(file.choose(), sep="\t", stringsAsFactors = TRUE, fill=TRUE, quote="", header=FALSE, select=c(1,2,3,4,5))
# setwd("~/Documents/Eichler Lab/Weekly plans/Ape_cen_assembly/CHM1/HORorganization/StV plots/v21") # set to folder of your files
chm1 = fread(file.choose(), sep="\t", stringsAsFactors = TRUE, fill=TRUE, quote="", header=FALSE, select=c(1,2,3,4,5))
# setwd("~/Documents/Eichler Lab/Weekly plans/HGSVC3/HORorganization/") # set to folder of your files
hgsvc3 = fread(file.choose(), sep="\t", stringsAsFactors = TRUE, fill=TRUE, quote="", header=FALSE, select=c(1,2,3,4,5))
cols =  c("chr", "start", "stop", "hor", "strand")
colnames(chm13) <- cols
colnames(chm1) <- cols
colnames(hgsvc3) <- cols

#select the chr
chm13_select <- chm13 %>%
  filter(str_detect(chr, "chr21$"))
chm1_select <- chm1 %>%
  filter(str_detect(chr, "chr21$"))

#combine the CHM1 and HGSVC3 centromeres
chm13_select$chr <- gsub("chr", "chm13_chr", chm13_select$chr)
chm1_select$chr <- gsub("chr", "chm1_chr", chm1_select$chr)
chm13_chm1_hgsvc3 <- rbind(chm13_select, chm1_select, hgsvc3)

#determine distance between start and stop
chm13_chm1_hgsvc3$length <- chm13_chm1_hgsvc3$stop-chm13_chm1_hgsvc3$start

#calculate monomer size and round
chm13_chm1_hgsvc3$mer <- as.numeric(round(chm13_chm1_hgsvc3$length/170))

#filter monomers
#chm13_chm1_hgsvc3_subset = subset(chm13_chm1_hgsvc3, as.numeric(mer) >= 5) #(for chr10)
#chm13_chm1_hgsvc3 = subset(chm13_chm1_hgsvc3, as.numeric(mer) >= 5) #(for chr20)
#chm13_chm1_hgsvc3 = subset(chm13_chm1_hgsvc3, as.numeric(mer) >= 30) #(for chrY)
#chm13_chm1_hgsvc3 = subset(chm13_chm1_hgsvc3, as.numeric(mer) >= 4) #(for chr17)

#determine length of HOR array
df_length = chm13_chm1_hgsvc3 %>%
  group_by(chr) %>%
  #slice(3:n()) %>% #for chr8
  #filter(row_number() <= n()-1) %>% #for chr8
  summarise(HOR_array_length = max(stop) - min(start))

#write bed output
write.table(df_length, file='chr14_length.tsv', quote=FALSE, sep='\t', row.names = FALSE)

# filter for HORs that occur at least 20 times (10 times per haplotype)
min_HORs = 0
df_mer <-  chm13_chm1_hgsvc3 %>% group_by(mer) %>% filter(n()>min_HORs) #

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
df_final <- df_final %>%
  arrange(mer) # First sort by val. This sorts the dataframe but NOT the factor levels #larger HORs on top
df_final <- df_final %>%
  arrange(-mer) # First sort by val. This sorts the dataframe but NOT the factor levels #smaller HORs on top

# Just the summary track
ggplot(df_final) +
  #geom_segment(aes(x = start, xend = stop+1000, y = mer, yend = mer, color = mer), size=7) +
  geom_segment(aes(x = start, xend = stop+2000, y = chr, yend = chr, color = as.factor(mer)), size = 9) +
  #geom_segment(aes(x = 0, xend = 1000000, y = "scale", yend = "scale", color = "black"), size = 9) +
  scale_color_manual(values=myColors) +
  theme_classic() +
  xlab("Position") + ylab("HOR size") +
  guides(color = guide_legend(override.aes = list(size = 7)))
#facet_wrap(~chr, ncol=2, scales="free_y")

# save plot
ggsave(("cens_to_annotate2_smallerontop.png"), width=10, height=11)
ggsave(("cens_to_annotate2_largerontop.png"), width=10, height=11)
ggsave(("chr14_hgsvc3_smallerontop.png"), width=10, height=5)
ggsave(("chr14_hgsvc3_largerontop.png"), width=10, height=5)
ggsave(("chrX_hgsvc3_smallerontop.png"), width=10, height=13)
ggsave(("chrX_hgsvc3_largerontop.png"), width=10, height=13)
ggsave(("chrY_hgsvc3_smallerontop.png"), width=10, height=3)
ggsave(("chrY_hgsvc3_largerontop.png"), width=10, height=3)
ggsave(("chr22_hgsvc3_smallerontop.png"), width=10, height=11)
ggsave(("chr22_hgsvc3_largerontop.png"), width=10, height=11)
ggsave(("chr21_hgsvc3_smallerontop.png"), width=10, height=10)
ggsave(("chr21_hgsvc3_largerontop.png"), width=10, height=10)
ggsave(("chr20_hgsvc3_smallerontop.png"), width=10, height=12)
ggsave(("chr20_hgsvc3_largerontop.png"), width=10, height=12)
ggsave(("cens_to_annotate_smallerontop.png"), width=10, height=16)
ggsave(("cens_to_annotate_largerontop.png"), width=10, height=16)
ggsave(("chr19_hgsvc3_smallerontop.png"), width=10, height=15)
ggsave(("chr19_hgsvc3_largerontop.png"), width=10, height=15)
ggsave(("chr18_hgsvc3_smallerontop.png"), width=10, height=8)
ggsave(("chr18_hgsvc3_largerontop.png"), width=10, height=8)
ggsave(("chr17_hgsvc3_smallerontop.png"), width=10, height=15)
ggsave(("chr17_hgsvc3_largerontop.png"), width=10, height=15)
ggsave(("chr16_hgsvc3_smallerontop.png"), width=10, height=14)
ggsave(("chr16_hgsvc3_largerontop.png"), width=10, height=14)
ggsave(("chr15_hgsvc3_smallerontop.png"), width=10, height=17)
ggsave(("chr15_hgsvc3_largerontop.png"), width=10, height=17)
ggsave(("chr13_hgsvc3_smallerontop.png"), width=10, height=20)
ggsave(("chr13_hgsvc3_largerontop.png"), width=10, height=20)
ggsave(("chr12_hgsvc3_smallerontop.png"), width=10, height=17)
ggsave(("chr12_hgsvc3_largerontop.png"), width=10, height=17)
ggsave(("chr11_hgsvc3_smallerontop.png"), width=10, height=15)
ggsave(("chr11_hgsvc3_largerontop.png"), width=10, height=15)
ggsave(("chr10_hgsvc3_smallerontop.png"), width=10, height=6)
ggsave(("chr10_hgsvc3_largerontop.png"), width=10, height=6)
ggsave(("chr9_hgsvc3_smallerontop.png"), width=10, height=17)
ggsave(("chr9_hgsvc3_largerontop.png"), width=10, height=17)
ggsave(("chr8_hgsvc3_smallerontop.png"), width=10, height=16)
ggsave(("chr8_hgsvc3_largerontop.png"), width=10, height=16)
ggsave(("chr7_hgsvc3_smallerontop.png"), width=10, height=12)
ggsave(("chr7_hgsvc3_largerontop.png"), width=10, height=12)
ggsave(("chr6_hgsvc3_smallerontop.png"), width=10, height=15)
ggsave(("chr6_hgsvc3_largerontop.png"), width=10, height=15)
ggsave(("chr5_hgsvc3_smallerontop.png"), width=10, height=15)
ggsave(("chr5_hgsvc3_largerontop.png"), width=10, height=15)
ggsave(("chr4_hgsvc3_smallerontop.png"), width=10, height=11)
ggsave(("chr4_hgsvc3_largerontop.png"), width=10, height=11)
ggsave(("chr3_hgsvc3_smallerontop.png"), width=10, height=8)
ggsave(("chr3_hgsvc3_largerontop.png"), width=10, height=8)
ggsave(("chr2_hgsvc3_smallerontop.png"), width=10, height=17)
ggsave(("chr2_hgsvc3_largerontop.png"), width=10, height=17)
ggsave(("chr1_hgsvc3_smallerontop.png"), width=10, height=4)
ggsave(("chr1_hgsvc3_largerontop.png"), width=10, height=4)
ggsave(("misassembledcontigs_hgsvc3_smallerontop.png"), width=10, height=3)
ggsave(("misassembledcontigs_hgsvc3_largerontop.png"), width=10, height=3)
