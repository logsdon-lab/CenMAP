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

plot_single_ctg <- function(ctg, df_rm_sat_out, df_humas_hmmer_stv_out, height = 10) {
  p <- ggplot(data = df_rm_sat_out[order(df_rm_sat_out$region), ] %>% filter(chr == ctg)) +
      geom_segment(
        aes(x = start2, y = chr, xend = stop2, yend = chr, color = region),
        alpha = 1,
        size = height
      ) +
      scale_color_manual(values = get_rm_sat_annot_colors()) +
      scale_x_continuous(breaks = scales::pretty_breaks(n = 12)) +
      # New colorscale.
      new_scale_color() +
      geom_segment(
        data = df_humas_hmmer_stv_out %>% filter(chr == ctg),
        aes(x = start, xend = stop , y = chr, yend = chr, color = as.factor(mer)),
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
  return(p)
}

plot_all_ctgs <- function(df_rm_sat_out, df_humas_hmmer_stv_out, height=10) {
  p <- ggplot(data = df_rm_sat_out[order(df_rm_sat_out$region), ]) +
    geom_segment(
      aes(x = start2, y = chr, xend = stop2, yend = chr, color = region),
      alpha = 1,
      size = height
    ) +
    scale_color_manual(values = get_rm_sat_annot_colors()) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 12)) +
    # New colorscale.
    new_scale_color() +
    geom_segment(
      data = df_humas_hmmer_stv_out,
      aes(x = start, xend = stop, y = chr, yend = chr, color = as.factor(mer)),
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

  return(p)
}

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
