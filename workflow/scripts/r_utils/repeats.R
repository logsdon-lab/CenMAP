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
    "#DDDDDD",
    "#3A3A3A",
    "#3A3A3A",
    "#3A3A3A",
    "#3A3A3A",
    "#3A3A3A"
  )
  names(myColors) <- levels(as.factor(c(
    "asat",
    "bsat",
    "ct",
    "gsat",
    "hsat1A",
    "hsat1B",
    "hsat2",
    "hsat3"
  )))
  return(myColors)
}

plot_single_ctg <- function(ctg, df_rm_sat_out, df_humas_hmmer_stv_out, height = 10) {
  p <- ggplot(data = df_rm_sat_out[order(df_rm_sat_out$region), ] %>% filter(chr == ctg)) +
    geom_segment(
      aes(x = start2, y = chr, xend = stop2 + 1000, yend = chr, color = region),
      alpha = 1,
      size = height
    ) +
    scale_color_manual(values = get_rm_sat_annot_colors()) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 12)) +
    # New colorscale.
    new_scale_color() +
    geom_segment(
      data = df_humas_hmmer_stv_out %>% filter(chr == ctg),
      aes(x = start, xend = stop + 2000, y = chr, yend = chr, color = as.factor(mer)),
      size = height
    ) +
    scale_color_manual(values = get_humas_hmmer_stv_annot_colors()) +
    theme_classic() +
    # Remove everything.
    theme(
      axis.line.x = element_line(colour = "white"),
      axis.line.y = element_line(colour = "white"),
      axis.title.y = element_blank(),
      axis.ticks.y = element_blank(),
      legend.position = "none"
    ) +
    guides(color = guide_legend(override.aes = list(size = 6))) +
    scale_x_continuous(
      labels = scales::unit_format(scale = 1e-6, accuracy=0.1, unit="")
    ) +
    xlab("Position (Mbp)")
  return(p)
}

plot_all_ctgs <- function(df_rm_sat_out, df_humas_hmmer_stv_out, height = 10) {
  p <- ggplot(data = df_rm_sat_out[order(df_rm_sat_out$region), ]) +
    geom_segment(
      aes(x = start2, y = chr, xend = stop2 + 1000, yend = chr, color = region),
      alpha = 1,
      size = height
    ) +
    scale_color_manual(values = get_rm_sat_annot_colors()) +
    # New colorscale.
    new_scale_color() +
    geom_segment(
      data = df_humas_hmmer_stv_out,
      aes(x = start, xend = stop + 2000, y = chr, yend = chr, color = as.factor(mer)),
      size = height
    ) +
    scale_color_manual(values = get_humas_hmmer_stv_annot_colors()) +
    theme_classic() +
    theme(
      legend.position = "bottom",
      legend.key.width = unit(1, "cm"),
      legend.title = element_blank(),
      axis.line.y = element_line(colour = "white")
    ) +
    scale_x_continuous(
      labels = scales::unit_format(scale = 1e-6, accuracy=0.1, unit="")
    ) +
    guides(color = guide_legend(override.aes = list(size = 6))) +
    xlab("Position (Mbp)")

  return(p)
}

read_multiple_repeatmasker_sat_input <- function(input_file) {
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

read_one_repeatmasker_sat_input <- function(input_file) {
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
    mutate(chr = str_replace(chr, "cen", "chr")) %>%
    filter(!str_detect(chr, "^chr[0-9XY]+$")) %>%
    # Correct for version and different naming of chr. ex. chm1_cen1v8 -> chm1_chr1
    mutate(chr = str_remove(chr, "v\\d+"))

  # reorder rows so that live arrays are plotted on top
  df$region <- factor(df$region, levels = c("ct", "asat", "bsat", "gsat", "hsat1A", "hsat1B", "hsat2", "hsat3"), ordered = T)

  # Adjust for reverse complemented regions.
  # Set to always start at 0.
  df <- df %>%
    mutate(
      ctg_start = as.integer(replace_na(str_extract(chr, ":(\\d+)-", 1), 0)),
      ctg_stop = as.integer(replace_na(str_extract(chr, "-(\\d+)$", 1), 0))
    ) %>%
    # Because coords originate from rm output, already correctly oriented.
    mutate(
      start2 = start - ctg_start,
      stop2 = stop - ctg_start
    )

  return(df)
}

read_one_humas_hmmer_input <- function(
    input_chr,
    hor_filter = 0) {
  cols_to_take <- seq(6)
  monomer_len <- 170
  cols <- c("chr", "start", "stop", "hor", "0", "strand")

  df_samples <- fread(input_chr,
    sep = "\t",
    stringsAsFactors = TRUE,
    fill = TRUE, quote = "",
    header = FALSE, select = cols_to_take
  )
  colnames(df_samples) <- cols

  chr_name <- str_extract(df_samples$chr, "(chr[\\dXY]+)")[[1]]

  # determine distance between start and stop
  df_samples$length <- df_samples$stop - df_samples$start

  # calculate monomer size and round
  df_samples$mer <- as.numeric(round(df_samples$length / monomer_len))

  # filter monomers
  df_samples <- switch(chr_name,
    "chr10" = subset(df_samples, as.numeric(mer) >= 5),
    "chr20" = subset(df_samples, as.numeric(mer) >= 5),
    "chrY" = subset(df_samples, as.numeric(mer) >= 30),
    "chr17" = subset(df_samples, as.numeric(mer) >= 4),
    df_samples
  )

  # filter for HORs that occur at least 20 times (10 times per haplotype)
  df_stv <- df_samples %>%
    group_by(mer) %>%
    filter(n() > hor_filter)

  df_stv <- df_stv %>%
    mutate(
      ctg_start = as.integer(replace_na(str_extract(chr, ":(\\d+)-", 1), 0)),
      ctg_stop = as.integer(replace_na(str_extract(chr, "-(\\d+)$", 1), 0))
    ) %>%
    mutate(
      start2 = start - ctg_start,
      stop2 = stop - ctg_start
    ) %>%
    mutate(start=start2, stop=stop2)

   return(df_stv)
}

read_multiple_humas_hmmer_input <- function(
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

  return(df_stv)
}

read_one_cdr_input <- function(input_cdr) {
  df_cdr <- fread(
    input_cdr,
    header = FALSE,
    select = c(1:3),
    col.names = c("chr", "start", "stop")
  )
  if (nrow(df_cdr) == 0) {
    return(NA)
  }

  df_cdr <- df_cdr %>%
    mutate(
      ctg_start = as.integer(replace_na(str_extract(chr, ":(\\d+)-", 1), 0)),
      ctg_stop = as.integer(replace_na(str_extract(chr, "-(\\d+)$", 1), 0))
    ) %>%
    mutate(
      start2 = start - ctg_start,
      stop2 = stop - ctg_start
    )

  return(df_cdr)
}

read_one_methyl_bed_input <- function(input_methyl) {
  df_methyl_binned <- fread(
    input_methyl,
    header = FALSE,
    select = c(1:4),
    col.names = c("chr", "start", "stop", "meth_prob")
  )
  if (nrow(df_methyl_binned) == 0) {
    return(NA)
  }
  df_methyl_binned <- df_methyl_binned %>%
    mutate(
      ctg_start = as.integer(replace_na(str_extract(chr, ":(\\d+)-", 1), 0)),
      ctg_stop = as.integer(replace_na(str_extract(chr, "-(\\d+)$", 1), 0))
    ) %>%
    mutate(
      start2 = start - ctg_start,
      stop2 = stop - ctg_start
    ) %>%
    select(chr, start2, stop2, meth_prob)

  return(df_methyl_binned)
}
