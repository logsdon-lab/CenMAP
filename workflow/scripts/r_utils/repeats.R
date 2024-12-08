get_humas_hmmer_stv_annot_colors <- function() {
  myColors <- c(
    "#A8275C", "#9AC78A", "#CC8FC1", "#3997C6", "#8882C4", "#8ABDD6", "#096858", "#45B4CE", "#AFA7D8", "#A874B5", "#3F66A0",
    "#D66C54", "#BFDD97", "#AF5D87", "#E5E57A", "#ED975D", "#F9E193", "#93430C", "#E5D1A1", "#A1B5E5", "#9F68A5", "#81B25B",
    "#F4DC78", "#7EC0B3", "#B23F73", "#8CC49F", "#893F89", "#6565AA"
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

add_cdr_to_plot <- function(plot, df_cdr) {
  if (!(typeof(df_cdr) == "logical"))  {
    plot <- plot +
      geom_segment(
        data = df_cdr,
        aes(x = start2, y = chr, xend = stop2, yend = chr),
        linewidth = 1,
        position = position_nudge(y=0.3),
        colour = "black",
      )
  }
  return(plot)
}

add_stv_ort_to_plot <- function(plot, df_stv_ort) {
  if (!(typeof(df_stv_ort) == "logical")) {
    plot <- plot +
      geom_segment(
        data=df_stv_ort,
        aes(
          x = start2,
          xend = stop2,
          y = chr,
          yend = chr,
        ),
        position = position_nudge(y=-0.45),
        size = 1.5,
        lineend = "butt",
        linejoin = "mitre",
        # Point arrow to last position. See https://rdrr.io/r/grid/arrow.html
        arrow = arrow(length = unit(0.2, "cm"), ends = "last", type = "closed")
      )
  }
  return(plot)
}

plot_single_ctg <- function(ctg, df_rm_sat_out, df_humas_hmmer_stv_out, df_cdr, df_stv_ort, height = 10) {
  plot <- ggplot(data = df_rm_sat_out[order(df_rm_sat_out$region), ] %>% filter(chr == ctg)) +
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

  if (!(typeof(df_cdr) == "logical"))  {
    df_cdr_ctg <- df_cdr %>% filter(chr == ctg)
  } else {
    df_cdr_ctg <- df_cdr
  }
  if (!(typeof(df_stv_ort) == "logical"))  {
    df_stv_ort_ctg <- df_stv_ort %>% filter(chr == ctg)
  } else {
    df_stv_ort_ctg <- df_stv_ort
  }
  plot <- add_cdr_to_plot(plot, df_cdr_ctg)
  plot <- add_stv_ort_to_plot(plot, df_stv_ort_ctg)
  return(plot)
}

plot_all_ctgs <- function(df_rm_sat_out, df_humas_hmmer_stv_out, df_cdr, df_stv_ort, height = 8) {
  plot <- ggplot(data = df_rm_sat_out[order(df_rm_sat_out$region), ]) +
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
      axis.line.y = element_line(colour = "white"),
      axis.title.y = element_blank()
    ) +
    scale_x_continuous(
      labels = scales::unit_format(scale = 1e-6, accuracy=0.1, unit="")
    ) +
    guides(color = guide_legend(override.aes = list(size = 6))) +
    xlab("Position (Mbp)")

  plot <- add_cdr_to_plot(plot, df_cdr)
  plot <- add_stv_ort_to_plot(plot, df_stv_ort)
  return(plot)
}

read_multiple_repeatmasker_sat_input <- function(input_file) {
  cols <- c("chr", "start", "stop", "region", "value", "strand", "start2", "stop2", "rgb")
  df <- fread(
    input_file,
    select = c(1:9),
    stringsAsFactors = TRUE,
    fill = TRUE,
    sep = "\t",
    quote = "",
    header = FALSE,
    col.names = cols
  )

  if (nrow(df) == 0) {
    df <- data.table(matrix(ncol=length(cols), nrow=0))
    colnames(df) <- cols
  }
  # Filter duplicated chm13 rows.
  df <- df %>%
    mutate(chr=str_replace(chr, "cen", "chr")) %>%
    filter(!str_detect(chr, "^chr[0-9XY]+$")) %>%
    # Correct for version and different naming of chr. ex. chm1_cen1v8 -> chm1_chr1
    mutate(chr=str_remove(chr, "v\\d+")) %>%
    mutate(
      ctg_name = str_extract(chr, "^(.*?):|^(.*?)$", 1)
    ) %>%
    select(-chr) %>%
    mutate(chr=ctg_name)

  # reorder rows so that live arrays are plotted on top
  df$region <- factor(df$region, levels = c("ct", "asat", "bsat", "gsat", "hsat1A", "hsat1B", "hsat2", "hsat3"), ordered = T)

  return(df)
}

read_one_repeatmasker_sat_input <- function(input_file) {
  # read in BED file
  cols <- c("chr", "start", "stop", "region", "value", "strand", "start2", "stop2", "rgb")
  df <- fread(
    input_file,
    select = c(1:9),
    stringsAsFactors = TRUE,
    fill = TRUE,
    sep = "\t",
    quote = "",
    header = FALSE,
    col.names = cols
  )
  if (nrow(df) == 0) {
    df <- data.table(matrix(ncol=length(cols), nrow=0))
    colnames(df) <- cols
  }
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
      ctg_start = as.integer(replace_na(str_extract(chr, ":(\\d+)-", 1), "0")),
      ctg_stop = as.integer(replace_na(str_extract(chr, "-(\\d+)$", 1), "0"))
    ) %>%
    # Because coords originate from rm output, already correctly oriented.
    mutate(
      start2 = start - ctg_start,
      stop2 = stop - ctg_start
    ) %>%
    mutate(
      ctg_name = str_extract(chr, "^(.*?):|^(.*?)$", 1)
    ) %>%
    select(-chr) %>%
    mutate(chr=ctg_name)

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
    header = FALSE, select = cols_to_take, col.names = cols
  )
  if (nrow(df_samples) == 0) {
    df_samples <- data.table(matrix(ncol=length(cols), nrow=0))
    colnames(df_samples) <- cols
  }

  chr_name <- first(str_extract(df_samples$chr, "(chr[\\dXY]+)"))
  if (length(chr_name) == 0) {
    chr_name = c("")
  }
  # determine distance between start and stop
  df_samples$length <- df_samples$stop - df_samples$start

  # calculate monomer size and round
  df_samples$mer <- as.numeric(round(df_samples$length / monomer_len))

  # filter monomers
  df_samples <- switch(chr_name,
    "chr10" = subset(df_samples, as.numeric(mer) >= 5),
    "chr20" = subset(df_samples, as.numeric(mer) >= 5),
    "chrY" = subset(df_samples, as.numeric(mer) >= 20),
    "chr17" = subset(df_samples, as.numeric(mer) >= 4),
    df_samples
  )

  # filter for HORs that occur at least 20 times (10 times per haplotype)
  df_stv <- df_samples %>%
    group_by(mer) %>%
    filter(n() > hor_filter)

  df_stv <- df_stv %>%
    mutate(
      ctg_start = as.integer(replace_na(str_extract(chr, ":(\\d+)-", 1), "0")),
      ctg_stop = as.integer(replace_na(str_extract(chr, "-(\\d+)$", 1), "0"))
    ) %>%
    mutate(
      start2 = start - ctg_start,
      stop2 = stop - ctg_start
    ) %>%
    mutate(start=start2, stop=stop2) %>%
    mutate(
      ctg_name = str_extract(chr, "^(.*?):|^(.*?)$", 1)
    ) %>%
    select(-chr) %>%
    mutate(chr=ctg_name)

   return(df_stv)
}

read_multiple_humas_hmmer_input <- function(
    input_chr,
    chr_name,
    hor_filter = 0
) {
  cols_to_take <- seq(5)
  monomer_len <- 170
  cols <- c("chr", "start", "stop", "hor", "strand")

  samples <- fread(input_chr,
    sep = "\t",
    stringsAsFactors = TRUE,
    fill = TRUE, quote = "",
    header = FALSE, select = cols_to_take, col.names = cols
  )
  if (nrow(samples) == 0) {
    samples <- data.table(matrix(ncol=length(cols), nrow=0))
    colnames(samples) <- cols
  }

  # determine distance between start and stop
  samples$length <- samples$stop - samples$start

  # calculate monomer size and round
  samples$mer <- as.numeric(round(samples$length / monomer_len))

  # filter monomers
  samples <- switch(chr_name,
    "chr10" = subset(samples, as.numeric(mer) >= 5),
    "chr20" = subset(samples, as.numeric(mer) >= 5),
    "chrY" = subset(samples, as.numeric(mer) >= 30),
    "chr17" = subset(samples, as.numeric(mer) >= 4),
    samples
  )

  # filter for HORs that occur at least 20 times (10 times per haplotype)
  df_stv <- samples %>%
    group_by(mer) %>%
    filter(n() > hor_filter) %>%
    mutate(
      ctg_name = str_extract(chr, "^(.*?):|^(.*?)$", 1)
    ) %>%
    select(-chr) %>%
    mutate(chr=ctg_name)

  return(df_stv)
}

read_one_cdr_input <- function(input_cdr) {
  if (is.na(input_cdr)) {
    return(NA)
  }
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
      ctg_start = as.integer(replace_na(str_extract(chr, ":(\\d+)-", 1), "0")),
      ctg_stop = as.integer(replace_na(str_extract(chr, "-(\\d+)$", 1), "0"))
    ) %>%
    mutate(
      start2 = start - ctg_start,
      stop2 = stop - ctg_start
    ) %>%
    mutate(
      ctg_name = str_extract(chr, "^(.*?):|^(.*?)$", 1)
    ) %>%
    select(-chr) %>%
    mutate(chr=ctg_name)

  return(df_cdr)
}

read_one_methyl_bed_input <- function(input_methyl) {
   if (is.na(input_methyl)) {
    return(NA)
  }
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
      ctg_start = as.integer(replace_na(str_extract(chr, ":(\\d+)-", 1), "0")),
      ctg_stop = as.integer(replace_na(str_extract(chr, "-(\\d+)$", 1), "0"))
    ) %>%
    mutate(
      start2 = start - ctg_start,
      stop2 = stop - ctg_start
    ) %>%
    select(chr, start2, stop2, meth_prob) %>%
    mutate(
      ctg_name = str_extract(chr, "^(.*?):|^(.*?)$", 1)
    ) %>%
    select(-chr) %>%
    mutate(chr=ctg_name)

  return(df_methyl_binned)
}

read_one_hor_mon_ort_input <- function(input_hor_ort) {
   if (is.na(input_hor_ort)) {
    return(NA)
  }
  cols <- c("chr", "start", "stop", "strand")
  df_hor_ort <- fread(
    input_hor_ort,
    header = FALSE,
    select = c(1:4),
    col.names = cols
  )
  if (nrow(df_hor_ort) == 0) {
    df_hor_ort <- data.table(matrix(ncol=length(cols), nrow=0))
    colnames(df_hor_ort) <- cols
  }
  df_hor_ort <- df_hor_ort %>%
    mutate(
      ctg_start = as.integer(replace_na(str_extract(chr, ":(\\d+)-", 1), "0")),
      ctg_stop = as.integer(replace_na(str_extract(chr, "-(\\d+)$", 1), "0"))
    ) %>%
    mutate(
      # Reverse arrows if -
      # Add some bases to not visually clutter arrows.
      # Not true but this image is only for qualitative analysis.
      # TODO: This may fail on smaller contigs. Maybe group_by to find length and use some prop (1% maybe.)
      start2 = ifelse(
        strand == "+",
        start - ctg_start + 20000,
        stop - ctg_start - 20000
      ),
      stop2 = ifelse(
        strand == "+",
        stop - ctg_start - 30000,
        start - ctg_start + 30000
      )
    ) %>%
    select(chr, start2, stop2, strand) %>%
    mutate(
      ctg_name = str_extract(chr, "^(.*?):|^(.*?)$", 1)
    ) %>%
    select(-chr) %>%
    mutate(chr=ctg_name)

  return(df_hor_ort)
}
