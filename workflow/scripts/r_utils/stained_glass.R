make_scale <- function(vals) {
  comma(vals / 1e6)
}

make_k <- function(vals) {
  comma(vals / 1e3)
}

get_colors <- function(sdf) {
  bot <- floor(min(sdf$perID_by_events))
  top <- 100
  breaks <- seq(70, 100, 1)
  labels <- seq(length(breaks) - 1)
  return(cut(sdf$perID_by_events, breaks = seq(70, 100, 1), labels = seq(length(breaks) - 1), include.lowest = TRUE))
}

read_bedpe <- function(all.files) {
  l <- lapply(all.files, fread, sep = "\t")
  df <- rbindlist(l)

  df$discrete <- get_colors(df)
  if ("#query_name" %in% colnames(df)) {
    df$q <- df$`#query_name`
    df$q_st <- df$query_start
    df$q_en <- df$query_end
    df$r <- df$reference_name
    df$r_st <- df$reference_start
    df$r_en <- df$reference_end
  }
  window <- max(df$q_en - df$q_st)
  df$first_pos <- df$q_st / window
  df$second_pos <- df$r_st / window

  df <- df %>% mutate(
    new_discrete = case_when(
      perID_by_events > 0 & perID_by_events < 80 ~ "1-10",
      perID_by_events >= 80 & perID_by_events < 90 ~ "11-20",
      .default = as.character(discrete)
    )
  )
  return(df)
}

discrete_color_ranges <- c(
  "#4b3991", # 1-10
  "#2974af", # 11-20
  "#4a9da8", # 21
  "#57b894", # 22
  "#9dd893", # 23
  "#e1f686", # 24
  "#ffffb2", # 25
  "#fdda79", # 26
  "#fb9e4f", # 27
  "#ee5634", # 28
  "#c9273e", # 29
  "#8a0033" # 30
)

# get the lowest 0.1% of the data so we can not plot it
make_hist <- function(sdf) {
  bot <- quantile(sdf$perID_by_events, probs = 0.001)[[1]]
  p <- ggplot(sdf) +
    geom_histogram(aes(perID_by_events, fill = new_discrete), bins = 300) +
    scale_fill_manual(values = discrete_color_ranges) +
    theme_cowplot() +
    theme(
      legend.position = "none",
      text = element_text(size = 8),
      axis.text = element_text(size = 8)
    ) +
    scale_y_continuous(labels = make_k) +
    coord_cartesian(xlim = c(bot, 100)) +
    xlab("% identity") +
    ylab("# of alignments (k)")
  return(p)
}

diamond <- function(row) {
  side_length <- as.numeric(row["window"])
  x <- as.numeric(row["w"])
  y <- as.numeric(row["z"])
  base <- matrix(c(1, 0, 0, 1, -1, 0, 0, -1), nrow = 2) * sqrt(2) / 2
  trans <- (base * side_length) + c(x, y)
  df <- as.data.frame(t(trans))
  colnames(df) <- c("w", "z")
  df$new_discrete <- row["new_discrete"]
  df$group <- as.numeric(row["group"])
  return(df)
}

make_tri_df <- function(df) {
  df$w <- (df$first_pos + df$second_pos)
  df$z <- -df$first_pos + df$second_pos
  window <- max(df$q_en - df$q_st)
  tri_scale <- max(df$q_st) / max(df$w)
  df$window <- max(df$q_en - df$q_st) / tri_scale
  df$group <- seq(nrow(df))
  df_d <- rbindlist(apply(df, 1, diamond))
  df_d$x <- df_d$w * tri_scale
  df_d$y <- df_d$z * window
  return(df_d)
}

make_tri <- function(sdf, rname = "") {
  df_d <- make_tri_df(sdf)

  plt <- ggplot(df_d) +
    geom_polygon(aes(x = x, y = y, fill = new_discrete, group = group)) +
    theme_cowplot() +
    scale_fill_manual(values = discrete_color_ranges) +
    scale_x_continuous(labels = make_scale, limits = c(0, NA)) +
    scale_y_reverse(labels = make_scale, limits = c(NA, 0)) +
    xlab("Genomic position (Mbp)") +
    ylab("") +
    theme(
      legend.position = "none",
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.y = element_blank()
    ) +
    ggtitle(rname)
  return(plt)
}

make_dot <- function(sdf, rname = "") {
  max <- max(sdf$q_en, sdf$r_en)
  window <- max(sdf$query_end - sdf$query_start)
  plt <- ggplot(sdf) +
    geom_tile(aes(x = q_st, y = r_st, fill = new_discrete, height = window, width = window)) +
    theme_cowplot() +
    scale_fill_manual(values = discrete_color_ranges) +
    theme(legend.position = "none") +
    scale_x_continuous(labels = make_scale, limits = c(0, max)) +
    scale_y_continuous(labels = make_scale, limits = c(0, max)) +
    coord_fixed(ratio = 1) +
    facet_grid(r ~ q) +
    xlab("Genomic position (Mbp)") +
    ylab("") +
    ggtitle(rname)

  return(plt)
}


make_cen_plot <- function(rname, df_seq_ident, df_humas_hmmer_stv_out, df_rm_sat_out) {
  df_rname_seq_ident <- df_seq_ident %>% filter(q == rname & r == rname)

  # make the tri sequence identity plots
  df_d <- make_tri_df(df_rname_seq_ident)
  # make the histogram
  plot_hist <- make_hist(df_rname_seq_ident)
  grob_hist <- ggplotGrob(plot_hist)
  df_rm_sat_out <- df_rm_sat_out[order(df_rm_sat_out$region)] %>%
    filter(chr == rname)
  df_rm_sat_out_ct <- df_rm_sat_out %>% filter(region == "ct")
  df_rm_sat_out_other <- df_rm_sat_out %>% filter(region != "ct")

  segment_linewidth <- 10
  segment_y <- -200000
  ct_outline_edges_x <- 9000

  plot_ident_cen <- ggplot() +
    # Make larger ct segment as outline
    geom_segment(
      data = df_rm_sat_out_ct,
      aes(
        x = start2,
        y = segment_y,
        xend = stop2,
        yend = segment_y,
        color = "black"
      ),
      linewidth = segment_linewidth + 1
    ) +
    geom_segment(
      data = df_rm_sat_out_ct,
      aes(
        x = start2 + ct_outline_edges_x,
        y = segment_y,
        xend = stop2 - ct_outline_edges_x,
        yend = segment_y,
        color = region
      ),
      linewidth = segment_linewidth
    ) +
    geom_segment(
      data = df_rm_sat_out_other,
      aes(x = start2, y = segment_y, xend = stop2, yend = segment_y, color = region),
      linewidth = segment_linewidth
    ) +
    scale_color_manual(values = get_rm_sat_annot_colors()) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 12)) +
    # New colorscale.
    new_scale_color() +
    geom_segment(
      data = df_humas_hmmer_stv_out %>% filter(chr == rname),
      aes(
        x = start,
        xend = stop,
        y = segment_y,
        yend = segment_y,
        color = as.factor(mer)
      ),
      linewidth = segment_linewidth
    ) +
    scale_color_manual(values = get_humas_hmmer_stv_annot_colors()) +
    # New colorscale.
    new_scale_color() +
    geom_polygon(
      df_d,
      mapping = aes(x = x, y = y, fill = new_discrete, group = group),
      show.legend = FALSE
    ) +
    scale_fill_manual(values = discrete_color_ranges) +
    scale_x_continuous(labels = make_scale, limits = c(0, NA)) +
    # Reorient so stainedglass plot is on the bottom
    scale_y_reverse(labels = make_scale) +
    theme_classic() +
    theme(
      legend.position = "right",
      legend.title = element_blank(),
      legend.background = element_blank(),
      legend.key.size = unit(1.0, "cm"),
      legend.key = element_rect(color = "black", fill = NA),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.y = element_blank(),
    ) +
    xlab("Position (Mbp)")

  # ranges for inset hist
  xmax <- max(df_rname_seq_ident$q_en, df_rname_seq_ident$r_en)
  ylim <- max(df_d$y)
  final_plot <- plot_ident_cen +
    annotation_custom(
      grob_hist,
      xmin = 0, xmax = xmax / 3,
      ymin = -ylim, ymax = -ylim + (ylim / 3.5)
    ) +
    ggtitle(rname)

  return(final_plot)
}
