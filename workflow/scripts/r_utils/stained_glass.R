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
  ctg_start <- min(df$query_start)
  ctg_end <- max(df$query_end)

  df$discrete <- get_colors(df)
  if ("#query_name" %in% colnames(df)) {
    df$q <- df$`#query_name`
    df <- df %>%
      mutate(
        q_st = query_start - ctg_start,
        r_st = reference_start - ctg_start,
        q_en = query_end - ctg_start,
        r_en = reference_end - ctg_start
      )
    df$r <- df$reference_name
  }
  window <- max(df$q_en - df$q_st)
  df$first_pos <- df$q_st / window
  df$second_pos <- df$r_st / window

  df <- df %>%
    mutate(
      new_discrete = case_when(
        perID_by_events > 0 & perID_by_events < 90 ~ "1-10",
        perID_by_events >= 90 & perID_by_events < 97.5 ~ "11-20",
        # TODO: There should be a way to do this automatically.
        #  100 - 95.0 = 5.0 / 10 = 0.5 increment
        #  100 - 97.5 = 2.5 / 10 = 0.25 increment
        perID_by_events >= 97.5 & perID_by_events < 97.75 ~ "21",
        perID_by_events >= 97.75 & perID_by_events < 98.0 ~ "22",
        perID_by_events >= 98.0 & perID_by_events < 98.25 ~ "23",
        perID_by_events >= 98.25 & perID_by_events < 98.5 ~ "24",
        perID_by_events >= 98.5 & perID_by_events < 98.75 ~ "25",
        perID_by_events >= 98.75 & perID_by_events < 99.0 ~ "26",
        perID_by_events >= 99.0 & perID_by_events < 99.25 ~ "27",
        perID_by_events >= 99.25 & perID_by_events < 99.5 ~ "28",
        perID_by_events >= 99.5 & perID_by_events < 99.75 ~ "29",
        perID_by_events >= 99.75 & perID_by_events < 100.0 ~ "30",
        .default = as.character(discrete)
      )
    ) %>%
    mutate(
      ctg_name = str_extract(r, "^(.*?):|^(.*?)$", 1),
      or = r
    ) %>%
    select(-c(r, q)) %>%
    mutate(r=ctg_name, q=ctg_name)

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
    xlab("Average Nucleotide Identity (%)") +
    ylab("Number of Intervals (k)")
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

make_cen_plot <- function(rname, df_seq_ident, df_humas_hmmer_stv_out, df_rm_sat_out, df_cdr, df_hor_ort, df_methyl_binned) {
  df_rname_seq_ident <- df_seq_ident %>% filter(q == rname & r == rname)

  # make the tri sequence identity plots
  df_d <- make_tri_df(df_rname_seq_ident)

  # make the histogram
  plot_hist <- make_hist(df_rname_seq_ident)

  # Filter data.
  # Get centromeric transition regions separately to outline.
  df_rm_sat_out <- df_rm_sat_out[order(df_rm_sat_out$region)] %>%
    filter(chr == rname)
  df_humas_hmmer_stv_out <- df_humas_hmmer_stv_out %>%
    filter(chr == rname)
  df_rm_sat_out_ct <- df_rm_sat_out %>%
    filter(region == "ct" & chr == rname)
  df_rm_sat_out_other <- df_rm_sat_out %>%
    filter(region != "ct" & chr == rname)

  segment_linewidth <- 10
  contig_len <- max(df_rname_seq_ident$q_en) - min(df_rname_seq_ident$q_st)
  # Calculated adjustment factor (y-px / 3.5mb) for segment y position.
  segment_y_adj_factor <- -0.07

  segment_y <- segment_y_adj_factor * contig_len
  ct_outline_edges_x <- 5000

  plot_ident_cen <- ggplot()

  # Also plot methyl binned if present.
  if (!typeof(df_methyl_binned) == "logical") {
    plot_methyl_binned <- ggplot() +
      geom_area(
        data = df_methyl_binned,
        aes(x = start2, y = meth_prob),
        fill = "black",
      ) +
      scale_x_continuous(expand=c(0,0)) +
      theme_cowplot() +
      theme(
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
      ) +
      # Negative margins??? ggplot >:(
      theme(
        plot.margin = margin(t = 20, b = -10, l = 7, r = 35)
      ) +
      ylab("Mean CpG\nmethylation (%)") +
      theme(
        axis.title.y = element_text(size = 12)
      )
  } else {
    plot_methyl_binned <- NA
  }

  # If df_cdr and methyl bins are provided, add lines for CDR on top.
  if (!(typeof(df_cdr) == "logical") && (!(typeof(plot_methyl_binned) == "logical")))  {
    plot_methyl_binned <- plot_methyl_binned +
      geom_segment(
        data = df_cdr,
        aes(x = start2, y = 100, xend = stop2, yend = 100),
        linewidth = 1,
        colour = "black",
      )
  }
  plot_ident_cen <- plot_ident_cen +
    # Add HOR orientation segments.
    geom_segment(
      data=df_hor_ort,
      aes(
        x = start2,
        xend = stop2,
        y = segment_y * 2,
        yend = segment_y * 2,
      ),
      linewidth = 1.5,
      lineend = "butt",
      linejoin = "mitre",
      # Point arrow to last position. See https://rdrr.io/r/grid/arrow.html
      arrow = arrow(length = unit(0.25, "cm"), ends = "last", type = "closed")
    ) +
    # Make larger ct segment as outline
    geom_segment(
      data = df_rm_sat_out_ct,
      aes(
        x = start2,
        y = segment_y,
        xend = stop2 + 1000,
        yend = segment_y,
      ),
      linewidth = segment_linewidth + 0.5
    ) +
    # Need to change colorscale as outline gets desaturated.
    scale_color_manual(values = c("black")) +
    new_scale_color() +
    geom_segment(
      data = df_rm_sat_out_ct,
      # Adjust x and y values
      # to be slightly smaller than the outline
      aes(
        x = start2 + ct_outline_edges_x,
        y = segment_y,
        xend = (stop2 - ct_outline_edges_x) + 1000,
        yend = segment_y,
        color = region
      ),
      key_glyph = "rect",
      linewidth = segment_linewidth
    ) +
    geom_segment(
      data = df_rm_sat_out_other,
      aes(
        x = start2,
        y = segment_y,
        xend = stop2 + 1000,
        yend = segment_y,
        color = region,
      ),
      key_glyph = "rect",
      linewidth = segment_linewidth
    ) +
    guides(
      color = guide_legend(nrow = 2, override.aes = list(linewidth=0.5))
    ) +
    scale_color_manual(values = get_rm_sat_annot_colors()) +
    labs(color="Sequence composition") +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 12)) +
    # New colorscale for hor monomers
    new_scale_color() +
    geom_segment(
      data = df_humas_hmmer_stv_out %>% filter(chr == rname),
      aes(
        x = start,
        xend = stop + 2000,
        y = segment_y,
        yend = segment_y,
        color = as.factor(mer),
      ),
      key_glyph = "rect",
      linewidth = segment_linewidth
    ) +
    guides(
      color = guide_legend(nrow = 5, override.aes = list(linewidth=0.5)),
    ) +
    scale_color_manual(values = get_humas_hmmer_stv_annot_colors()) +
    labs(color="Alpha-satellite HOR monomers")

  # New colorscale for stainedglass.
  plot_ident_cen <- plot_ident_cen +
    geom_polygon(
      df_d,
      mapping = aes(x = x, y = y, fill = new_discrete, group = group),
      show.legend = FALSE
    ) +
    scale_y_reverse() +
    scale_fill_manual(values = discrete_color_ranges) +
    scale_x_continuous(labels = make_scale, limits = c(0, NA)) +
    theme_cowplot() +
    # Adjust legend.
    theme(
      legend.position = "right",
      legend.background = element_blank(),
      legend.key = element_rect(color = "black", fill = NA),
    ) +
    # Adjust axes
    theme(
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.y = element_blank()
    ) +
    xlab("Position (Mbp)") +
    theme(
      axis.title.x = element_text(size = 12)
    )

  if (!typeof(df_methyl_binned) == "logical") {
    plot_ident_cen <- plot_ident_cen +
      # Expand margins to account for y-axis elements on methylation plot.
      theme(
        plot.margin = margin(t = 0, b = 0, l = 30, r = 7)
      )
  }

  plot_cen_legend <- get_legend(
    plot_ident_cen +
    theme(
      legend.justification="center",
      legend.box.just = "left",
    )
  )
  plot_ident_cen <- plot_ident_cen + theme(legend.position = "none")

  grid_bottom_row <- cowplot::plot_grid(
    plotlist = list(NULL, plot_cen_legend, plot_hist, NULL),
    nrow = 1,
    rel_widths = c(0.1, 1, 1, 0.1),
    labels = NULL
  )

  plot_title <- ggdraw() +
    draw_label(
      rname,
      fontface = 'bold',
      x = 0,
      hjust = 0
    ) +
    theme(
      # add margin on the left of the drawing canvas,
      # so title is aligned with left edge of first plot
      plot.margin = margin(7, 0, -7, 60)
    )

  if (!typeof(plot_methyl_binned) == "logical") {
    plot <- cowplot::plot_grid(
      plot_title, plot_methyl_binned, plot_ident_cen, grid_bottom_row,
      nrow=4,
      ncol=1,
      rel_heights = c(0.05, 0.3, 1.2, 0.8),
      labels = NULL
    )
  } else {
    plot <- cowplot::plot_grid(
      plot_title, plot_ident_cen, grid_bottom_row,
      nrow=3,
      ncol=1,
      rel_heights = c(0.05, 1.2, 0.8),
      labels = NULL
    )
  }

  return(plot)
}
