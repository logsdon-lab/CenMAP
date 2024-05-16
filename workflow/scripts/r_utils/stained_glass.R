
make_scale = function(vals){
  comma(vals/1e6)
}

make_k = function(vals){
  comma(vals/1e3)
}

get_colors = function(sdf){
  bot = floor(min(sdf$perID_by_events)); top = 100
  breaks = seq(70,100,1)
  labels = seq(length(breaks)-1)
  return( cut(sdf$perID_by_events, breaks=seq(70,100,1), labels=seq(length(breaks)-1), include.lowest=TRUE)  )

}

read_bedpe <- function(all.files){
  l <- lapply(all.files, fread, sep="\t")
  df <- rbindlist( l )

  df$discrete = get_colors(df)
  if("#query_name" %in% colnames(df)){
    df$q = df$`#query_name`
    df$q_st = df$query_start
    df$q_en = df$query_end
    df$r = df$reference_name
    df$r_st = df$reference_start
    df$r_en = df$reference_end
  }
  window=max(df$q_en-df$q_st)
  df$first_pos = df$q_st/window
  df$second_pos = df$r_st/window

  df <- df %>% mutate(
    new_discrete=case_when(
      perID_by_events > 0 & perID_by_events < 80 ~ "1-10",
      perID_by_events >= 80 & perID_by_events < 90 ~ "11-20",
      .default = as.character(discrete)
    )
  )
  return(df)
}

discrete_color_ranges <- c(
  "#4b3991", #1-10
  "#2974af", #11-20
  "#4a9da8", #21
  "#57b894", #22
  "#9dd893", #23
  "#e1f686", #24
  "#ffffb2", #25
  "#fdda79", #26
  "#fb9e4f", #27
  "#ee5634", #28
  "#c9273e", #29
  "#8a0033"  #30
)

# get the lowest 0.1% of the data so we can not plot it
make_hist = function(sdf){
  bot = quantile(sdf$perID_by_events, probs=0.001)[[1]]
  p <- ggplot(sdf) +
    geom_histogram(aes(perID_by_events, fill = new_discrete), bins = 300) +
    scale_fill_manual(values = discrete_color_ranges) +
    theme_cowplot() +
    theme(legend.position = "none") +
    scale_y_continuous(labels=make_k) +
    coord_cartesian(xlim = c(bot, 100))+
    xlab("% identity") +
    ylab("# of alignments (thousands)")
  p
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
  df
}

make_tri = function(sdf, rname=""){
  sdf$w = (sdf$first_pos + sdf$second_pos)
  sdf$z = -sdf$first_pos + sdf$second_pos
  window = max(sdf$q_en - sdf$q_st)
  tri_scale =  max(sdf$q_st) / max(sdf$w)
  sdf$window <- max(sdf$q_en - sdf$q_st) / tri_scale
  sdf$group <- seq(nrow(sdf))
  df_d <- rbindlist(apply(sdf, 1, diamond))

  ggplot(df_d) +
    geom_polygon(aes(x = w * tri_scale, y = z * window , fill = new_discrete, group = group)) +
    theme_cowplot() +
    scale_fill_manual(values = discrete_color_ranges) +
    scale_x_continuous(labels=make_scale, limits = c(0,NA)) +
    scale_y_continuous(labels=make_scale, limits = c(0,NA)) +
    xlab("Genomic position (Mbp)") +
    ylab("") +
    theme(
      legend.position = "none",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank()
      ) +
    ggtitle(rname)
}

make_dot = function(sdf, rname=""){
  max = max(sdf$q_en, sdf$r_en)
  window = max(sdf$query_end - sdf$query_start)
  ggplot(sdf) +
    geom_tile(aes(x = q_st, y = r_st, fill = new_discrete, height=window, width=window)) +
    theme_cowplot() +
    scale_fill_manual(values = discrete_color_ranges) +
    theme(legend.position = "none") +
    scale_x_continuous(labels=make_scale, limits = c(0,max)) +
    scale_y_continuous(labels=make_scale, limits = c(0,max)) +
    coord_fixed(ratio = 1) +
    facet_grid( r ~ q )+
    xlab("Genomic position (Mbp)") +
    ylab("") +
    ggtitle(rname)
}


make_plots <- function(r_name) {
  sdf = copy(df[r == r_name & q == r_name])

  # make the plots
  p_lone = make_tri(sdf, rname=r_name)
  scale = 1 / 2.25
  p_hist = make_hist(sdf)

  # setup the space
  dir.create(glue("{OUT}/pdfs/{r_name}/"))
  dir.create(glue("{OUT}/pngs/{r_name}/"))

  # save the plots
  p_hist <- p_hist +
      theme(text = element_text(size = 8), axis.text = element_text(size = 8))

  # ranges for inset hist
  mmax <- max(sdf$q_en, sdf$r_en)
  build <- ggplot_build(p_lone)
  yr <- build$layout$panel_params[[1]]$y.range
  xmin <- -1 / 10 * mmax
  xmax <- mmax * 1 / 3.75
  ymin <- yr[2] * 1 / 2
  ymax <- yr[2] * 2 / 2
  # combine
  plot <- p_lone + annotation_custom(
    ggplotGrob(p_hist),
    xmin = xmin, xmax = xmax,
    ymin = ymin, ymax = ymax
  )
  ggsave(
    plot = plot,
    file = glue("{OUT}/pdfs/{r_name}/{PRE}__{r_name}__tri.pdf"),
    height = 12 * scale,
    width = 9
  )

  ggsave(
    plot = plot,
    file = glue("{OUT}/pngs/{r_name}/{PRE}__{r_name}__tri.png"),
    height = 12 * scale,
    width = 9,
    dpi = DPI
  )
  return (plot)
}
