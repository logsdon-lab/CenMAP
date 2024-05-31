library(dplyr)
library(ggplot2)
library(data.table)
library(tidyr)
library(stringr)
library(argparser)


axis_lbl_suffix <- "of centromeres completely and accurately assembled"
read_cen_counts <- function(fpath, label) {
  df <- fread(
    fpath,
    sep = "\t",
    col.names = c("sample", "cnt", "mean")
  )
  df <- df %>% mutate(lbl = label)
  return(df)
}

# https://stackoverflow.com/a/68009651
reorder_where <- function (x, by, where, fun = mean, ...) {
  xx <- x[where]
  byby <- by[where]
  byby <- tapply(byby, xx, FUN = fun, ...)[x]
  reorder(x, byby)
}

plot_cen_counts <- function(df, means, lbl_colors, plt_title = NULL, hline_colors = NULL, reorder_lbl = NULL) {
  if (is.null(reorder_lbl)) {
    plt <- ggplot(df, aes(fill = lbl, y = mean, x = sample))
  } else {
    plt <- ggplot(df, aes(fill = lbl, y = mean, x = reorder_where(sample,-mean, lbl==reorder_lbl)))
  }
  if (is.null(hline_colors)) {
    hline_colors <- lbl_colors
  }

  plt <- plt +
    geom_bar(position = position_dodge2(), stat = "identity", color = "black", linewidth=0.1) +
    scale_fill_manual("Label", values = lbl_colors) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title.x = element_blank()
    ) +
    ylab(paste("%", axis_lbl_suffix)) +
    scale_y_continuous(
      # Prevent spacing between bottom of bar and x-axis.
      # https://stackoverflow.com/a/63544058
      expand = expansion(mult = c(0, .1)),
      limits = c(0, 100),
      # Secondary axis for number of centromeres
      sec.axis = sec_axis(
        ~ . * 0.46,
        name = paste("#", axis_lbl_suffix),
        breaks = c(seq(0, 40, by = 10), 46)
      )
    ) +
    geom_hline(
      yintercept = means,
      linetype = "dashed",
      color = hline_colors
    ) +
    annotate(
      "label",
      x = length(unique(df$sample)) - 5, y = 100,
      label = paste0("Mean = ", round(means, 1), "%", collapse = "\n"),
      size = 5
    ) +
    theme(legend.position = "bottom")

  if (!is.null(plt_title)) {
    plt <- plt + ggtitle(plt_title)
  }
  return(plt)
}

p <- arg_parser("Plot centromere counts by sample.")
p <- add_argument(
  p, "--input",
  help = "Input centromere HOR array lengths.",
  type = "character"
)
p <- add_argument(
  p, "--label",
  help = "Label for samples.",
  type = "character",
  default = NULL
)
p <- add_argument(
  p, "--color",
  help = "Bar plot color",
  type = "character",
  default = "#DA8B26"
)
p <- add_argument(
  p, "--output",
  help = "Output plot file.",
  type = "character"
)
argv <- parse_args(p)

if (is.na(argv$label)) {
  lbl <- "N/A"
  plt_title <- "Centromeres Correct Assembled"
} else {
  lbl <- argv$label
  plt_title <- paste0("Centromeres Correct Assembled for ", lbl)
}

df <- read_cen_counts(argv$input, lbl) %>%
  filter(sample != "all") %>%
  select(-c(cnt))

plt <- plot_cen_counts(
  df,
  c(mean(df$mean)),
  argv$color,
  reorder_lbl = lbl,
  hline_colors = c("black"),
  plt_title = plt_title
) + theme(legend.position = "none")

ggsave(argv$output, plot = plt, dpi = 300, width = 16, height = 8)
