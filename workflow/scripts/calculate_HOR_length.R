library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(argparser, quietly = TRUE)


p <- arg_parser("Estimate HOR array length from HumAS-HMMER output.")
p <- add_argument(p, "--input",
  help = "Input bed file made from HumAS-HMMER output.",
  type = "character"
)
p <- add_argument(p, "--output",
  help = "Output bed file with chr, start, stop, and len columns.",
  type = "character"
)
p <- add_argument(p, "--bp_jump_thr",
  help = "Base pair jump threshold to group by",
  type = "numeric",
  default = 20000
)

argv <- parse_args(p)
argv$input <- "data/annotations/AS-HOR-vs-chm13_cens_v18.correctcoords.stv_row.all.bed"


cols <- c("chr", "start", "stop", "hor", "strand")
df <- fread(
  argv$input,
  sep = "\t",
  stringsAsFactors = TRUE,
  fill = TRUE, quote = "",
  header = FALSE, select = c(1, 2, 3, 4, 6)
)
colnames(df) <- cols

hor_array_lengths <- list()
for (chr_name in unique(df$chr)) {
  chr_name = "chr8"
  df_chr <- df %>% filter(chr == chr_name) %>% mutate(idx=row_number())

  df_bp_jumps <- df_chr %>%
    # Take only live ALRs.
    filter(str_detect(hor, "L")) %>%
    mutate(diff = start - lag(start)) %>%
    filter(diff > argv$bp_jump_thr)

  # No large gaps. One HOR array.
  if (nrow(df_bp_jumps) == 0) {
    hor_array_lengths[[chr_name]] <- df_chr %>%
      filter(str_detect(hor, "L")) %>%
      summarize(
        chr_name = chr_name,
        start_pos = min(start),
        stop_pos = max(stop),
        len = max(stop) - min(start)
      )
    next
  }
  df_ranges <- rbind(
    data.frame(
      start_pos = c(min(df_chr$start)),
      stop_pos = c(first(df_bp_jumps$start))
    ),
    df_bp_jumps %>%
      mutate(
        start_pos = if_else(
          row_number() == 1,
          start,
          stop
        ),
        stop_pos = if_else(
          row_number() == n(),
          max(df_chr$stop),
          lead(start)
        )
      ) %>%
      select(start_pos, stop_pos)
  )
  elem_cnt <- 0
  ranges <- list()
  for (i in seq_len(nrow(df_bp_jumps))) {
    prev_range <- df_bp_jumps[i - 1]
    curr_range <- df_bp_jumps[i]
    next_range <- df_bp_jumps[i + 1]

    if (nrow(prev_range) == 0) {
      ranges[elem_cnt] <- c(min(df_chr$start), curr_range[["start"]])
      elem_cnt <- elem_cnt + 1
    }
    if (nrow(next_range) == 0) {
      ranges[elem_cnt] <- c(curr_range[["stop"]], max(df_chr$stop))
      elem_cnt <- elem_cnt + 1
    }

    ranges[elem_cnt] <- c(curr_range[["start"]], curr_range[["stop"]])
    elem_cnt <- elem_cnt + 1
  }

  # Each row should be contiguous.
  # If distance between prev row greater than 10bp, filtered out.

  # for (i in seq_len(nrow(df_ranges))) {
  #   start_pos <- df_ranges[[i]][1]
  #   stop_pos <- df_ranges[[i]][2]

  #   df_chr_section <- df_chr %>%
  #     filter(start >= start_pos & stop < stop_pos) %>%
  #     mutate(diff=replace_na(start - lag(stop), 0))
  #   print(df_chr_section)
  # }

  hor_array_lengths[[chr_name]] <- df_ranges %>%
    rowwise() %>%
    summarise(
      chr_name = chr_name,
      # Reset position with bed.
      start_pos = df_chr %>%
        filter(
          start >= start_pos & stop < stop_pos & replace_na(start - lag(stop), 0) < 10
        ) %>%
        pull(start) %>%
        first(),
      stop_pos = df_chr %>%
        filter(
          start >= start_pos & stop < stop_pos & replace_na(start - lag(stop), 0) < 10
        ) %>%
        pull(stop) %>%
        last(),
    ) %>%
    mutate(len = stop_pos - start_pos) %>%
    filter(
      len != 0
    )
}

df_final <- bind_rows(hor_array_lengths) %>%
  # Extract chr and assign value to sort by.
  mutate(chr_lbl = str_extract(chr_name, "chr([0-9XY]+)", group = 1)) %>%
  mutate(chr_lbl = case_when(
    chr_lbl == "X" ~ 23,
    chr_lbl == "Y" ~ 24,
    .default = as.numeric(chr_lbl)
  )) %>%
  arrange(chr_lbl) %>%
  select(-chr_lbl)

write.table(
  df_final,
  file = ifelse(
    is.na(argv$output),
    # Why are you like this R?
    # https://stackoverflow.com/a/66610704
    "",
    argv$output
  ),
  quote = FALSE,
  sep = "\t",
  col.names = FALSE,
  row.names = FALSE
)
