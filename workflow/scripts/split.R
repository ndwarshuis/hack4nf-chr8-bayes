library(tidyverse)

exp_df <- readr::read_tsv(
  snakemake@input[[1]],
  ## "../../results/depmap/dose_mfi.tsv.gz",
  col_types = cols(
    dose = "d",
    viability = "d",
    cell_line = "f",
    name = "f",
    .default = "-"
  )
) %>%
  ## slice_head(n = 1000) %>%
  group_by(cell_line, name) %>%
  mutate(id = cur_group_id()) %>%
  ungroup() 

output <- snakemake@output[[1]]

# make sure id is first column
exp_df %>%
  select(id, cell_line, name) %>%
  unique() %>%
  readr::write_tsv(output)

exp_lst <- exp_df %>%
  select(-cell_line, -name) %>%
  group_by(id) %>%
  group_walk(
    ~ readr::write_tsv(
      .x,
      file.path(dirname(output), sprintf("%s.tsv", .y$id))
    )
  )
