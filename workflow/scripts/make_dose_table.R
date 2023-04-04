library(tidyverse)

treatment_info <- readr::read_csv(
  snakemake@input[["treat"]],
  ## "../../resources/depmap/treatments.csv",
  col_types = cols(column_name = "f",
                   broad_id = "f",
                   screen_id = "f",
                   dose = "d",
                   detection_plate = "f",
                   perturbation_type = "f",
                   compound_plate = "f",
                   name = "f",
                   well = "f",
                   .default = "-")) %>%
  rename(condition_id = column_name)

pooling_info <- readr::read_csv(
  snakemake@input[["pool"]],
  ## "../../resources/depmap/pool_info.csv",
  col_types = cols(row_name = "f",
                   screen_id = "f",
                   detection_pool = "f",
                   pool_id = "f",
                   .default = "-")) %>%
  rename(cell_line = row_name)

mfi_df <- readr::read_csv(
  snakemake@input[["mfi"]],
  col_types = cols("...1" = "f", .default = "d")
) %>%
  rename(cell_line = "...1") %>%
  pivot_longer(cols = -cell_line, names_to = "condition_id", values_to = "MFI") %>%
  mutate(condition_id = factor(condition_id))


chr8_df <- readr::read_csv(
  snakemake@input[["cna"]],
  col_types = cols("...1" = "f", "8q" = "d", .default = "-")) %>%
  rename(cell_line = "...1",
         "gain" = "8q")

cell_df <- readr::read_csv(
  snakemake@input[["cell"]],
  col_types = cols("row_name" = "c",
                   "primary_tissue" = "f",
                   .default = "-")) %>%
  rename(cell_line = "row_name")

data_df <- mfi_df %>%
  left_join(treatment_info, by = "condition_id") %>%
  left_join(pooling_info, by = c("cell_line", "screen_id")) %>%
  filter(is.finite(MFI))

barcodes <- data_df %>%
  filter(pool_id == "CP01") %>%
  group_by(condition_id) %>%
  summarize(barcode_MFI = median(MFI))

ndata_df <- data_df %>%
  filter(pool_id != "CP01") %>%
  left_join(barcodes, by = "condition_id") %>%
  mutate(nMFI = MFI / barcode_MFI) %>%
  select(-MFI, -barcode_MFI)

ctl_df <- ndata_df %>%
  filter(perturbation_type != "experimental_treatment") %>%
  group_by(cell_line, screen_id, detection_plate, perturbation_type) %>%
  summarize(nMFI = median(nMFI), .groups = "drop") %>%
  pivot_wider(id_cols = c(cell_line, screen_id, detection_plate),
              names_from = perturbation_type,
              values_from = nMFI)

filter(ndata_df, perturbation_type == "experimental_treatment") %>%
  left_join(ctl_df, by = c("cell_line", "screen_id", "detection_plate")) %>%
  mutate(viability = (nMFI - positive_control) / (vehicle_control - positive_control)) %>%
  filter(is.finite(viability)) %>%
  select(cell_line, condition_id, screen_id, detection_plate, name, dose, viability) %>%
  left_join(chr8_df, by = "cell_line") %>%
  left_join(cell_df, by = "cell_line") %>%
  readr::write_tsv(snakemake@output[[1]])

