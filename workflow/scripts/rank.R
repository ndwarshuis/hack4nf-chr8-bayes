library(tidyverse)
library(furrr)

plan(multicore, workers = snakemake@threads)

ref_id <- snakemake@params$refid

id_df <- readr::read_tsv(
  snakemake@input$ids,
  col_types = cols(
    id = "i",
    gain = "i",
    .default = "c",
  )
)

param_df <- readr::read_tsv(
  snakemake@input$summary,
  col_types = cols(id = "i", .default = "d")
) %>%
  left_join(id_df, by = "id") %>%
  select(-gain, -primary_tissue)

sims_df <- readr::read_tsv(
  snakemake@input$sims,
  col_types = cols(id = "i", .default = "d")
)

# filter out id's that have not converged
good_param_df <- param_df %>%
  filter(alpha_Rhat < 1.1) %>%
  filter(beta_Rhat < 1.1) %>%
  filter(kappa_Rhat < 1.1) %>%
  filter(sigma_Rhat < 1.1)

alpha_sim <- sims_df %>%
  select("1.alpha", id) %>%
  as.matrix()

sample_sim <- function(this_id) {
  sample(alpha_sim[alpha_sim[,"id"] == this_id, "1.alpha"], 1000, replace = TRUE)
}

ref_df <- good_param_df %>%
  filter(cell_line == ref_id) %>%
  rename_with(~ paste0(.x, "_ref"), c(-cell_line, -name)) %>%
  select(-cell_line) %>%
  mutate(alpha_dist_ref = future_map(
    id_ref,
    sample_sim,
    .options = furrr_options(seed = 123, stdout = FALSE)
  ))

test_df <- good_param_df %>%
  filter(!cell_line == ref_id) %>%
  mutate(alpha_dist = future_map(
    id,
    sample_sim,
    .options = furrr_options(seed = 123, stdout = FALSE)
  ))

full_test_df <- ref_df %>%
  left_join(test_df, by = "name") %>%
  filter(!is.na(cell_line)) %>%
  mutate(diff = map2(alpha_dist_ref, alpha_dist, ~ .x - .y))

# only keep drugs with all cell lines intact
full_drugs <- full_test_df %>%
  group_by(name) %>%
  tally() %>%
  filter(n == max(n)) %>%
  pull(name)

rank_df <- full_test_df %>%
  select(name, cell_line, diff) %>%
  filter(name %in% full_drugs) %>%
  unnest(diff) %>%
  # ignore cell_line to weight each equally
  group_by(name) %>%
  summarize(risk = mean(diff < 0),
            diff = mean(diff)) %>%
  arrange(desc(diff), risk)

rank_df %>%
  readr::write_tsv(snakemake@output[["table"]])

rank_df %>%
  filter(risk < 0.1) %>%
  ggplot(aes(risk, diff)) +
  geom_point() +
  labs(x = "Risk (p[Diff EC50 < 0])",
       y = "Diff EC50 (log10[Ref] - log10[Test])")
ggsave(snakemake@output[["chart"]])
