library(tidyverse)
library(R2jags)
library(furrr)

extract_summary_col <- function(fit, name) {
  set_names(fit$BUGSoutput$summary[, name], ~ sprintf("%s_%s", .x, name))
}

# TODO add r^2
train_model <- function(df, key) {
  mod_spec <- function() {
    # Likelihood:
    for (i in 1:N) {
      b_resp[i] ~ dnorm(mu[i], tau) # tau is precision (1 / variance)
      mu[i] <- 1 / (1 + 2^(-kappa * (b_conc[i] - alpha)))
    }
    # Priors (totally made up...ssshh)
    alpha ~ dnorm(au, a_sd) # ec50
    au ~ dnorm(0, 1)
    as ~ dnorm(0.01, 5)
    a_sd <- 1 / (as * as)
    kappa ~ dgamma(kr, ks) # slope
    kr ~ dnorm(2, 0.01)
    ks ~ dnorm(1, 0.01)
    # response
    sigma ~ dunif(0.01, 0.1) # standard deviation
    tau <- 1 / (sigma * sigma) # sigma^2 doesn't work in JAGS
  }

  init_values <- function() {
    list(
      alpha = rnorm(1),
      kappa = dgamma(2, 1),
      sigma = runif(0.1)
    )
  }

  jags_data <- list(
    b_resp = df$viability,
    b_conc = df$dose,
    N = nrow(df)
  )

  fit <- jags(
    data = jags_data,
    inits = init_values,
    parameters.to.save = c("alpha", "kappa", "sigma"),
    model.file = mod_spec,
    n.chains = 3,
    n.iter = 12000,
    n.burnin = 2000,
    n.thin = 10,
    DIC = F
  )
  list(
    summary = c(
      extract_summary_col(fit, "mean"),
      extract_summary_col(fit, "sd"),
      extract_summary_col(fit, "Rhat"),
      extract_summary_col(fit, "n.eff")
    ),
    sims = fit$BUGSoutput$sims.array
  )
}

plan(multisession, workers = snakemake@threads)
## plan(multisession, workers = 8)

exp_df <- readr::read_tsv(
  snakemake@input[[1]],
  ## "../../results/depmap/dose_mfi.tsv.gz",
  col_types = cols(
    dose = "d",
    viability = "d",
    gain = "d",
    cell_line = "f",
    name = "f",
    .default = "-"
  )
) %>%
  group_by(cell_line, name) %>%
  mutate(id = cur_group_id()) %>%
  filter(gain == snakemake@wildcards[["gain"]]) %>%
  ungroup()

exp_df %>%
  select(cell_line, name, id) %>%
  unique() %>%
  readr::write_tsv(snakemake@output[["ids"]])

exp_lst <- exp_df %>%
  select(-cell_line, -name) %>%
  ## slice_head(n = 1000) %>%
  group_by(id) %>%
  group_map(~ list(df = .x, id = .y$id)) %>%
  ## head() %>%
  future_map(~ c(id = .x$id, train_model(.x$df)),
             .options = furrr_options(seed = 123, stdout = FALSE))

exp_lst %>%
  map_dfr(~ c(list(id = .x$id), .x$summary)) %>%
  readr::write_tsv(snakemake@output[["summary"]])

exp_lst %>%
  map_dfr(~ c(list(id = .x$id), as_tibble(.x$sims))) %>%
  readr::write_tsv(snakemake@output[["sims"]])
