library(tidyverse)
library(R2jags)
library(furrr)


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
  fit
}

extract_summary_col <- function(x, name) {
  set_names(x$fit$BUGSoutput$summary[, name], ~ sprintf("%s_%s", .x, name))
}

plan(multicore, workers = snakemake@threads)

exp_df <- readr::read_tsv(
  snakemake@input[[1]],
  col_types = cols(
    dose = "d",
    viability = "d",
    gain = "d",
    .default = "c"
  )
)

exp_lst <- exp_df %>%
  filter(gain == 1) %>%
  group_by(cell_line, name) %>%
  group_map(~ list(df = .x, cell_line = .y$cell_line, drug = .y$name)) %>%
  head() %>%
  future_map(~ c(.x, list(fit = train_model(.x$df))),
             .options = furrr_options(seed = 123, stdout = FALSE))

exp_lst %>%
  map_dfr(~ c(
    list(
      cell_line = .x$cell_line,
      drug = .x$drug
    ),
    extract_summary_col(.x, "mean"),
    extract_summary_col(.x, "sd"),
    extract_summary_col(.x, "Rhat"),
    extract_summary_col(.x, "n.eff")
  )) %>%
  pivot_longer(cols = c(-cell_line, -drug)) %>%
  separate(name, c("var", "stat"), "_") %>%
  readr::write_tsv(snakemake@output[["summary"]])

exp_lst %>%
  map_dfr(~ c(
    list(
      cell_line = .x$cell_line,
      drug = .x$drug
    ),
    as_tibble(.x$fit$BUGSoutput$sims.array)
  )) %>%
  pivot_longer(cols = c(-cell_line, -drug)) %>%
  separate(name, c("chain", "var"), "\\.") %>%
  readr::write_tsv(snakemake@output[["sims"]])
