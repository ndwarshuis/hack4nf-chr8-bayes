con <- file(snakemake@log[[1]])
sink(con, append=TRUE)
sink(con, append=TRUE, type="message")

library(tidyverse)
library(R2jags)

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

exp_df <- readr::read_tsv(
  snakemake@input[[1]],
  col_types = cols(
    dose = "d",
    viability = "d",
    .default = "-"
  )
)

mod <- train_model(exp_df)

tibble(id = mod$id, mod$summary) %>%
  readr::write_tsv(snakemake@output[["summary"]])

as_tibble(mod$sims) %>%
  mutate(id = mod$id) %>%
  readr::write_tsv(snakemake@output[["sims"]])
