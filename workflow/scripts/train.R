library(tidyverse)
library(R2jags)
library(furrr)

cachedir <- snakemake@params[["cachedir"]]
logdir <- file.path(cachedir, "log")
sumdir <- file.path(cachedir, "sums")
simdir <- file.path(cachedir, "sims")

walk(
  c(logdir, sumdir, simdir),
  ~ dir.create(.x, showWarnings = FALSE, recursive = TRUE)
)

extract_summary_col <- function(fit, name) {
  set_names(fit$BUGSoutput$summary[, name], ~ sprintf("%s_%s", .x, name))
}

# TODO add r^2
train_model <- function(df) {

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
    b_resp = df[,"viability"],
    b_conc = df[,"dose"],
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

write_model <- function(df, .id) {
  # set up logging locally
  logfile <- file.path(logdir, sprintf("%s.txt", .id))
  file.create(logfile)
  con <- file(logfile)
  sink(con, append=TRUE)
  sink(con, append=TRUE, type="message")

  mod <- train_model(df)

  as_tibble(as.list(mod$summary)) %>%
    mutate(id = .id) %>%
    readr::write_tsv(file.path(sumdir, sprintf("%s.tsv", .id)))

  as_tibble(mod$sims) %>%
    mutate(id = .id) %>%
    readr::write_tsv(file.path(simdir, sprintf("%s.tsv", .id)))

  sink(NULL, append=TRUE)
  sink(NULL, append=TRUE, type="message")
  close(con)
}

plan(multicore, workers = snakemake@threads)

id_df <- readr::read_tsv(
  snakemake@input[["ids"]],
  lazy = TRUE,
  ## n_max = 100,
  col_types = cols(
    id = "d",
    .default = "-"
  )
)

all_ids <- id_df$id

# ASSUME summary and sim matrix will be in sync
current_ids <- list.files(sumdir, full.names = TRUE) %>%
  map(~ as.integer(str_extract(basename(.x), "[0-9]+")))

# TODO use md5s here to really figure out what needs to be (re)done?
ids_to_train <- setdiff(all_ids, current_ids)

message(sprintf("ids to be trained: %s\n", str_c(ids_to_train, collapse = ", ")))

dose_df <- readr::read_tsv(
  snakemake@input[["doses"]],
  lazy = TRUE,
  col_types = cols(
    id = "i",
    dose = "d",
    viability = "d",
    .default = "-"
  )
) %>%
  right_join(tibble(id = ids_to_train), by = "id") %>%
  nest(.by = id) %>%
  mutate(data = map(data, as.matrix))

future_walk2(
  dose_df$data,
  as.list(dose_df$id),
  write_model,
  .options = furrr_options(seed = 123, stdout = FALSE)
)

# cat all cache files to satisfy snakemake
list.files(sumdir, full.names = TRUE) %>%
  map_dfr(~ readr::read_tsv(.x, col_types = "d")) %>%
  readr::write_tsv(snakemake@output[["summary"]])

list.files(simdir, full.names = TRUE) %>%
  map_dfr(~ readr::read_tsv(.x, col_types = "d")) %>%
  readr::write_tsv(snakemake@output[["sims"]])

# remove all cache files assuming we succeeded
list.files(sumdir, full.names = TRUE) %>%
  c(list.files(simdir, full.names = TRUE)) %>%
  walk(file.remove)
