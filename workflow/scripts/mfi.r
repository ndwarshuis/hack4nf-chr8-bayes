library(tidyverse)
library(R2jags)

# TODO add r^2
train_model <- function(df, key) {
    mod_spec <- function(){
        # Likelihood:
        for (i in 1:N){
            b_resp[i] ~ dnorm(mu[i], tau) # tau is precision (1 / variance)
            mu[i] <- 1 / (1 + 2 ^ (-kappa * (b_conc[i] - alpha)))
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

    init_values <- function(){
      list(alpha = rnorm(1),
           kappa = dgamma(2, 1),
           sigma = runif(0.1))
    }

    jags_data <- list(b_resp = df$viability,
                     b_conc = df$dose,
                     N = nrow(df))

    fit <- jags(data = jags_data,
                inits = init_values,
                parameters.to.save = c("alpha", "kappa", "sigma"),
                model.file = mod_spec,
                n.chains = 3,
                n.iter = 12000,
                n.burnin = 2000,
                n.thin = 10,
                DIC = F)
    fit
}

treatment_info <- readr::read_csv(
  "secondary-screen-replicate-treatment-info.csv",
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
  "secondary-screen-pooling-info.csv",
  col_types = cols(row_name = "f",
                   screen_id = "f",
                   detection_pool = "f",
                   pool_id = "f",
                   .default = "-")) %>%
  rename(cell_line = row_name)

mfi_df <- readr::read_csv(
  "secondary-screen-mfi.csv",
  col_types = cols("...1" = "f", .default = "d")
) %>%
  rename(cell_line = "...1") %>%
  pivot_longer(cols = -cell_line, names_to = "condition_id", values_to = "MFI") %>%
  mutate(condition_id = factor(condition_id))


chr8_df <- readr::read_csv(
  "8q_Arm_level_CNAs.csv",
  col_types = cols("DepMap ID" = "f",
                   "8q  Arm-level CNAs" = "d",
                   "Lineage" = "f",
                   .default = "-")) %>%
  rename("cell_line" = "DepMap ID",
         "gain" = "8q  Arm-level CNAs")

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

exp_df <- filter(ndata_df, perturbation_type == "experimental_treatment") %>%
  left_join(ctl_df, by = c("cell_line", "screen_id", "detection_plate")) %>%
  mutate(viability = (nMFI - positive_control) / (vehicle_control - positive_control)) %>%
  filter(is.finite(viability)) %>%
  select(cell_line, condition_id, screen_id, detection_plate, name, dose, viability) %>%
  left_join(chr8_df, by = "cell_line") %>%
  ## filter(Lineage %in% c("Skin", "Fibroblast", "Peripheral Nervous System"))
  filter(Lineage %in% c("Peripheral Nervous System"))

exp_lst <- exp_df %>%
  group_by(cell_line, name) %>%
  group_map(~ list(df = .x, cell_line = .y$cell_line, drug = .y$name)) %>%
  head() %>%
  map(~ c(.x, list(fit = train_model(.x$df))))

extract_summary_col <- function(x, name) {
  set_names(x$fit$BUGSoutput$summary[,name], ~ sprintf("%s_%s", .x, name))
}

summary_df <- exp_lst %>%
  map_dfr(~ c(list(cell_line = .x$cell_line,
                   drug = .x$drug),
              extract_summary_col(.x, "mean"),
              extract_summary_col(.x, "sd"),
              extract_summary_col(.x, "Rhat"),
              extract_summary_col(.x, "n.eff")))
