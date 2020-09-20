library(tidyverse)
library(rstan)
library(lubridate)
library(readxl)
library(here)
library(cmdstanr)
library(posterior)

options(mc.cores = parallel::detectCores())
parallel:::setDefaultClusterOptions(setup_strategy = "sequential")

Make_EpiEstim_Model <- function(fit_date = Sys.Date(), warmup = 500, iters = 500, chains = 4) {
    
    d <- read_csv("https://docs.google.com/spreadsheets/d/1xgDhtejTtcyy6EN5dbDp5W3TeJhKFRRgm6Xk0s0YFeA/export?format=csv&gid=1788393542") %>%
        select(date = Dagsetning, local = Innanlands_Smit, imported = Innflutt_Smit) %>% 
        mutate(date = ymd(date),
               total = local + imported) %>% 
        filter(date >= ymd("2020-02-28"))
    
    # shape <- 4.5
    # rate <- 0.6
    shape <- 1.54
    rate <- 0.28
    
    
    total <- d$total
    
    N_days <- length(total)
    
    local <- d$local
    total <- d$total
    N_days <- length(local)
    
    stan_data <- list(
        N_days = N_days,
        SI_shape = 1.54,
        SI_rate = 0.28,
        pi_shape = 2,
        pi_rate = 0.4,
        local = local,
        total = total
    )
    
    mod <- cmdstan_model(here("Stan", "epidemiaDiagnosis.stan"))
    
    fit <- mod$sample(
        data = stan_data, 
        show_messages = FALSE, 
        chains = chains, 
        parallel_chains = chains,
        iter_sampling = iters,
        iter_warmup = warmup,
        max_treedepth = 15,
        init = 0,
        refresh = 100
    )
    
    fit$save_object(file = here("Results", "Models", "EpiEstim", str_c("Model_", fit_date, ".rds")))
    
    
    
}