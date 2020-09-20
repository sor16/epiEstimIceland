library(tidyverse)
library(rstan)
library(lubridate)
library(readxl)
library(here)
library(cmdstanr)
library(posterior)

options(mc.cores = parallel::detectCores())
parallel:::setDefaultClusterOptions(setup_strategy = "sequential")

Make_EpiEstim_Model <- function(fit_date = Sys.Date(),use_quarantine=T, warmup = 500, iters = 500, chains = 4) {
    
    d <- read_csv("https://docs.google.com/spreadsheets/d/1xgDhtejTtcyy6EN5dbDp5W3TeJhKFRRgm6Xk0s0YFeA/export?format=csv&gid=1788393542") %>%
        select(date = Dagsetning, local = Innanlands_Smit,border_1=Landamaeri_Smit_1,border_2=Landamaeri_Smit_2, imported = Innflutt_Smit,prop_quarantine=Hlutf_Sottkvi,num_quarantine=Fjoldi_Sottkvi) %>% 
        mutate(date = ymd(date),
               total = local + imported,
               prop_quarantine=if_else(use_quarantine & total!=0,(num_quarantine+border_1+border_2)/total,0)) %>% 
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
    prop_quarantine <- d$prop_quarantine

    
    stan_data <- list(
        N_days = N_days,
        SI_shape = 1.54,
        SI_rate = 0.28,
        local = local,
        total = total,
        prop_quarantine=prop_quarantine
    )
    
    mod <- cmdstan_model(here("Stan", "EpiEstim.stan"))
    
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
    
    fit$save_object(file = here("Results", "Models", "EpiEstim", str_c("Model_", fit_date, if_else(use_quarantine,"_w_quarantine","") ,".rds")))
    
    
    
}
