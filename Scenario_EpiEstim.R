library(cmdstanr)
library(rstan)
library(posterior)
library(broom.mixed)
library(tidyverse)
library(cowplot)
library(scales)
library(lubridate)
library(plotly)
library(cmdstanr)
library(posterior)
library(tidybayes)
library(here)
library(gganimate)
library(data.table)

#### Help functions ####
make_preds <- function(d,...) {
    for (t in seq(N_days, nrow(d))) {
        idx_range <- (t - 1):(t - length(SI)+1)
        d$lambda[t] <- t((1-d$prop_quarantine[idx_range])*d$mu_hat[idx_range]) %*% head(SI,-1)
        d$mu_hat[t] <- d$mu_hat[t]+d$lambda[t] * d$R[t] ## add artificially imported infections already stored in mu_hat[t]
    }
    return(d)
}

get_SI_vec <- function(N_days,shape=1.54,rate=0.28){
    SI_dat <- tibble(t = seq(1, N_days)) %>% 
        mutate(
            p = case_when(
                TRUE ~ pgamma(t + 0.5, shape = shape, rate = rate) - pgamma(t - 0.5, shape = shape, rate = rate)
            ),
            p = p / sum(p)
        )
    return(SI_dat$p)
}

calculate_lambda <- function(total,SI,prop_quarantine){
    N_days <- length(total)
    lambda <- numeric(N_days)
    for (t in 2:N_days) {
        lambda[t] <- t(head((1-prop_quarantine)*d$total, t - 1)) %*% tail(rev(SI), t - 1) / sum(tail(rev(SI), t - 1))
    }
    return(lambda)
}


scenario <- function(m,d,R_draws,R_fun,prop_imported,future_prop_quarantine=rep(0,pred_days-1),num_border_tests=2000,pred_days=42,use_quarantine=T){
    N_iter <- max(R_draws$iter)
    N_days <- nrow(d)
    last_R <- R_draws %>% 
              filter(day == max(day)) %>% 
              .$R
    future_R <- crossing(day = max(R_draws$day) + seq_len(pred_days),iter = seq_len(N_iter)) %>% 
                group_by(iter) %>% 
                mutate(R = R_fun(last_R[iter],day-N_days)) %>%
                       #R = as.numeric(R) + cumsum(rnorm(n(), sd = 0.02))) %>% 
                ungroup()
    
    plot_dat <- R_draws %>% 
                bind_rows(future_R) %>% 
                group_by(iter) %>% 
                mutate(prop_quarantine=c(d$prop_quarantine,future_prop_quarantine),
                       lambda = c(d$lambda,rep(0,pred_days-1)),
                       mu_hat = R * lambda+(day>N_days)*as.numeric(rbinom(n(),size=num_border_tests,prob = prop_imported))) %>% 
                # filter(iter %in% 1:100) %>% 
                group_by(iter) %>% 
                group_modify(make_preds) %>% 
                ungroup %>% 
                mutate(y_hat = rnbinom(n(), mu = mu_hat, size = m$phi[iter])) %>% 
                pivot_longer(c(-iter, -day, -lambda, -mu_hat)) %>% 
                group_by(day, name) %>% 
                summarise(lower_50 = quantile(value, 0.25),
                          upper_50 = quantile(value, 0.75),
                          lower_60 = quantile(value, 0.2),
                          upper_60 = quantile(value, 0.8),
                          lower_70 = quantile(value, 0.15),
                          upper_70 = quantile(value, 0.85),
                          lower_80 = quantile(value, 0.1),
                          upper_80 = quantile(value, 0.9),
                          lower_90 = quantile(value, 0.05),
                          upper_90 = quantile(value, 0.95),
                          lower_95 = quantile(value, 0.025),
                          upper_95 = quantile(value, 0.975)) %>% 
                pivot_longer(c(-day, -name), names_to = c("which", "prob"), names_sep = "_") %>% 
                pivot_wider(names_from = which, values_from = value) %>% 
                mutate(prob = parse_number(prob),
                       date = ymd("2020-02-28") + day)
    return(plot_dat)
}
