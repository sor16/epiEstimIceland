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
make_preds <- function(draws,N_days, SI, pop=356991, ...) {
    for (t in seq(N_days, nrow(draws))) {
        idx_range <- (t - 1):(t - length(SI)+1)
        # Population of Iceland
        pop = 356991
        draws$lambda[t] <- t((1-draws$prop_quarantine[idx_range])*draws$mu_hat[idx_range]) %*% head(SI,-1)
        draws$lambda_q[t] <- t(draws$prop_quarantine[idx_range]*draws$mu_hat[idx_range]) %*% head(SI,-1)
        # S represents the proportion of people exposed to the virus
        S_t = if_else(sum(draws$mu_hat) < 356991, (1-sum(draws$mu_hat)/356991),0)
        draws$mu_hat[t] <- S_t*(draws$lambda[t] * draws$R[t] + draws$lambda_q[t] * draws$R_q[t])# + draws$mu_hat[t])## add artificially imported infections already stored in mu_hat[t]
    }
    
    return(draws)
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

calculate_lambda <- function(cases,SI,p_vec){
    N_days <- length(cases)
    lambda <- numeric(N_days)
    for (t in 2:N_days) {
        lambda[t] <- t(head(p_vec*cases, t - 1)) %*% tail(rev(SI), t - 1) / sum(tail(rev(SI), t - 1))
    }
    return(lambda)
}

scenario <- function(m,d,R_fun,prop_imported=0,future_prop_quarantine=rep(0,pred_days-1),num_border_tests=2000,pred_days=42,use_quarantine=T, fit_date=Sys.Date(), future_border_cases=rep(0,pred_days-1)) {
    d <- d %>% filter(date < fit_date)
    N_days <- nrow(d)
    SI <- get_SI_vec(N_days)
    lambda_vec <- calculate_lambda(d$total,SI,1-d$prop_quarantine)
    lambda_q_vec <- calculate_lambda(d$total,SI,d$prop_quarantine) 
    R_q_draws <- spread_draws(m,log_R_q) %>%
                    mutate(iter = row_number()) %>%
                    mutate(R_q=exp(log_R_q)) %>%
                    select(iter,R_q)
    R_draws <- spread_draws(m, R[day]) %>% 
                group_by(day) %>% 
                mutate(iter = row_number()) %>%
                ungroup() %>% 
                select(iter, day, R) %>%
                inner_join(R_q_draws,by='iter')
    N_iter <- max(R_draws$iter)
    last_R <- R_draws %>% 
              filter(day == max(day)) %>% 
              .$R
    last_R_q <- R_draws %>% 
                filter(day == max(day)) %>% 
                .$R_q
    future_R <- crossing(day = max(R_draws$day) + seq_len(pred_days),iter = seq_len(N_iter)) %>% 
                group_by(iter) %>% 
                mutate(R = R_fun(last_R[iter],day-N_days),
                       R_q = last_R_q[iter]) %>%
                       #R = as.numeric(R) + cumsum(rnorm(n(), sd = 0.02))) %>% 
                ungroup()
    
    posterior_dat <- R_draws %>% 
                bind_rows(future_R) %>% 
                group_by(iter) %>% 
                mutate(prop_quarantine=c(d$prop_quarantine,future_prop_quarantine),
                       lambda = c(lambda_vec, rep(0,pred_days-1)),
                       lambda_q = c(lambda_q_vec, rep(0,pred_days-1)),
                       mu_hat = (day<N_days)*(R * lambda + R_q*lambda_q)) %>% #*as.numeric(rbinom(n(),size=num_border_tests,prob = prop_imported))) %>% # 
                # filter(iter %in% 1:100) %>% 
                group_modify(make_preds, N_days=N_days,SI=SI) %>% 
                ungroup %>% 
                mutate(y_hat = rnbinom(n(), mu = mu_hat, size = m$phi[iter])) %>% 
                pivot_longer(c(-iter, -day, -lambda, -mu_hat))
    return(posterior_dat)
}

summarize_posterior_dat <- function(posterior_dat) {
    plot_dat <- posterior_dat %>% 
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
