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

theme_set(theme_classic(base_size = 12) + 
              theme(legend.position = "none"))




mod <- read_rds(here("Results", "Models", "EpiEstim", str_c("Model_", Sys.Date(), ".rds")))
m <- mod$draws() %>% as_draws_df
rm(mod)

d <- read_csv("https://docs.google.com/spreadsheets/d/1xgDhtejTtcyy6EN5dbDp5W3TeJhKFRRgm6Xk0s0YFeA/export?format=csv&gid=1788393542") %>%
    select(date = Dagsetning, local = Innanlands_Smit, imported = Innflutt_Smit) %>% 
    mutate(date = ymd(date),
           total = local + imported) %>% 
    filter(date >= ymd("2020-02-28"))


icelandic_dates <- function(x) {
    months <- c("janúar", "febrúar", "mars", "apríl", "maí", "júní", 
                "júlí", "ágúst", "september", "október", "nóvember", "desember")
    
    paste0(mday(x), ". ", months[month(x)])
}

p1 <- d %>% 
    rename(Local = local, Imported = imported) %>% 
    pivot_longer(c(Local, Imported)) %>% 
    ggplot(aes(date, value, fill = name)) +
    geom_col(width = 1, alpha = 0.8) +
    geom_vline(xintercept = ymd("2020-03-16"), lty = 2) +
    geom_vline(xintercept = ymd("2020-03-24"), lty = 2) +
    geom_vline(xintercept = ymd("2020-05-04"), lty = 2) +
    geom_vline(xintercept = ymd("2020-06-15"), lty = 2) +
    geom_vline(xintercept = ymd("2020-07-31"), lty = 2) +
    geom_vline(xintercept = ymd("2020-08-19"), lty = 2) +
    scale_fill_brewer(type = "qual", palette = "Set1") +
    scale_x_date(breaks = c(ymd(c("2020-03-01", 
                                  "2020-03-16", "2020-03-24",
                                  "2020-05-04",
                                  "2020-06-15",
                                  "2020-07-31",
                                  "2020-08-16",
                                  "2020-09-01"))),
                 date_labels = "%B %d",
                 expand = expansion(add = 0),
                 limits = c(ymd("2020-02-27"), Sys.Date() + 1)) +
    scale_y_continuous(breaks = pretty_breaks(5), expand = expansion(mult = 0), limits = c(0, 110)) +
    labs(subtitle = "Incidence") +
    theme(axis.title = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          plot.margin = margin(5, 5, 5, 5),
          legend.position = c(0.28, 0.8), 
          legend.title = element_blank())

p2 <- spread_draws(m, R[day]) %>% 
    group_by(day) %>% 
    summarise(lower_50 = quantile(R, 0.25),
              upper_50 = quantile(R, 0.75),
              lower_60 = quantile(R, 0.2),
              upper_60 = quantile(R, 0.8),
              lower_70 = quantile(R, 0.15),
              upper_70 = quantile(R, 0.85),
              lower_80 = quantile(R, 0.1),
              upper_80 = quantile(R, 0.9),
              lower_90 = quantile(R, 0.05),
              upper_90 = quantile(R, 0.95),
              lower_95 = quantile(R, 0.025),
              upper_95 = quantile(R, 0.975)) %>% 
    pivot_longer(c(-day), names_to = c("which", "prob"), names_sep = "_") %>% 
    pivot_wider(names_from = which, values_from = value) %>% 
    mutate(prob = parse_number(prob),
           date = ymd("2020-02-28") + day) %>% 
    ggplot(aes(date, ymin = lower, ymax = upper)) +
    geom_ribbon(aes(fill = factor(-prob)), alpha = 0.7) +
    geom_vline(xintercept = ymd("2020-03-16"), lty = 2) +
    geom_vline(xintercept = ymd("2020-03-24"), lty = 2) +
    geom_vline(xintercept = ymd("2020-05-04"), lty = 2) +
    geom_vline(xintercept = ymd("2020-06-15"), lty = 2) +
    geom_vline(xintercept = ymd("2020-07-31"), lty = 2) +
    geom_vline(xintercept = ymd("2020-08-19"), lty = 2) +
    geom_hline(yintercept = 1, lty = 2) +
    scale_fill_brewer() +
    scale_x_date(breaks = c(ymd(c("2020-03-01", "2020-04-01",
                                  "2020-05-01", "2020-06-01", "2020-07-01",
                                  "2020-08-01"))),
                 date_labels = "%B %d",
                 expand = expansion(add = 0),
                 limits = c(ymd("2020-02-27"), Sys.Date() + 1)) +
    scale_y_continuous(breaks = pretty_breaks(8)) +
    ggtitle(label = waiver(),
            subtitle = latex2exp::TeX("$R_t$")) +
    theme(axis.title = element_blank(),
          plot.margin = margin(5, 5, 5, 8))



plot_grid(p1, p2, ncol = 1)


p3 <- spread_draws(m, y_hat[day]) %>% 
    group_by(day) %>% 
    summarise(lower_50 = quantile(y_hat, 0.25),
              upper_50 = quantile(y_hat, 0.75),
              lower_60 = quantile(y_hat, 0.2),
              upper_60 = quantile(y_hat, 0.8),
              lower_70 = quantile(y_hat, 0.15),
              upper_70 = quantile(y_hat, 0.85),
              lower_80 = quantile(y_hat, 0.1),
              upper_80 = quantile(y_hat, 0.9),
              lower_90 = quantile(y_hat, 0.05),
              upper_90 = quantile(y_hat, 0.95),
              lower_95 = quantile(y_hat, 0.025),
              upper_95 = quantile(y_hat, 0.975)) %>% 
    pivot_longer(c(-day), names_to = c("which", "prob"), names_sep = "_") %>% 
    pivot_wider(names_from = which, values_from = value) %>% 
    mutate(prob = parse_number(prob),
           date = ymd("2020-02-28") + day) %>% 
    ggplot(aes(date, ymin = lower, ymax = upper)) +
    geom_ribbon(aes(fill = factor(-prob)), alpha = 0.7) +
    geom_point(
        data = d, aes(x = date, y = local), inherit.aes = F
    ) +
    geom_vline(xintercept = ymd("2020-03-16"), lty = 2) +
    geom_vline(xintercept = ymd("2020-03-24"), lty = 2) +
    geom_vline(xintercept = ymd("2020-05-04"), lty = 2) +
    geom_vline(xintercept = ymd("2020-06-15"), lty = 2) +
    geom_vline(xintercept = ymd("2020-07-31"), lty = 2) +
    scale_fill_brewer() +
    scale_x_date(date_breaks = "month",
                 date_labels = "%B %d",
                 expand = expansion(add = 0),
                 limits = c(ymd("2020-02-27"), Sys.Date() + 1)) +
    scale_y_continuous(breaks = pretty_breaks(8), expand = expansion(mult = 0), limits = c(0, 145)) +
    labs(subtitle = "Retrofitted local cases") +
    theme(axis.title = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          plot.margin = margin(5, 5, 5, 11),
          legend.title = element_blank()) 

p4 <- spread_draws(m, R[day]) %>% 
    group_by(day) %>% 
    summarise(lower_50 = quantile(R, 0.25),
              upper_50 = quantile(R, 0.75),
              lower_60 = quantile(R, 0.2),
              upper_60 = quantile(R, 0.8),
              lower_70 = quantile(R, 0.15),
              upper_70 = quantile(R, 0.85),
              lower_80 = quantile(R, 0.1),
              upper_80 = quantile(R, 0.9),
              lower_90 = quantile(R, 0.05),
              upper_90 = quantile(R, 0.95),
              lower_95 = quantile(R, 0.025),
              upper_95 = quantile(R, 0.975)) %>% 
    pivot_longer(c(-day), names_to = c("which", "prob"), names_sep = "_") %>% 
    pivot_wider(names_from = which, values_from = value) %>% 
    mutate(prob = parse_number(prob),
           date = ymd("2020-02-28") + day) %>% 
    ggplot(aes(date, ymin = lower, ymax = upper)) +
    geom_ribbon(aes(fill = factor(-prob)), alpha = 0.7) +
    geom_hline(yintercept = 1, lty = 2) +
    geom_vline(xintercept = ymd("2020-03-16"), lty = 2) +
    geom_vline(xintercept = ymd("2020-03-24"), lty = 2) +
    geom_vline(xintercept = ymd("2020-05-04"), lty = 2) +
    geom_vline(xintercept = ymd("2020-06-15"), lty = 2) +
    geom_vline(xintercept = ymd("2020-07-31"), lty = 2) +
    scale_fill_brewer() +
    scale_x_date(date_breaks = "month",
                 date_labels = "%B %d",
                 expand = expansion(add = 0),
                 limits = c(ymd("2020-02-27"), Sys.Date() + 1)) +
    scale_y_continuous(breaks = pretty_breaks(8)) +
    ggtitle(label = waiver(),
            subtitle = latex2exp::TeX("$R_t$")) +
    theme(axis.title = element_blank(),
          plot.margin = margin(5, 5, 5, 14))

plot_grid(p3, p4, ncol = 1)

plot_dat <- spread_draws(m, R[day]) %>% 
    group_by(day) %>% 
    summarise(under_one = mean(R < 1)) %>% 
    mutate(date = ymd("2020-02-28") + day)

p4b <- plot_dat %>% 
    ggplot(aes(date, under_one)) +
    geom_area(fill = "#08519c", alpha = 1/8) +
    geom_area(data = plot_dat %>% mutate(under_one = ifelse(under_one < 15/16, 0, under_one)),
              fill = "#08519c", alpha = 1/8) +
    geom_area(data = plot_dat %>% mutate(under_one = ifelse(under_one < 13/16, 0, under_one)),
              fill = "#08519c", alpha = 1/8) +
    geom_area(data = plot_dat %>% mutate(under_one = ifelse(under_one < 11/16, 0, under_one)),
              fill = "#08519c", alpha = 1/8) +
    geom_area(data = plot_dat %>% mutate(under_one = ifelse(under_one < 9/16, 0, under_one)),
              fill = "#08519c", alpha = 1/8) +
    geom_area(data = plot_dat %>% mutate(under_one = ifelse(under_one < 7/16, 0, under_one)),
              fill = "#08519c", alpha = 1/8) +
    geom_area(data = plot_dat %>% mutate(under_one = ifelse(under_one < 5/16, 0, under_one)),
              fill = "#08519c", alpha = 1/8) +
    geom_area(data = plot_dat %>% mutate(under_one = ifelse(under_one < 3/16, 0, under_one)),
              fill = "#08519c", alpha = 1/8) +
    geom_area(data = plot_dat %>% mutate(under_one = ifelse(under_one < 14/16, 0, under_one)),
              fill = "#08519c", alpha = 1/8) +
    geom_area(data = plot_dat %>% mutate(under_one = ifelse(under_one < 12/16, 0, under_one)),
              fill = "#08519c", alpha = 1/8) +
    geom_area(data = plot_dat %>% mutate(under_one = ifelse(under_one < 10/16, 0, under_one)),
              fill = "#08519c", alpha = 1/8) +
    geom_area(data = plot_dat %>% mutate(under_one = ifelse(under_one < 8/16, 0, under_one)),
              fill = "#08519c", alpha = 1/8) +
    geom_area(data = plot_dat %>% mutate(under_one = ifelse(under_one < 6/16, 0, under_one)),
              fill = "#08519c", alpha = 1/8) +
    geom_area(data = plot_dat %>% mutate(under_one = ifelse(under_one < 4/16, 0, under_one)),
              fill = "#08519c", alpha = 1/8) +
    geom_area(data = plot_dat %>% mutate(under_one = ifelse(under_one < 2/16, 0, under_one)),
              fill = "#08519c", alpha = 1/8) +
    geom_area(data = plot_dat %>% mutate(under_one = ifelse(under_one < 1/16, 0, under_one)),
              fill = "#08519c", alpha = 1/8) +
    # geom_col(aes(fill = -under_one), width = 1) +
    geom_vline(xintercept = ymd("2020-03-16"), lty = 2) +
    geom_vline(xintercept = ymd("2020-03-24"), lty = 2) +
    geom_vline(xintercept = ymd("2020-05-04"), lty = 2) +
    geom_vline(xintercept = ymd("2020-06-15"), lty = 2) +
    geom_vline(xintercept = ymd("2020-07-31"), lty = 2) +
    scale_fill_distiller() +
    scale_x_date(date_breaks = "month",
                 date_labels = "%B %d",
                 expand = expansion(add = 0),
                 limits = c(ymd("2020-02-27"), Sys.Date() + 1)) +
    scale_y_continuous(breaks = pretty_breaks(5), labels = label_percent(), 
                       limits = c(0, 1.01), expand = expansion(mult = 0)) +
    ggtitle(label = waiver(),
            subtitle = latex2exp::TeX("$P(R_t < 1)$")) +
    theme(axis.title = element_blank(),
          plot.margin = margin(5, 5, 5, 3))

plot_grid(p3, p4, p4b, ncol = 1)

spread_draws(m, R[day], phi) %>% 
    ungroup %>% 
    mutate(infected_by_one = rnbinom(n(), mu = 1 * R, size = phi)) %>% 
    group_by(day) %>% 
    summarise(lower_50 = quantile(infected_by_one, 0.25),
              upper_50 = quantile(infected_by_one, 0.75),
              lower_60 = quantile(infected_by_one, 0.2),
              upper_60 = quantile(infected_by_one, 0.8),
              lower_70 = quantile(infected_by_one, 0.15),
              upper_70 = quantile(infected_by_one, 0.85),
              lower_80 = quantile(infected_by_one, 0.1),
              upper_80 = quantile(infected_by_one, 0.9),
              lower_90 = quantile(infected_by_one, 0.05),
              upper_90 = quantile(infected_by_one, 0.95),
              lower_95 = quantile(infected_by_one, 0.025),
              upper_95 = quantile(infected_by_one, 0.975)) %>% 
    pivot_longer(c(-day), names_to = c("which", "prob"), names_sep = "_") %>% 
    pivot_wider(names_from = which, values_from = value) %>% 
    mutate(prob = parse_number(prob),
           date = ymd("2020-02-28") + day) %>% 
    ggplot(aes(date, ymin = lower, ymax = upper)) +
    geom_ribbon(aes(fill = factor(-prob)), alpha = 0.7) +
    geom_hline(yintercept = 1, lty = 2) +
    geom_vline(xintercept = ymd("2020-03-16"), lty = 2) +
    geom_vline(xintercept = ymd("2020-03-24"), lty = 2) +
    geom_vline(xintercept = ymd("2020-05-04"), lty = 2) +
    geom_vline(xintercept = ymd("2020-06-15"), lty = 2) +
    geom_vline(xintercept = ymd("2020-07-31"), lty = 2) +
    scale_fill_brewer() +
    scale_x_date(date_breaks = "month",
                 date_labels = "%B %d",
                 expand = expansion(add = 0),
                 limits = c(ymd("2020-02-27"), Sys.Date() + 1)) +
    scale_y_continuous(breaks = pretty_breaks(8), limits = c(0, NA), expand = expansion(mult = 0)) +
    ggtitle(label = waiver(),
            subtitle = latex2exp::TeX("$R_t$")) +
    theme(axis.title = element_blank(),
          plot.margin = margin(5, 5, 5, 14))


R_draws <- spread_draws(m, R[day]) %>% 
    group_by(day) %>% 
    mutate(iter = row_number()) %>%
    ungroup %>% 
    select(iter, day, R)

last_R <- R_draws %>% 
    filter(day == max(day)) %>% 
    .$R

pred_days <- 91
N_iter <- 2000
future_R <- crossing(day = max(R_draws$day) + seq_len(pred_days),
                     iter = seq_len(N_iter)) %>% 
    group_by(iter) %>% 
    mutate(R = last_R[iter] + 0.9 * scale(day - max(R_draws$day) - 1, center = F) - 0.5 * scale(day - max(R_draws$day), center = F)^2,
           R = as.numeric(R) + cumsum(rnorm(n(), mean = - 0.01 * (R - 0.8), sd = 0.01))) %>% 
    ungroup

future_R %>% 
    ggplot(aes(day, R, group = iter)) +
    geom_line(alpha = 0.1) +
    geom_hline(yintercept = 1, lty = 2)



R_draws %>% 
    bind_rows(future_R) %>% 
    ggplot(aes(day, R, group = iter)) +
    geom_line(alpha = 0.01) +
    geom_hline(yintercept = 1, lty = 2)

shape <- 1.54
rate <- 0.28
N_days <- nrow(d)
SI_dat <- tibble(t = seq(1, N_days)) %>% 
    mutate(
        p = case_when(
            TRUE ~ pgamma(t + 0.5, shape = shape, rate = rate) - pgamma(t - 0.5, shape = shape, rate = rate)
        ),
        p = p / sum(p)
    )

SI <- SI_dat$p


lambda <- numeric(N_days)

for (t in 2:N_days) {
    lambda[t] <- t(head(d$total, t - 1)) %*% tail(rev(SI), t - 1) / sum(tail(rev(SI), t - 1))
}

lambda <- c(lambda, rep(0, pred_days))



make_preds <- function(d, ...) {
    for (t in seq(N_days, nrow(d))) {
        d$lambda[t] <- t(d$mu_hat[(t - 1):(t - length(SI) + 1)]) %*% head(SI, -1)
        d$mu_hat[t] <- d$lambda[t] * d$R[t]
    }
    
    d
}



plot_dat <- R_draws %>% 
    bind_rows(future_R) %>% 
    group_by(iter) %>% 
    mutate(lambda = lambda[-1],
           mu_hat = R * lambda) %>% 
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

p5 <- plot_dat %>% 
    filter(name == "y_hat") %>% 
    ggplot(aes(date, ymin = lower, ymax = upper)) +
    geom_ribbon(aes(fill = factor(-prob)), alpha = 0.7) +
    geom_point(data = d %>% rename(y_hat = local) %>% pivot_longer(c(y_hat)),
               inherit.aes = F, aes(x = date, y = value)) +
    geom_vline(xintercept = Sys.Date(), lty = 2) +
    scale_x_date(date_breaks = "month", 
                 date_labels = "%B %d",
                 limits = c(ymd("2020-02-27"), Sys.Date() + 1 + pred_days), 
                 expand = expansion(add = 0)) +
    scale_y_continuous(expand = expansion(mult = 0.01)) +
    scale_fill_brewer() +
    labs(subtitle = "New local cases") +
    theme(axis.title = element_blank()) +
    theme(axis.title = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          plot.margin = margin(5, 5, 5, 5),
          legend.title = element_blank()) 

p6 <- plot_dat %>% 
    filter(name == "R") %>% 
    ggplot(aes(date, ymin = lower, ymax = upper)) +
    geom_ribbon(aes(fill = factor(-prob)), alpha = 0.7) +
    geom_hline(yintercept = 1, lty = 2) +
    geom_vline(xintercept = Sys.Date(), lty = 2) +
    scale_x_date(date_breaks = "month", date_labels = "%B %d",
                 limits = c(ymd("2020-02-27"), Sys.Date() + 1 + pred_days), 
                 expand = expansion(add = 0)) +
    scale_y_continuous(expand = expansion(mult = 0.01), breaks = pretty_breaks(8)) +
    scale_fill_brewer() +
    ggtitle(label = waiver(),
            subtitle = latex2exp::TeX("$R_t$")) +
    theme(axis.title = element_blank(),
          plot.margin = margin(5, 5, 5, 8))

plot_grid(p5, p6, ncol = 1)





