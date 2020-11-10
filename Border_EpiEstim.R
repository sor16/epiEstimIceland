source("Scenario_EpiEstim.R")

## Estimate for R

fill_up_dates <- ymd(c('2020-12-01', '2020-11-01'))

R_estim_mod <- read_rds(here("Results", "EpiEstim", str_c("R.estimate.pred_", fit_date,".rds")))
R_estim_m <- R_estim_mod$draws() %>% as_draws_df

R_estim_p <- spread_draws(R_estim_m, R[day]) %>% 
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
           date = Sys.Date() + day) %>% 
    ggplot() +
    geom_ribbon(aes(date, ymin = lower, ymax = upper,fill = factor(-prob)), alpha = 0.7) +
    geom_hline(yintercept = 1,col='#BC3C29FF') +
    geom_vline(xintercept = as.numeric(Sys.Date()), lty = 2) + 
    scale_fill_brewer() +
    scale_x_date(breaks = sort(c(fill_up_dates, head(intervention_dat$date,-1), Sys.Date())),
                 labels = icelandic_dates,
                 expand = expansion(add = 0),
                 limits = c(Sys.Date()+1, ymd(end_date)+60)) +
    scale_y_continuous(breaks = pretty_breaks(8),expand = c(0,0)) +
    ggtitle(label = waiver(),
            subtitle = latex2exp::TeX("Þróun smitstuðulsins ($R_t$) innanlands utan sóttkvíar")) +
    theme(axis.title = element_blank(),
          plot.margin = margin(5, 5, 5, 14))
R_estim_p
ggsave(here('Results','Figures', 'Predictions', paste0('R_estimate_',Sys.Date(),'.png')),height=4.5,width=10,device='png')


## Estimate for new cases

Ro = 0.5
Rq = 0.5

R_draws <- spread_draws(m, R[day]) %>% 
    group_by(day) %>% 
    mutate(iter = row_number()) %>%
    ungroup() %>% 
    select(iter, day, R) %>%
    mutate(R = Ro, R_q=Rq)

posterior_dat <- R_draws %>% filter(day <= 60) %>%
    group_by(iter) %>% 
    mutate(prop_quarantine=c(d$prop_quarantine),
           lambda = c(lambda_vec, rep(0,pred_days-1)),
           lambda_q = c(lambda_q_vec, rep(0,pred_days-1)),
           pred_day = if_else(day>N_days, 1,0),
           mu_hat = R * lambda + R_q*lambda_q + d$border_1) %>% #(day>N_days)*as.numeric(rbinom(n(),size=num_border_tests,prob = prop_imported))) %>% # 
    # filter(iter %in% 1:100) %>% 
    group_modify(make_preds, N_days=N_days,SI=SI) %>% 
    ungroup %>% 
    mutate(y_hat = rnbinom(n(), mu = mu_hat, size = m$phi[iter])) %>% 
    pivot_longer(c(-iter, -day, -lambda, -mu_hat)) %>% summarize_posterior_dat

pred_lp <- plot_dat %>% filter(name=='y_hat') %>%
    ggplot(aes(date, ymin = lower, ymax = upper)) +
    geom_ribbon(aes(fill = factor(-prob)), alpha = 0.7) +
    geom_point(data = d, aes(x = date, y = total), inherit.aes = F) +
    geom_vline(xintercept = as.numeric(intervention_dat$date), lty = 2) +     
    scale_fill_brewer() +
    scale_x_date(limits = c(ymd("2020-02-28"), ymd('2020-02-28') + 60), 
                 expand = expansion(add = 0),
                 breaks = sort(c(fill_up_dates, intervention_dat$date)),
                 labels = icelandic_dates) +
    scale_y_continuous(breaks = pretty_breaks(8), 
                       expand = expansion(mult = 0)) +
    labs(subtitle = "Innlend dagleg smit") +
    theme(axis.title = element_blank(),
          plot.margin = margin(5, 5, 5, 11),
          legend.title = element_blank()) 
pred_lp