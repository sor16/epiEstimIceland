library(cmdstanr)
library(rstan)
library(posterior)
library(broom.mixed)
library(tidyverse)
library(cowplot)
library(scales)
library(lubridate)
library(cmdstanr)
library(posterior)
library(tidybayes)
library(here)
library(zoo)
library(googlesheets4)
library(gargle)

source("Scenario_EpiEstim.R")

theme_set(theme_classic(base_size = 12) + theme(legend.position = "none"))

make_preds <- function(draws, SI,pos=1,start_day,...) {
    pop = 356991
    for (t in start_day:nrow(draws)) {
        idx_range <- (t-1):1
        draws$lambda[t] <- t((1-draws$prop_quarantine[idx_range])*draws$mu_hat[idx_range]) %*% head(SI,t-1)
        draws$lambda_q[t] <- t(draws$prop_quarantine[idx_range]*draws$mu_hat[idx_range]) %*% head(SI,t-1)
        #draws$lambda_i[t] <- t(draws$imported[idx_range]) %*% head(SI,t-1)
        draws$lambda_i[t] <- t(draws$imported[(t-1):1]) %*% SI[pos:(t-1+pos-1)]
        
        draws$mu_hat[t] <- (1-sum(draws$mu_hat)/pop)*(draws$lambda[t] * draws$R[t] + draws$lambda_q[t] * draws$R_q[t] + draws$lambda_i[t] * draws$R_i[t]) #+ draws$mu_hat[t]## add artificially imported infections already stored in mu_hat[t]
    }
    return(draws)
}

make_preds_smitgat <- function(draws, SI,start_day,...) {
    pop = 356991
    for (t in start_day:nrow(draws)) {
        idx_range <- (t-1):1
        draws$lambda[t] <- t((1-draws$prop_quarantine[idx_range])*draws$mu_hat[idx_range]) %*% head(SI,t-1)
        draws$lambda_q[t] <- t(draws$prop_quarantine[idx_range]*draws$mu_hat[idx_range]) %*% head(SI,t-1)

        draws$lambda_i[t] <- t(draws$infected_first_test[max(1,(t-5)):(t-1)]) %*% head(SI,min(5,t-1))
        if(t > 6)
            draws$lambda_i2[t] <- t(draws$imported[(t-1):6]) %*% tail(head(SI,t-1),t-6)
        
        draws$mu_hat[t] <- (1-sum(draws$mu_hat)/pop)*(draws$lambda[t] * draws$R[t] + draws$lambda_q[t] * draws$R_q[t] + draws$lambda_i[t] * draws$R_i[t]) + draws$lambda_i2[t] * draws$R_i2[t] #+ draws$mu_hat[t]## add artificially imported infections already stored in mu_hat[t]
    }
    return(draws)
}

icelandic_dates <- function(x) {
    months <- c("janúar", "febrúar", "mars", "apríl", "maí", "júní", 
                "júlí", "ágúst", "september", "október", "nóvember", "desember")
    
    paste0(mday(x), ". ", months[month(x)])
}

## Fetch dispersion data from today's R model
fit_date <- Sys.Date()
mod <- read_rds(here("Results", "EpiEstim", str_c("tot.Q.late.border_", '2020-12-09',".rds")))
m <- mod$draws() %>% as_draws_df

# Set fixed parameters
N_iter <- 2000
#N_days <- 30
pred_days <- 60
mu <- 0
prop_q = 0
R_i1 = 1
R_i2 = 0.5

scenarios <- crossing(R = c(2), pred_days = 60, prevalence=c("05","1","2"), type=c('PCRv', 'smitgat', "one_test", "two_test_q"), tourists=c(200,500,1000,1500))

future_prop_quarantine = rep(prop_q,pred_days)

delay = 7
for(j in 1:(pred_days-delay)) {
    future_prop_quarantine[j+delay] = 0.5*(2/(1+exp(-j*12/42))-1)
}

tmp = tibble(ind = 1:60, pq = future_prop_quarantine) %>% mutate(date = ymd('2020-01-31')+ind)
ggplot(tmp, aes(date,pq)) + geom_line() + coord_cartesian(ylim=c(0,1)) +
    scale_x_date(breaks=ymd(c('2020-02-01', '2020-03-01')), labels = icelandic_dates)+
    scale_y_continuous(expand=c(0,0), breaks=c(0,0.5,1), labels=c('0', '0.5', '1')) + ylab("") + xlab("") +
    ggsave(here('Results', 'Figures', 'Other', 'sigmoid.png'),height=3,width=5,device='png')

SI <- get_SI_vec(100)

# Generate infections for N_days days and then estimate imported cases added on for pred_days days
pre_scenario_dat <- crossing(day = 1:30,iter = seq_len(N_iter)) %>% 
                    group_by(iter) %>% 
                    mutate(R = 2,
                           R_q = m$R_q[iter],
                           R_i = R_i1,
                           R_i2 = R_i2,
                           prop_quarantine=prop_q,
                           prop_imported = 0,
                           total = rnbinom(n(), mu = mu, size = m$phi[iter]),
                           lambda = calculate_lambda(total, SI, 1-prop_quarantine),
                           lambda_q = calculate_lambda(total, SI, prop_quarantine),
                           lambda_i = 0,
                           lambda_i2 = 0,
                           imported = 0,
                           infected_first_test = 0,
                           mu_hat = total)

for(i in 1:nrow(scenarios)) {
    s = scenarios[i,]
    imported_infections <- read_csv(here("Simulation_Data",paste0("01_", s$tourists, if_else(s$type == "PCRv", "_PCRv_", "_smitgat_"), s$prevalence, ".csv")), col_types = cols()) %>%
                           mutate(day = as.numeric(difftime(date, ymd('2021-01-01'), units = "days"))) %>%
                           filter(day <= s$pred_days+30) %>%
                           #select(sleppa, day) %>%
                           rename(imported = sleppa) %>%
                           mutate(iter = rep(1:(n()/60), each=s$pred_days)) %>%
                           filter(iter <= N_iter)
    
    future_dat <- crossing(day = seq_len(s$pred_days)+30,iter = seq_len(N_iter)) %>% 
        group_by(iter) %>% 
        mutate(R = s$R,
               R_q = m$R_q[iter],
               R_i = R_i1,
               R_i2 = R_i2,
               prop_quarantine=future_prop_quarantine,
               prop_imported = 0,
               total = 0,
               lambda = 0,
               lambda_q = 0,
               lambda_i = 0,
               lambda_i2 = 0,
               mu_hat = 0) %>%
        ungroup() %>%
        mutate(imported = if(s$type == "one_test") c(imported_infections$test2_smitadir[6:nrow(imported_infections)], rep(0,5)) else imported_infections$imported, 
               infected_first_test = c(imported_infections$test2_smitadir[6:nrow(imported_infections)], rep(0,5)))
    
    scenario_dat <- pre_scenario_dat %>% bind_rows(future_dat)
    
    if(s$type == 'smitgat') {
        posterior_dat <- scenario_dat %>%
            group_by(iter) %>%
            group_modify(make_preds_smitgat,SI=SI, start_day=31) %>%
            ungroup %>%
            mutate(y_hat = rnbinom(n(), mu = mu_hat, size = m$phi[iter]))
        save(posterior_dat, file=here('shiny2' ,'data', paste0("border_", s$type, "_prevalence_", s$prevalence,"_tourists_",s$tourists, ".Rdata")))
        print(paste('Scenario', i))
    } else {
        posterior_dat <- scenario_dat %>%
            group_by(iter) %>%
            group_modify(make_preds,SI=SI, pos=if_else(s$type=="one_test", 1, if_else(s$type=="two_test_q", 6, 4)), start_day=31) %>%
            ungroup %>%
            group_by(day, iter) %>% 
            mutate(y_hat = if(day>=31) rnbinom(n(), mu = mu_hat, size = m$phi[iter]) else total)
        save(posterior_dat, file=here('shiny2' ,'data', paste0("border_", s$type, "_prevalence_", s$prevalence,"_tourists_",s$tourists, ".Rdata")))
        print(paste('Scenario', i))
    }
}

## Create images for text

theme_set(theme_classic(base_size = 12))# + 
#theme(legend.position = "none"))

plot_results <- function(type, prevalence, tourists,legend=F) {
    load(here('shiny2', 'data', paste0("border_", type, "_prevalence_", prevalence,"_tourists_", tourists, ".Rdata")))
    plot_dat <- posterior_dat %>% group_by(day) %>% 
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
        mutate(prob = parse_number(prob), date = ymd("2020-03-01")+day-31)
    
    med = posterior_dat %>% group_by(day) %>% summarise(median = median(y_hat)) %>% select(median)
    med = med$median
    start_date = Sys.Date()
    pop = 356991
    red_day = 100
    for(t in 45:90) {
        sum = 0
        for(i in 1:14) {
            sum = sum + med[t-i]
        }
        if(sum/(pop/100000) >= 150) {
            red_day = t
            break
        }
    }
    
    over_20_day = 100
    for(t in 30:90) {
        if(med[t] >= 20) {
            over_20_day = t
            break
        }
    }
    
    type_label =""
    if(type == 'PCRv')  type_label = "tvær skimanir, PCR vottorð"
    else if (type == 'smitgat') type_label = "tvær skimanir og ferðamannasmitgát"
    else if (type == 'one_test') type_label = "ein skimun"
    else if (type == 'two_test_q') type_label = "tvær skimanir, sóttkví á milli skimana"
    
    fills <- tibble(color = c("red", "orange"), xint = c(ymd("2020-02-01")-30+red_day, ymd("2020-02-01")-30+over_20_day), text=c("Ísland á rauðum lista ECDC", "Yfir 20 dagleg smit (miðgildi)"))
    plot_dat %>% ggplot(aes(date, ymin = lower, ymax = upper)) +
        geom_ribbon(aes(fill = factor(-prob)), alpha = 0.7, show.legend = F) +
        scale_fill_brewer() +
        geom_vline(aes(color="Ísland á rauðum lista ECDC", xintercept=ymd("2020-03-01")-30+red_day), lty=2)+
        geom_vline(aes(color="Yfir 20 dagleg smit (miðgildi)", xintercept=ymd("2020-03-01")-30+over_20_day), lty=2)+
        scale_color_manual(values=fills$color) + 
        labs(subtitle = paste0("Daglegur fjöldi innanlandssmita miðað við ",tourists," ferðamenn,", if_else(prevalence == "05", "0.5", prevalence),"% algengi, R=2 og ", type_label)) +
        scale_x_date(breaks=c(ymd(c('2020-03-01', if_else((red_day < 85 || red_day > 90) && (over_20_day < 85 || over_20_day > 90),'2020-05-01', '')), ymd("2020-03-01")-30+red_day, ymd("2020-03-01")-30+over_20_day)), expand=c(0,0), labels=icelandic_dates, limits = c(ymd("2020-02-28"), ymd("2020-05-03"))) + 
        scale_y_continuous(breaks=pretty_breaks(8), expand=c(0,0)) +
        coord_cartesian(ylim=c(0,150)) + 
        theme(axis.title = element_blank(),
              plot.margin = margin(5, 5, 5, 11),
              legend.title = element_blank(), legend.position=if_else(legend==T, "bottom", "none"), legend.direction="horizontal") 
}
    
plt <- plot_grid(plot_results("one_test", "05", 1000), 
          plot_results("smitgat", "05", 1000),
          plot_results("PCRv", "05", 1000),
          plot_results("two_test_q", "05", 1000),
          plot_results("one_test", "1", 1000), 
          plot_results("smitgat", "1", 1000),
          plot_results("PCRv", "1", 1000),
          plot_results("two_test_q", "1", 1000),
          plot_results("one_test", "2", 1000), 
          plot_results("smitgat", "2", 1000),
          plot_results("PCRv", "2", 1000),
          plot_results("two_test_q", "2", 1000),
          ncol = 4, nrow=3, align='h')  

plot_results("one_test", "05", 1000, legend = T)

title <- ggdraw() + draw_label("Spá yfir 60 daga tímabil miðað við 1000 ferðamenn á dag", fontface='bold')
plot_grid(title, plt, ncol=1, rel_heights=c(0.1, 1)) +
    ggsave(here("Results", "Figures","Other", paste0("Border_scenarios.png")), device = "png", height = 12, width = 22)  


for(i in 1:nrow(scenarios)) {
    plot_results(scenarios[i,]$type, scenarios[i,]$prevalence, scenarios[i,]$tourists, legend=T) + 
        ggsave(here("Results", "Figures", "Other", paste0("COVID_19_svidsmynd_",case_when(scenarios[i,]$type=="one_test" ~ "ein_skimun_", scenarios[i,]$type=="PCRv" ~ "PCRv_", scenarios[i,]$type=="two_test_q" ~ "tvaer_skimanir_sottkvi_", scenarios[i,]$type=="smitgat" ~ "smitgat_"), scenarios[i,]$prevalence, "_algengi_", scenarios[i,]$tourists, "_ferdamenn.png")),
               device = "png", height = 8, width = 16)
}

 # p <- posterior_dat %>% 
#     group_by(iter) %>%
#     mutate(cumc = cumsum(y_hat)) %>%
#     ungroup %>%
#     group_by(day) %>% 
#     summarise(lower_50 = quantile(cumc, 0.25),
#               upper_50 = quantile(cumc, 0.75),
#               lower_60 = quantile(cumc, 0.2),
#               upper_60 = quantile(cumc, 0.8),
#               lower_70 = quantile(cumc, 0.15),
#               upper_70 = quantile(cumc, 0.85),
#               lower_80 = quantile(cumc, 0.1),
#               upper_80 = quantile(cumc, 0.9),
#               lower_90 = quantile(cumc, 0.05),
#               upper_90 = quantile(cumc, 0.95),
#               lower_95 = quantile(cumc, 0.025),
#               upper_95 = quantile(cumc, 0.975)) %>% 
#     pivot_longer(c(-day), names_to = c("which", "prob"), names_sep = "_") %>% 
#     pivot_wider(names_from = which, values_from = value) %>% 
#     mutate(prob = parse_number(prob)) %>% 
#     ggplot(aes(day, ymin=lower, ymax=upper)) +
#     geom_ribbon(aes(fill = factor(-prob)), alpha = 0.7) +
#     scale_fill_brewer() +
#     labs(subtitle = paste0("Smit yfir 60 daga tímabil miðað við ", s$algengi,"% algengi, R=",s$R, " og ", test_label, " landamæraskimun")) +
#     scale_x_continuous(breaks=pretty_breaks(6), expand=c(0,0), limits=c(0,14)) + 
#     scale_y_continuous(breaks=pretty_breaks(8), expand=c(0,0)) +
#     coord_cartesian(ylim=c(0,if_else(s$tests == 0, 50, 50))) + 
#     theme(axis.title = element_blank(),
#           plot.margin = margin(5, 5, 5, 11),
#           legend.title = element_blank())
# p
# ## Estimate for new cases
# 
# posterior_dat <- R_draws %>% filter(day <= 60) %>%
#     group_by(iter) %>% 
#     mutate(prop_quarantine=c(d$prop_quarantine),
#            lambda = c(lambda_vec, rep(0,pred_days-1)),
#            lambda_q = c(lambda_q_vec, rep(0,pred_days-1)),
#            pred_day = if_else(day>N_days, 1,0),
#            mu_hat = R * lambda + R_q*lambda_q + d$border_1) %>% #(day>N_days)*as.numeric(rbinom(n(),size=num_border_tests,prob = prop_imported))) %>% # 
#     # filter(iter %in% 1:100) %>% 
#     group_modify(make_preds, N_days=N_days,SI=SI) %>% 
#     ungroup %>% 
#     mutate(y_hat = rnbinom(n(), mu = mu_hat, size = m$phi[iter])) %>% 
#     pivot_longer(c(-iter, -day, -lambda, -mu_hat)) %>% summarize_posterior_dat
# 
# pred_lp <- plot_dat %>% filter(name=='y_hat') %>%
#     ggplot(aes(date, ymin = lower, ymax = upper)) +
#     geom_ribbon(aes(fill = factor(-prob)), alpha = 0.7) +
#     geom_point(data = d, aes(x = date, y = total), inherit.aes = F) +
#     geom_vline(xintercept = as.numeric(intervention_dat$date), lty = 2) +     
#     scale_fill_brewer() +
#     scale_x_date(limits = c(ymd("2020-02-28"), ymd('2020-02-28') + 60), 
#                  expand = expansion(add = 0),
#                  breaks = sort(c(fill_up_dates, intervention_dat$date)),
#                  labels = icelandic_dates) +
#     scale_y_continuous(breaks = pretty_breaks(8), 
#                        expand = expansion(mult = 0)) +
#     labs(subtitle = "Innlend dagleg smit") +
#     theme(axis.title = element_blank(),
#           plot.margin = margin(5, 5, 5, 11),
#           legend.title = element_blank()) 
# pred_lp
# 
# ## Estimate for R
# 
# fill_up_dates <- ymd(c('2020-12-01', '2020-11-01'))
# 
# R_estim_mod <- read_rds(here("Results", "EpiEstim", str_c("R.estimate.pred_", fit_date,".rds")))
# R_estim_m <- R_estim_mod$draws() %>% as_draws_df
# 
# R_estim_p <- spread_draws(R_estim_m, R[day]) %>% 
#     group_by(day) %>%
#     summarise(lower_50 = quantile(R, 0.25),
#               upper_50 = quantile(R, 0.75),
#               lower_60 = quantile(R, 0.2),
#               upper_60 = quantile(R, 0.8),
#               lower_70 = quantile(R, 0.15),
#               upper_70 = quantile(R, 0.85),
#               lower_80 = quantile(R, 0.1),
#               upper_80 = quantile(R, 0.9),
#               lower_90 = quantile(R, 0.05),
#               upper_90 = quantile(R, 0.95),
#               lower_95 = quantile(R, 0.025),
#               upper_95 = quantile(R, 0.975)) %>% 
#     pivot_longer(c(-day), names_to = c("which", "prob"), names_sep = "_") %>% 
#     pivot_wider(names_from = which, values_from = value) %>% 
#     mutate(prob = parse_number(prob),
#            date = Sys.Date() + day) %>% 
#     ggplot() +
#     geom_ribbon(aes(date, ymin = lower, ymax = upper,fill = factor(-prob)), alpha = 0.7) +
#     geom_hline(yintercept = 1,col='#BC3C29FF') +
#     geom_vline(xintercept = as.numeric(Sys.Date()), lty = 2) + 
#     scale_fill_brewer() +
#     scale_x_date(breaks = sort(c(fill_up_dates, head(intervention_dat$date,-1), Sys.Date())),
#                  labels = icelandic_dates,
#                  expand = expansion(add = 0),
#                  limits = c(Sys.Date()+1, ymd(end_date)+60)) +
#     scale_y_continuous(breaks = pretty_breaks(8),expand = c(0,0)) +
#     ggtitle(label = waiver(),
#             subtitle = latex2exp::TeX("Þróun smitstuðulsins ($R_t$) innanlands utan sóttkvíar")) +
#     theme(axis.title = element_blank(),
#           plot.margin = margin(5, 5, 5, 14))
# R_estim_p
# ggsave(here('Results','Figures', 'Predictions', paste0('R_estimate_',Sys.Date(),'.png')),height=4.5,width=10,device='png')