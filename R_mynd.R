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
library(zoo)

theme_set(theme_classic(base_size = 12) + 
              theme(legend.position = "none"))
icelandic_dates <- function(x) {
    months <- c("janúar", "febrúar", "mars", "apríl", "maí", "júní", 
                "júlí", "ágúst", "september", "október", "nóvember", "desember")
    
    paste0(mday(x), ". ", months[month(x)])
}
fit_date <- Sys.Date()

fill_up_dates <- ymd(c('2020-03-01','2020-04-01','2020-06-01','2020-07-01'))

height=6.7
diff=0.5
extra_diff=2
start_date = '2020-02-28'
end_date = Sys.Date()

intervention_dat <- tibble(date=ymd(c("2020-03-16","2020-03-24","2020-05-04","2020-06-15","2020-07-31","2020-08-19",'2020-09-07', '2020-09-18')),
                           lab=c('Samkomur,takmarkaðar,við 100 manns','Samkomur,takmarkaðar,við 20 manns','Fjöldatakmarkanir,rýmkaðar'
                                 ,'Landamæraskimun,hefst','Samkomur,takmarkaðar,við 100 manns','Sóttkví milli,tveggja skimanna,á landamærunum'
                                 ,'Eins metra regla,í stað tveggja,metra reglu', 'Skemmtistöðum,og krám á,höfuðborgarsvæðinu,lokað tímabundið'),
                           lab_pos_x=date,
                           lab_pos_y=height)
intervention_dat$lab_pos_x[1] <- ymd('2020-03-01')
intervention_dat$lab_pos_x[nrow(intervention_dat)-1] <- ymd('2020-08-22')
intervention_dat$lab_pos_y[nrow(intervention_dat)-1] <- height+extra_diff
intervention_dat$lab_pos_x[nrow(intervention_dat)] <- ymd('2020-09-18')+1
intervention_dat$lab_pos_y[nrow(intervention_dat)] <- height+extra_diff

intervention_lab_dat <- separate_rows(intervention_dat,lab,sep=',') %>% 
                        group_by(date) %>% 
                        mutate(lab_pos_y=seq(lab_pos_y[1],lab_pos_y[1]-diff*(n()-1),by=-diff))

mod <- read_rds(here("Results", "EpiEstim", str_c("tot.Q.late.border_", fit_date,".rds")))
m <- mod$draws() %>% as_draws_df
rm(mod)

breaks_type='month'

p <- spread_draws(m, R[day]) %>% 
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
           date = ymd(start_date) + day) %>% 
    ggplot() +
    geom_ribbon(aes(date, ymin = lower, ymax = upper,fill = factor(-prob)), alpha = 0.7) +
    geom_hline(yintercept = 1,col='#BC3C29FF') +
    geom_label(data = intervention_lab_dat, aes(lab_pos_x, lab_pos_y,label=lab), hjust = 0, label.size = NA, fill = 'white') +
    geom_vline(xintercept = as.numeric(intervention_dat$date), lty = 2) + 
    scale_fill_brewer() +
    scale_x_date(breaks = sort(c(fill_up_dates, intervention_dat$date)),
                 labels = icelandic_dates,
                 expand = expansion(add = 0),
                 limits = c(ymd(start_date) - 1, ymd(end_date) + 1)) +
    scale_y_continuous(breaks = pretty_breaks(8),expand = c(0,0)) +
    ggtitle(label = waiver(),
            subtitle = latex2exp::TeX("Þróun smitstuðulsins ($R_t$) innanlands utan sóttkvíar")) +
    theme(axis.title = element_blank(),
          plot.margin = margin(5, 5, 5, 14))

ggsave(paste0('~/epiEstimIceland/Results/Figures/','Rt_',Sys.Date(),'.pdf'),height=4.5,width=19,device='pdf')
ggsave(paste0('~/epiEstimIceland/Results/Figures/','Rt_',Sys.Date(),'.png'),height=4.5,width=19,device='png')

d <- read_csv("https://docs.google.com/spreadsheets/d/1xgDhtejTtcyy6EN5dbDp5W3TeJhKFRRgm6Xk0s0YFeA/export?format=csv&gid=1788393542",col_types=cols()) %>%
    select(date = Dagsetning, local = Innanlands_Smit,border_1=Landamaeri_Smit_1,border_2=Landamaeri_Smit_2, imported = Innflutt_Smit,prop_quarantine=Hlutf_Sottkvi,num_quarantine=Fjoldi_Sottkvi) %>% 
    mutate(date = ymd(date),
           total = if_else(date >= ymd("2020-07-23"),local,local + imported),
           prop_quarantine=if_else(total!=0,num_quarantine/total,0)) %>% 
    filter(date >= ymd(start_date) & date <= ymd(end_date))

lp<- d %>%
    mutate(rmean = c(rep(0,6),rollmean(local, 7))) %>% 
    ggplot(aes(date, local)) +
    geom_point(
        data = d, aes(x=date, y=local), inherit.aes = F
    ) +
    geom_line(aes(y=rmean),color = '#0072B5FF') +
    geom_vline(xintercept = as.numeric(intervention_dat$date), lty = 2) + 
    scale_x_date(breaks = sort(c(fill_up_dates, intervention_dat$date)),
                 labels = icelandic_dates,
                 expand = expansion(add = 0),
                 limits = c(ymd(start_date) - 1, ymd(end_date) + 1)) +
    scale_y_continuous(breaks = pretty_breaks(8), expand = expansion(mult = 0), limits = c(0, 145)) +
    labs(subtitle = "Dagleg greind smit innanlands") +
    theme(axis.title = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          plot.margin = margin(5, 5, 5, 22),
          legend.title = element_blank())

plot_grid(p, lp, ncol = 1, align='v') +
    ggsave(here("Results", "Figures", paste0("R_and_local_cases_",Sys.Date(),".png")), device = "png", 
           height = 4.5 * 2, width = 19)
