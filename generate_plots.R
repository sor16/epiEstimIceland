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
icelandic_dates <- function(x) {
    months <- c("janúar", "febrúar", "mars", "apríl", "maí", "júní", 
                "júlí", "ágúst", "september", "október", "nóvember", "desember")
    
    paste0(mday(x), ". ", months[month(x)])
}
fit_date <- Sys.Date()
intervention_dates <- ymd(c("2020-03-16","2020-03-24","2020-05-04","2020-06-15","2020-07-31","2020-08-19"))

model_table <- tibble(model_name=c('tot','2.wave_plus.3.wave','1.wave','2.wave','3.wave'),
                      start_date=c('2020-02-28','2020-07-23','2020-02-28','2020-07-23','2020-09-11'),
                      end_date=c(Sys.Date(),Sys.Date(),'2020-05-04','2020-09-01',Sys.Date()))
Q_model_table <- tibble(model_name=c('tot.Q','tot.Q.late.border'),
                        start_date=c('2020-02-28','2020-02-28'),
                        end_date=c(Sys.Date(),Sys.Date()))
model_table <- bind_rows(model_table,Q_model_table)

path <- here("Results", "EpiEstim")
models <- list.files(path)
models <- models[grepl(fit_date,models)]
for(file in models){
    mod <- read_rds(paste0(path,'/',file))
    m <- mod$draws() %>% as_draws_df
    rm(mod)
    current_model_table <- filter(model_table,model_name==gsub('_.*','',file))
    current_intervention_dates <- intervention_dates[intervention_dates>current_model_table$start_date & intervention_dates<current_model_table$end_date]
    if(ymd(current_model_table$end_date)-ymd(current_model_table$start_date) < 60){
        breaks_type='week'
    }else{
        breaks_type='month'
    }
    
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
               date = ymd(current_model_table$start_date) + day) %>% 
        ggplot(aes(date, ymin = lower, ymax = upper)) +
        geom_ribbon(aes(fill = factor(-prob)), alpha = 0.7) +
        geom_hline(yintercept = 1, lty = 2)
        for(i_date in current_intervention_dates){
            p <- p+geom_vline(xintercept = i_date, lty = 2)
        }
        p <- p + 
            scale_fill_brewer() +
            scale_x_date(date_breaks = breaks_type,
                         labels = icelandic_dates,
                         expand = expansion(add = 0),
                         limits = c(ymd(current_model_table$start_date)-1, ymd(current_model_table$end_date)+1)) +
            scale_y_continuous(breaks = pretty_breaks(8),expand = c(0,0)) +
            ggtitle(label = waiver(),
                    subtitle = latex2exp::TeX("$R_t$")) +
            theme(axis.title = element_blank(),
                  plot.margin = margin(5, 5, 5, 14))
        ggsave(here('Results','Figures',str_c(current_model_table$model_name,'_',fit_date,'.png')),width=12,height=6,device='png')
}

### fyrir covid.hi.is
#file='tot.Q.late.border_2020-09-24.rds'
# name='tot.Q.late.border_2020-09-24'
# height=6
# diff=0.5
# extra_diff=2
name='tot_2020-09-24'
height=3.5
diff=0.25
extra_diff=1

mod <- read_rds(paste0(path,'/',name,'.rds'))
m <- mod$draws() %>% as_draws_df
rm(mod)
current_model_table <- filter(model_table,model_name==gsub('_.*','',file))
fill_up_dates <- ymd(c('2020-03-01','2020-04-01','2020-06-01','2020-07-01'))
intervention_dat <- tibble(date=ymd(c("2020-03-16","2020-03-24","2020-05-04","2020-06-15","2020-07-31","2020-08-19",'2020-09-07')),
                           lab=c('Samkomur,takmarkaðar,við 100 manns','Samkomur,takmarkaðar,við 20 manns','Fjöldatakmarkanir,rýmkaðar','Landamæraskimun,hefst','Samkomur,takmarkaðar,við 100 manns','Sóttkví milli,tveggja skimanna,á landamærunum','Eins metra regla,í stað tveggja, metra reglu'),
                           lab_pos_x=date,
                           lab_pos_y=height)
intervention_dat$lab_pos_x[1] <- ymd('2020-03-01')
intervention_dat$lab_pos_x[nrow(intervention_dat)] <- ymd('2020-08-22')
intervention_dat$lab_pos_y[nrow(intervention_dat)] <- height+extra_diff
intervention_lab_dat <- separate_rows(intervention_dat,lab,sep=',') %>% group_by(date) %>% mutate(lab_pos_y=seq(lab_pos_y[1],lab_pos_y[1]-diff*(n()-1),by=-diff))
if(ymd(current_model_table$end_date)-ymd(current_model_table$start_date) < 60){
    breaks_type='week'
}else{
    breaks_type='month'
}

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
           date = ymd(current_model_table$start_date) + day) %>% 
    ggplot() +
    geom_ribbon(aes(date, ymin = lower, ymax = upper,fill = factor(-prob)), alpha = 0.7) +
    geom_hline(yintercept = 1,col='#BC3C29FF') +
    geom_label(data=intervention_lab_dat,aes(lab_pos_x,lab_pos_y,label=lab),hjust=0,label.size=NA,fill='white')
for(i_date in intervention_dat$date){
    p <- p+geom_vline(xintercept = i_date, lty = 2)
}
p <- p + 
    scale_fill_brewer() +
    scale_x_date(breaks = sort(c(fill_up_dates,intervention_dat$date)),
                 labels = icelandic_dates,
                 expand = expansion(add = 0),
                 limits = c(ymd(current_model_table$start_date)-1, ymd(current_model_table$end_date)+1)) +
    scale_y_continuous(breaks = pretty_breaks(8),expand = c(0,0)) +
    ggtitle(label = waiver(),
            subtitle = latex2exp::TeX("Þróun smitstuðulsins ($R_t$) innanlands utan sóttkvíar")) +
    theme(axis.title = element_blank(),
          plot.margin = margin(5, 5, 5, 14))

ggsave(paste0('~/Documents/COVID-19/epiEstimIceland/Results/Figures/',name,'.pdf'),height=4.5,width=16,device='pdf')

