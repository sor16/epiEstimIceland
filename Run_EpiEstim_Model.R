library(tidyverse)
library(rstan)
library(lubridate)
library(readxl)
library(here)
library(cmdstanr)
library(posterior)

options(mc.cores = parallel::detectCores())
parallel:::setDefaultClusterOptions(setup_strategy = "sequential")
model_table <- tibble(model_name=c('tot','2_wave_plus_3_wave','1_wave','2_wave','3_wave'),
                 start_date=c('2020-02-28','2020-07-23','2020-02-28','2020-07-23','2020-09-11'),
                 end_date=c(Sys.Date(),Sys.Date(),'2020-05-04','2020-09-01',Sys.Date()))
Q_model_table <- tibble(model_name=c('tot.Q','tot.Q.late.border'),
                                       start_date=c('2020-02-28','2020-02-28'),
                                       end_date=c(Sys.Date(),Sys.Date()))
for(i in 1:nrow(model_table)){
    Make_EpiEstim_Model(model_name=model_table$model_name[i],
                        start_date=model_table$start_date[i],
                        end_date=model_table$end_date[i])
}

for(i in 1:nrow(Q_model_table)){
    Make_EpiEstim_Q_Model(model_name=Q_model_table$model_name[i],
                        start_date=Q_model_table$start_date[i],
                        end_date=Q_model_table$end_date[i])
}



source(here('generate_plots.R'))


