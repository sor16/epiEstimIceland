library(tidyverse)
library(rstan)
library(lubridate)
library(readxl)
library(here)
library(cmdstanr)
library(posterior)

source('Make_EpiEstim_Q_Model.R')

options(mc.cores = parallel::detectCores())
parallel:::setDefaultClusterOptions(setup_strategy = "sequential")

Make_EpiEstim_Q_Model(model_name='tot.Q.late.border',
                        start_date='2020-02-28',
                        end_date=Sys.Date())
