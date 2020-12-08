library(tidyverse)
library(rstan)
library(lubridate)
library(readxl)
library(here)
library(cmdstanr)
library(posterior)
library(readr)
library(googlesheets4)

options(gargle_oauth_email = "karirogg@gmail.com")
gs4_auth(email = "karirogg@gmail.com")

#set_cmdstan_path('~/software/cmdstan')
set_cmdstan_path('/nfs/bin/cmdstan-2.24.0') 

source('Make_EpiEstim_Q_Model.R')

options(mc.cores = parallel::detectCores())
parallel:::setDefaultClusterOptions(setup_strategy = "sequential")

Make_EpiEstim_Q_Model(model_name='tot.Q.late.border',
                      start_date='2020-02-28',
                      end_date=Sys.Date(),
                      pred_d = NULL)
