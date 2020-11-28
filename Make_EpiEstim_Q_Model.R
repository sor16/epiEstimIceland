Make_EpiEstim_Q_Model <- function(model_name, end_date = Sys.Date(),start_date=NULL, warmup = 500, iters = 500, chains = 4, pred_d=NULL) {
    if(is.null(start_date)){
        start_date <- "2020-02-28"
    }
    
    #This data is read in to seperate imported and local cases for the first two months
    pre_dat <- read_csv("https://docs.google.com/spreadsheets/d/1xgDhtejTtcyy6EN5dbDp5W3TeJhKFRRgm6Xk0s0YFeA/export?format=csv&gid=1788393542",col_types=cols()) %>%
         select(date = Dagsetning, local = Innanlands_Smit,border_1=Landamaeri_Smit_1,border_2=Landamaeri_Smit_2, imported = Innflutt_Smit,prop_quarantine=Hlutf_Sottkvi,num_quarantine=Fjoldi_Sottkvi) %>%
         mutate(date = ymd(date),
               total = if_else(date >= ymd("2020-07-23"),local,local + imported),
               prop_quarantine=if_else(total!=0,num_quarantine/total,0)) %>%
         filter(date >= ymd('2020-02-28') & date <= ymd('2020-04-08'))
    
    diagnosis_dat <- read_sheet("https://docs.google.com/spreadsheets/d/1gbn2CUSPFExQ4mdHsC6Iq671N46ogS8G4mGgISb4x4g/edit#gid=1720508877", range = "FjoldiSmitaEftirDogum") %>%
        rename(date = ...1,
               lsh = "Sýkla- og veirufræðideild LSH",
               ie = "Íslensk erfðagreining",
               border = "Landamæraskimun") %>%
        mutate(date=dmy(date))
    quarantine_dat <-  read_sheet("https://docs.google.com/spreadsheets/d/1gbn2CUSPFExQ4mdHsC6Iq671N46ogS8G4mGgISb4x4g/edit#gid=1720508877", range = "SottkviVidVeikindiAnLandamaera") %>%
        head(-1) %>% 
        rename(prop_quarantine = "Eru í sóttkví við greiningu",
               date = "Dags.") %>%
        mutate(prop_quarantine = prop_quarantine/100, date=dmy(date)) %>%
        select(date, prop_quarantine) 

    d <- inner_join(diagnosis_dat, quarantine_dat, by="date") %>%
        mutate(local = lsh+ie,
               total = if_else(date >= ymd("2020-07-23"),local,local + border),
               local= c(pre_dat$local, local[date > ymd('2020-04-08')]),
               imported = c(pre_dat$imported, border[date > ymd('2020-04-08')]),
               num_quarantine = round(prop_quarantine*local)) %>%
        select(date, local, border, total, imported, prop_quarantine, num_quarantine)

    if(!is.null(pred_d)) {
        pred_d <- pred_d %>% mutate(date=ymd(end_date)+day) %>% select(date, local,border_1,border_2, imported, prop_quarantine, num_quarantine, total)
        d <- d %>% bind_rows(pred_d)
    }
    
    # shape <- 4.5
    # rate <- 0.6
    shape <- 1.54
    rate <- 0.28
    
    total <- d$total
    
    N_days <- length(total)
    
    local <- d$local
    total <- d$total
    N_days <- length(local)
    prop_quarantine <- d$prop_quarantine
    
    stan_data <- list(
        N_days = N_days,
        SI_shape = 1.54,
        SI_rate = 0.28,
        local = local,
        total = total,
        prop_quarantine=prop_quarantine
    )
    
    mod <- cmdstan_model(here("Stan", "EpiEstim_Q.stan"))
    
    fit <- mod$sample(
        data = stan_data, 
        show_messages = FALSE, 
        chains = chains, 
        parallel_chains = chains,
        iter_sampling = iters,
        iter_warmup = warmup,
        max_treedepth = 15,
        init = 0,
        refresh = 10
    )
    
    fit$save_object(file = here("Results", "EpiEstim", str_c(model_name,'_', end_date,".rds")))
    
}
