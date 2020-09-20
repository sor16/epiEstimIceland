date <- as.Date('2020-09-09')
mod <- read_rds(here("Results", "Models", "EpiEstim", str_c("Model_",date,'_w_quarantine','',".rds")))
m <- mod$draws() %>% as_draws_df
d <- read_csv("https://docs.google.com/spreadsheets/d/1xgDhtejTtcyy6EN5dbDp5W3TeJhKFRRgm6Xk0s0YFeA/export?format=csv&gid=1788393542",col_types=cols()) %>%
    select(date = Dagsetning, local = Innanlands_Smit,border_1=Landamaeri_Smit_1,border_2=Landamaeri_Smit_2,
           imported = Innflutt_Smit,prop_quarantine=Hlutf_Sottkvi, num_quarantine=Fjoldi_Sottkvi) %>% 
    mutate(date = ymd(date),
           total = local + imported,
           prop_quarantine_and_border=if_else(total!=0,(num_quarantine+border_1+border_2)/total,0)) %>% 
    filter(date >= ymd('2020-02-28') & date<ymd('2020-09-09'))
SI <- get_SI_vec(nrow(d))
d <- mutate(d,lambda=calculate_lambda(total,SI,prop_quarantine))
R_draws <- spread_draws(m, R[day]) %>% 
    group_by(day) %>% 
    mutate(iter = row_number()) %>%
    ungroup %>% 
    select(iter, day, R)

future_R <- function(R_t,t) {
    R_t
}

theme_set(theme_classic(base_size = 12) + 
              theme(legend.position = "none"))
pred_days <- 42

#prop_imported = c(0.00106502,0.001929385,0.002662551)
prop_imported = c(0,0.001929385,0.002662551)

future_prop_quarantine = matrix(0,ncol=4,nrow=41)

pred_days = 42
for(i in 1:4) {
    future_prop_quarantine[,i] = rep(0,pred_days-1)
    
    if(i<4) {
        delay = (i-1)*7
        for(j in 1:(pred_days-1-delay)) {
            future_prop_quarantine[j+delay, i] = 2/(1+exp(-j*12/42))-1
        }
    }
}

name_border <- c('second_test', 'first_test', 'no_test')
name_quarantine <- c('now', 'after_one_week', 'after_two_weeks', 'none')

for(imp in 1:length(prop_imported)) {
    for(quar in 1:ncol(future_prop_quarantine)) {
        plot_dat <- scenario(m=m,d=d, R_draws=R_draws, R_fun = future_R, future_prop_quarantine=future_prop_quarantine[,quar], prop_imported=prop_imported[imp],pred_days=42)
        save(plot_dat, file=here('Results', 'shiny_data', str_c('border_', name_border[imp], '_quarantine_', name_quarantine[quar], '.Rdata')))
    }
}

