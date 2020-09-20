border_control = read_excel('/Users/kari/Downloads/border.xlsx', sheet='Sheet1') %>% rename(border_1_sample = border_1, border_2_sample= border_2)

use_quarantine= T

ice_covid_dat <- read_csv("https://docs.google.com/spreadsheets/d/1xgDhtejTtcyy6EN5dbDp5W3TeJhKFRRgm6Xk0s0YFeA/export?format=csv&gid=1788393542",col_types=cols()) %>%
    select(date = Dagsetning, local = Innanlands_Smit,border_1=Landamaeri_Smit_1,border_2=Landamaeri_Smit_2,
           imported = Innflutt_Smit,prop_quarantine=Hlutf_Sottkvi, num_quarantine=Fjoldi_Sottkvi) %>% 
    mutate(date = ymd(date),
           total = local + imported,
           prop_quarantine_and_border=if_else(use_quarantine & total!=0,(num_quarantine+border_1+border_2)/total,0)) %>% 
    filter(date >= ymd("2020-02-28"))


#ein skimun = 15/6 - 12/7
# heimkoma = 13/7 - 19/8
# twofold = 20/8 - ...

ice_covid_dat = inner_join(border_control, ice_covid_dat, by='date')

two_test_ice = ice_covid_dat %>% filter(date>=ymd('2020-07-13'), date <= ymd('2020-08-18'))
two_test_all = ice_covid_dat %>% filter(date>=ymd('2020-08-19'))

border_2_prop_imported_ice = (sum(two_test_ice$border_1)+sum(two_test_ice$border_2))/sum(two_test_ice$border_1_sample, na.rm=T)
border_2_prop_imported_all = (sum(two_test_all$border_1)+sum(two_test_all$border_2))/sum(two_test_all$border_1_sample, na.rm=T)

#0.6 is the amount of tourists - 0.4 its counterpart
no_test_prop_imported = border_2_prop_imported_all
one_test_prop_imported = sum(two_test_all$border_1)/sum(two_test_all$border_1_sample, na.rm=T)
two_test_ice_prop_imported = border_2_prop_imported_all*0.4
two_test_all_prop_imported = 0                              
