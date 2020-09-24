library(shiny)
library(tidyverse)
library(lubridate)
library(here)
library(ggplot2)
library(writexl)
library(cowplot)
library(scales)
library(plotly)
library(latex2exp)

# future_R <- function(R_t,t) {
#     R_t/(log(t+20)-2)
# }
future_R <- function(R_t,t) {
    R_t
}

icelandic_dates <- function(x) {
    months <- c("janúar", "febrúar", "mars", "apríl", "maí", "júní", 
                "júlí", "ágúst", "september", "október", "nóvember", "desember")
    
    paste0(mday(x), ". ", months[month(x)])
}

#Setting up the Shiny Server
shinyServer(function(input, output, session) {
    get_scenario_dat<-eventReactive(input$go,{
        theme_set(theme_classic(base_size = 12) + 
                      theme(legend.position = "none"))

        load(here('shiny', 'data', str_c('border_', input$test, '_quarantine_', input$quarantine, '.Rdata')))
         plot_dat <- plot_dat %>% filter(prob %in% input$confidence)
        
         start_date = Sys.Date()
         if(input$xlim == 'from_start') start_date = ymd('2020-02-27')
         else if(input$xlim == 'from_second_wave') start_date = ymd('2020-07-20')
    
         p5 <- plot_dat %>%
             filter(name == "y_hat") %>%
             ggplot(aes(date, ymin = lower, ymax = upper)) +
             geom_ribbon(aes(fill = factor(-prob)), alpha = 0.7) +
             geom_point(data = d %>% rename(y_hat = local) %>% pivot_longer(c(y_hat)),
                        inherit.aes = F, aes(x = date, y = value)) +
             geom_vline(xintercept = date, lty = 2) +
             scale_x_date(date_breaks = "month",
                          labels = icelandic_dates,
                          limits = c(start_date, date + 1 + pred_days),
                          expand = expansion(add = 0)) +
             scale_y_continuous(expand = expansion(mult = 0.01),breaks = pretty_breaks(8)) +
             scale_fill_brewer() +
             labs(subtitle = "Ný smit") +
             theme(axis.title = element_blank(),
                   plot.margin = margin(5, 5, 5, 5),
                   legend.title = element_blank())

         p6 <- plot_dat %>%
             filter(name == "R") %>%
             ggplot(aes(date, ymin = lower, ymax = upper)) +
             geom_ribbon(aes(fill = factor(-prob)), alpha = 0.7) +
             geom_hline(yintercept = 1, lty = 2) +
             geom_vline(xintercept = date, lty = 2) +
             scale_x_date(date_breaks = "month",
                          labels = icelandic_dates,
                          limits = c(start_date, date + 1 + pred_days),
                          expand = expansion(add = 0)) +
             scale_y_continuous(expand = expansion(mult = 0.01), breaks = pretty_breaks(8)) +
             scale_fill_brewer() +
             ggtitle(label = waiver(),
                 subtitle = latex2exp::TeX("$R_t$")) +
             theme(axis.title = element_blank(),
                   plot.margin = margin(5, 5, 5, 8))

        plot_grid(p5, p6, ncol = 1)
    }) 
    
    output$epidemia_plot <- renderPlot({
        p = get_scenario_dat()
        if (is.null(p))
            return(NULL)
        p        
    })

    output$downloadPlot <- downloadHandler(
        filename = function() {'COVID_19_svidsmynd.pdf'},
        content = function(file) {
            p = get_scenario_dat()
            ggsave(file, plot = p,width=8, height=6, device = 'pdf')
        }
    )
    
    output$downloadCSV <- downloadHandler(
        filename = function() { "COVID_19_svidsmynd.xlsx"},
        content = function(file) {
            load(here('shiny', 'data', str_c('border_', input$test, '_quarantine_', input$quarantine, '.Rdata')))
            plot_dat <- plot_dat %>% ungroup() %>%
                select(date,name, prob, lower, upper) %>% 
                filter(name %in% c('R', 'y_hat')) %>% 
                mutate(name=if_else(name == 'y_hat', "Spáð smit", "R")) %>%
                rename(Dagsetning = date, Staerd=name, Nedri_mork=lower, Efri_mork=upper, Tegund_spabils=prob)
            write_xlsx(plot_dat, path = file)
        }
    )
})