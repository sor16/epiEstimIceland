library(shiny)
library(gridExtra)
library(shinythemes)

shinyUI(fluidPage( 
    theme = shinytheme("flatly"),
    
    navbarPage("COVID-19 Sviðsmyndir",
        tabPanel("Sviðsmyndir",
            sidebarLayout(
                sidebarPanel(
                    h3('Sviðsmyndir', style="font-weight:bold;"),
                    radioButtons("test", "Skimun við landamærin:",
                                 c("Engin skimun" = "no_test",
                                    "Ein skimun fyrir alla" = "first_test",
                                   "Tvær skimanir fyrir alla" = "second_test")),
                    radioButtons('quarantine', 'Smitrakning skilar árangri:',
                                 c('Strax' = 'now', 
                                   'Eftir eina viku' = 'after_one_week',
                                   'Eftir tvær vikur' = 'after_two_weeks',
                                   'Aldrei' = 'none'
                                   )
                                 ),
                    br(),
                    h3('Stillingar', style="font-weight:bold;"),
                    radioButtons('xlim', 'Tímabil:',
                                 c('Frá upphafi faraldurs' = 'from_start', 
                                   'Frá upphafi seinni bylgju' = 'from_second_wave',
                                   'Frá deginum í dag' = 'from_today'
                                 )
                                ),
                    checkboxGroupInput('confidence','Öryggismörk',
                                       c('50%' = 50,
                                         '60%' = 60,
                                         '70%' = 70,
                                         '80%' = 80,
                                         '90%' = 90,
                                         '95%' = 95),
                                       selected=c(50,60,70,80,90,95)
                    ),
                    br(),
                    actionButton('go',label="Teikna myndir"),
                    br(),
                    br(),
                    br(),
                    downloadButton('downloadPlot', label = "Hlaða niður mynd sem PDF"),
                    br(),
                    br(),
                    downloadButton('downloadCSV', label='Hlaða niður gögnum sem XLSX')
                ),
                
                mainPanel(
                    plotOutput('epidemia_plot',height = '800px')
                )
            )
        ),
        tabPanel('Um verkefnið',
            mainPanel(
                h4('Um verkefnið')
            )         
        )
    )
))
