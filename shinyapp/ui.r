library(shiny)
library(shinydashboard)
library(shinycssloaders)
library(dplyr)

library(shiny)
library(plotly)
library(knitr)
library(tictoc)
library(workflowr)
library(VennDiagram)
library(dplyr)
library(foreach)
library(doParallel)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(devtools)
library(ggsignif)
library(shinydashboard)
library(shinycssloaders)

title=tags$a(href="https://pritchardlabatpsu.github.io/pairwisecomparisons/",tags$img(src="logo4.png",height='45',width='200'), "")
# Define UI for application that plots random distributions
shinyUI(
  dashboardPage(skin="green",
  # pageWithSidebar(
  # Application title
  # headerPanel("Pairwise Comparisons of Conditional Selections Demo!"),
    dashboardHeader(title = title,titleWidth = 250),

  # Sidebar with a slider input for number of observations
  # sidebarPanel(
  dashboardSidebar(
    sliderInput("cohort_size","Please select the Cohort Size",1,10000,500),
    sliderInput("goi_abundance","Please select the abundance of the Gene of Interst in the cohort",1,1000,38),
    sliderInput("pc1_abundance","Please select the abundance of the first Positive Control gene in the cohort",1,1000,200),
    sliderInput("pc2_abundance","Please select the abundance of the second Positive Control gene in the cohort",1,1000,150),
    sliderInput("pc1_goi","Please select the overlap between the positive control and the gene of interest",1,100,10),
    sliderInput("pc1_pc2","Please select the overlap between the two positive controls",1,100,5),
    actionButton('update',"Apply",icon("refresh"),class = "btn btn-primary")
  ),


  # Show a plot of the generated distribution
  dashboardBody(
    fluidRow(
      column(width=4,
        h2("Introduction")
      )
    ),
    fluidRow(
      box(title="",solidHeader = T,width=8,collapsible = T,
          h4(htmlOutput('intro'))),
      infoBoxOutput("warning_info",width = 4)
    ),
    fluidRow(
      column(width=4,
             h2("Inputs")
      )
    ),
    fluidRow(
      column(width=4,
             box(
               title="Positive Control vs Gene of Interest",width = NULL, solidHeader = T,
               # h4(htmlOutput('fet_results_pc1_goi')),
               plotOutput("vennDiagram_pc1goi")
      )
      ),
      column(width=4,
             box(
               title="Positive Control 1 vs Positive Control 2",width = NULL, solidHeader = T,
               # h4(htmlOutput('fet_results_pc1_pc2')),
               plotOutput("vennDiagram_pc1pc2")
      )
      ),
      column(width=4,
             infoBoxOutput("abundance_info",width=NULL),
             infoBoxOutput("pc1goi_info",width=NULL),
             infoBoxOutput("pc1pc2_info",width=NULL)
      )
      ),
    fluidRow(
      column(width=4,
             h2("Results")
      )
    ),
    fluidRow(
      box(width=8,
             box(
               title="Pairwise Comparisons Results", width=NULL, solidHeader = T,
               h4(tags$b(htmlOutput('take_action'))),
               plotOutput('plotlyplot') %>% withSpinner(color="#0dc5c1")
             )
             ),
      column(width=4,
             infoBoxOutput("ks_info",width=NULL)
             )
    ),
    fluidRow(
      column(width=4,
             h2("Conditional Selection Decision Tree")
      )
    ),
    fluidPage(
      box(
        width = 12,
        h4(tags$b("The decision tree below illustrates the various scenarios that pairwise comparisons are able to explore")),
        tags$img(src="decision_tree_transparent.png",height=900,width=900)
      )
    )
  )
)
)
