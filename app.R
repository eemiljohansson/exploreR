###################
# app.R
# 
# Main controller. 
# Used to import your ui and server components; initializes the app.
###################
library(shiny)
library(shinydashboard)
library(ggplot2)
library(dplyr)
library(DT)  # For displaying data tables interactively
library(limma)
library(purrr)

source('./ui.R')
source('./server.R')
source('./functions.R')

# Run the application 
shinyApp(ui = ui, server = server)