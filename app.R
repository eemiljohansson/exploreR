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

source('./ui.R')
source('./server.R')

# Run the application 
shinyApp(ui = ui, server = server)