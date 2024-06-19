###################
# sidebar.R
# 
# Create the sidebar menu options for the ui.
###################
sidebar <- dashboardSidebar(
  sidebarMenu(
    menuItem("Info", tabName = "info", icon = icon("info-circle")),
    menuItem("Data", tabName = "data", icon = icon("database")),
    menuItem("Limma", tabName = "limma", icon = icon("box")),
    menuItem("Plot", tabName = "plot", icon = icon("chart-bar"),
      menuSubItem("Volcano Plot", tabName = "volcano", icon = icon("volcano")),
      menuSubItem("Boxplot", tabName = "boxplot", icon = icon("square")),
      menuSubItem("Scatterplot", tabName = "scatterplot", icon = icon("braille")),
      menuSubItem("Pie Chart", tabName = "piechart", icon = icon("chart-pie")),
      menuSubItem("Barplot", tabName = "barplot", icon = icon("chart-bar"))
  ),
    menuItem("Session", tabName = "session", icon = icon("sd-card"))
  )
)