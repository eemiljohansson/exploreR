###################
# sidebar.R
# 
# Create the sidebar menu options for the ui.
###################
sidebar <- dashboardSidebar(
  sidebarMenu(
    menuItem("Info", tabName = "info", icon = icon("info-circle")),
    menuItem("Data", tabName = "data", icon = icon("database")),
    menuItem("Plot", tabName = "plot", icon = icon("chart-bar")),
    menuItem("Session", tabName = "session", icon = icon("sd-card"))
  )
)