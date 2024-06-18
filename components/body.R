###################
# body.R
# 
# Create the body for the ui. 
# If you had multiple tabs, you could split those into their own
# components as well.
###################
body <- dashboardBody(
  tabItems(
    
    ###################
    # First tab content
    ###################
    
    tabItem(tabName = "info",
            h2("Information Page"),
            p("Details and information about the data and the analysis.")
    ),
    
    ####################
    # Second tab content
    ####################
    
    tabItem(tabName = "data",
            h2("Data Preview"),
            DT::dataTableOutput("dataTable")
    ),
    
    ####################
    # Third tab content
    ####################
    
    tabItem(tabName = "plot",
            fluidRow(
              box(title = "Controls", status = "primary", solidHeader = TRUE, 
                  sliderInput("slider", "Number of observations:", 1, 100, 50)),
              box(title = "Histogram", status = "primary", solidHeader = TRUE, 
                  plotOutput("plot1", height = 250))
            )
    ),
    
    ####################
    # Fourth tab content
    ####################
    
    tabItem(tabName = "session",
            h2("Session Info"),
            verbatimTextOutput("sessionInfo")
    )
  )
)