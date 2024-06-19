###################
# body.R
# 
# Create the body for the ui. 
# If you had multiple tabs, you could split those into their own
# components as well.
###################
body <- dashboardBody(
  tabItems(
    
    ######################################
    # First tab content - Information Page
    ######################################
    
    tabItem(tabName = "info",
            h2("Information Page"),
            p("The exploreR app is a simple Shiny app that allows you to explore Olink data."),
            uiOutput("howToList"),
            actionButton("addRow", "Add New Step", icon = icon("plus")),
            textAreaInput("newStep", "New Step", value = "", placeholder = "Type new instruction here..."),
            br(),
            p("Credentials:
                Emil Johansson, PhD at KTH Royal Institute of Technology
                Josefin Kenrick, PhD at KTH Royal Institute of Technology
                Jonathan Cedernaes, M.D. PhD at Uppsala University")
    ),
    
    ################################################
    # Second tab content - Data up-loader and viewer
    ################################################
    
    tabItem(tabName = "data",
            h2("Data Preview"),
            fluidRow(
              fileInput("file", "Choose CSV File",
                        accept = ".csv"),
              actionButton("view", "View My Data", icon = icon("eye"))
            ),
            DTOutput("dataTable"), # Output slot for the data table
    ),
    
    ####################
    # Limma part
    ####################
    
    tabItem(tabName = "limma",
            h2("Run Differential Expression Analysis"),
            fluidRow(
              fileInput("dataFile", "Choose CSV File in wide format",
                        accept = ".csv")),
            uiOutput("covariateSelection"),
            uiOutput("differentialVariableSelection"),
            actionButton("runAnalysis", "Run Analysis"),
            DTOutput("resultsTable"), # Output slot for limma results
            
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
    
    #######
    # Plots
    #######
    
    tabItem(tabName = "volcano", h2("Volcano Plot")),
    tabItem(tabName = "boxplot", h2("Boxplot")),
    tabItem(tabName = "scatterplot", h2("Scatterplot")),
    tabItem(tabName = "piechart", h2("Pie Chart")),
    tabItem(tabName = "barplot", h2("Barplot")),
    
    ####################
    # Fourth tab content
    ####################
    
    tabItem(tabName = "session",
            h2("Session Info"),
            verbatimTextOutput("sessionInfo")
    )
  )
)