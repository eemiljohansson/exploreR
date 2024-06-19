library(shiny)
library(DT)
library(dplyr)
library(progressr)
# Source the function script
source("run_limma.R")

# Define UI for the Shiny app
ui <- fluidPage(
  titlePanel("Differential Expression Analysis using limma"),
  
  sidebarLayout(
    sidebarPanel(
      fileInput("dataFile", "Choose CSV File", accept = ".csv"),
      uiOutput("covariateSelection"),
      uiOutput("differentialVariableSelection"),
      actionButton("runAnalysis", "Run Analysis"),
      br(),
      textOutput("progress")
    ),
    
    mainPanel(
      DTOutput("resultsTable")
    )
  )
)

# Define server logic for the Shiny app
server <- function(input, output, session) {
  # Reactive value to store uploaded data
  uploadedData <- reactiveVal(NULL)
  
  observeEvent(input$dataFile, {
    req(input$dataFile)
    data <- read.csv(input$dataFile$datapath, row.names = 1)
    uploadedData(data)
    
    # Generate UI for covariate selection based on uploaded data
    output$covariateSelection <- renderUI({
      if (!is.null(uploadedData())) {
        # Filter columns to include only those that start with an uppercase letter
        colnames <- colnames(uploadedData())
        covariateCols <- colnames[grepl("^[A-Z]", colnames) & colnames != "ID"]
        checkboxGroupInput("selectedCovariates", "Select Covariates for Design Matrix", 
                           choices = covariateCols, selected = covariateCols[2:4])
      }
    })
    
    # Generate UI for differential variable selection based on uploaded data
    output$differentialVariableSelection <- renderUI({
      if (!is.null(uploadedData())) {
        colnames <- colnames(uploadedData())
        covariateCols <- colnames[grepl("^[A-Z]", colnames)]
        checkboxGroupInput("differentialVariable", "Select Differential Variable", 
                    choices = covariateCols[1:4], selected = colnames[1])
      }
    })
  })
  
  observeEvent(input$runAnalysis, {
    req(uploadedData())
    req(input$selectedCovariates)
    req(input$differentialVariable)
    
    data <- uploadedData()
    covariates <- input$selectedCovariates
    differentialVariable <- input$differentialVariable
    
    # Separate expression data and design matrix based on user selection
    exprData <- data |> 
      select(-all_of(c(covariates, differentialVariable)))
    designMatrix <- data |>  
      select(all_of(covariates)) |> 
      mutate(!!differentialVariable := data[[differentialVariable]])
    
    # Store ID column separately for later use
    ids <- data$ID
    
    # Use withProgress to show progress bar
    output$progress <- renderText({
      "Running analysis..."
    })
    
    progress <- shiny::Progress$new()
    progress$set(message = "Running limma analysis", value = 0)
    on.exit(progress$close())
    
    # Simulate steps with progress
    progress$inc(0.3, detail = "Fitting model")
    fit <- lmFit(exprData, designMatrix)
    
    progress$inc(0.3, detail = "Computing statistics")
    fit <- eBayes(fit)
    
    progress$inc(0.4, detail = "Creating results table")
    results <- topTable(fit, adjust.method = "BH", number = Inf)
    
    # Output the results as a DataTable
    output$resultsTable <- renderDT({
      datatable(results, options = list(pageLength = 10))
    })
    
    output$progress <- renderText({
      "Analysis complete!"
    })
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
