library(shiny)
library(limma)
library(progressr)
library(DT)

# Define UI for the Shiny app
ui <- fluidPage(
  titlePanel("Differential Expression Analysis using limma"),
  
  sidebarLayout(
    sidebarPanel(
      fileInput("dataFile", "Choose CSV File", accept = ".csv"),
      uiOutput("differentialVariable"),
      uiOutput("covariateSelection"),
      actionButton("runAnalysis", "Run Limma Analysis"),
      br(),
      uiOutput("progress")
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
    
    # generate UI for differential variable selection
    output$differentialVariable <- renderUI({
      if (!is.null(uploadedData())) {
        colnames <- colnames(uploadedData())
        dECols <- colnames[grepl("^[A-Z]", colnames)]
        checkboxGroupInput("differentialVariable", "Select Comparison Variable", 
                           choices = dECols, selected = dECols[1])
      }
    })
    
    # Generate UI for covariate selection based on uploaded data
    output$covariateSelection <- renderUI({
      if (!is.null(uploadedData())) {
        # Filter columns to include only those that start with an uppercase letter
        colnames <- colnames(uploadedData())
        covariateCols <- colnames[grepl("^[A-Z]", colnames)]
        checkboxGroupInput("selectedCovariates", "Select Covariates for Design Matrix", 
                           choices = covariateCols, selected = covariateCols[2:4])
      }
    })
  })
  
  observeEvent(input$runAnalysis, {
    req(input$data)
    req(input$selectedCovariates)
    req(input$differentialVariable)
    
    data <- uploadedData()
    covariates <- input$selectedCovariates
    diff_Variable <- input$differentialVariable
  
    # Separate expression data and design matrix based on user selection
    exprData <- data[, 1:(ncol(data) - length(covariates))]
    exprData <- data |> 
      select(-ID)
    #designMatrix <- data[, covariates, drop = FALSE]
    
    # Ensure the design matrix is a matrix
    designMatrix <- model.matrix(~ ., data = designMatrix)
    
    # Use progressr to show progress bar
    handlers("shiny")
    output$progress <- renderUI({
      progressOutput("progressBar")
    })
    
    withProgressShiny({
      # Run the limma analysis
      fit <- withProgress({
        message("Running lmFit...")
        incProgress(0.4, detail = "Fitting linear model")
        lmFit(exprData, design = design)
      })
      
      fit <- withProgress({
        message("Running eBayes...")
        incProgress(0.4, detail = "Computing statistics")
        eBayes(fit)
      })
      
      results <- withProgress({
        message("Generating topTable...")
        incProgress(0.2, detail = "Creating results table")
        topTable(fit, adjust.method = "BH", number = Inf)
      })
    }, expr = {
      # Output the results as a DataTable
      output$resultsTable <- renderDT({
        datatable(exprData, options = list(pageLength = 10))
        
       })
     })
  })
}

# Function to handle progress in Shiny
withProgressShiny <- function(expr, message = NULL, detail = NULL, value = 0) {
  progressr::with_progress(expr, handlers = progressr::handler_shiny(id = "progressBar", message = message, detail = detail, value = value))
}

# Run the application 
shinyApp(ui = ui, server = server)
