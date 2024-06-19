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
    revisedData <- checkData(data)
    # Separate expression data and design matrix based on user selection
    exprData <- revisedData |>
      select(-all_of(c(covariates, differentialVariable)))
    designMatrix <- revisedData |>
      select(all_of(covariates)) |>
      mutate(!!differentialVariable := revisedData[[differentialVariable]])
    reviseddesignMatrix <- validateDesignMatrix(designMatrix)
    results <- run_limma_analysis(exprData, covariates, reviseddesignMatrix)
    # Output the results as a DataTable
    output$resultsTable <- renderDT({
      datatable(results, options = list(pageLength = 10))
    })
  })
}