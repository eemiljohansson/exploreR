###################
# server.R
# 
# For all your server needs 
###################
server <- function(input, output, session) {
  
  ######
  # Info 
  ######
  
  # Reactive value to store list of steps
  steps <- reactiveVal(c(
    "Upload your data in the data tab",
    "In data you can preview your file",
    "Select the plot tab to view what visualizations you can do",
    "Select the session tab to view your session info"
  ))
  
  # Observe when 'addRow' button is clicked
  observeEvent(input$addRow, {
    # Add the new step if it's not empty
    if (nchar(input$newStep) > 0) {
      steps(c(steps(), input$newStep))
      updateTextAreaInput(session, "newStep", value = "")  # Clear the input area
    }
  })
  
  # Render the list of instructions
  output$howToList <- renderUI({
    howToList <- steps()
    html <- tagList(h3("Pipeline:"))
    for (step in howToList) {
      html <- tagAppendChild(html, tags$p(tags$li(step)))
    }
    return(html)
  })
  
  ##############
  # Data preview
  ##############
  
  # Reactive value to store the dataset
  data <- reactiveVal()
  
  # Observe event for the button click
  observeEvent(input$view, {
    # Check if a file is uploaded
    req(input$file)
    
    # Read the CSV file
    inFile <- input$file
    if (!is.null(inFile)) {
      data(read.csv(inFile$datapath))
    }
  })
  
  # Render the DataTable
  output$dataTable <- renderDT({
    req(data())  # Make sure data is loaded
    datatable(data(), options = list(
      autoWidth = TRUE,  # Automatically adjust the width
      scrollX = TRUE     # Enable horizontal scrolling
    ), 
              class = 'nowrap')  # Keeps the text in a single line for better scrolling
  })
  
  ############
  # Limma part
  ############
  
  # Reactive value to store uploaded data
  uploadedData <- reactiveVal() # Store the uploaded data, NULL was removed
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
  
  #############
  # Sessioninfo
  #############
  
  # Display session info
  output$sessionInfo <- renderPrint({
    sessionInfo()
  })
  
}







