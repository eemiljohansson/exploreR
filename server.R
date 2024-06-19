###################
# server.R
# 
# For all your server needs 
###################
server <- function(input, output, session) {
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
  
  # Display session info
  output$sessionInfo <- renderPrint({
    sessionInfo()
  })
  
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
}







