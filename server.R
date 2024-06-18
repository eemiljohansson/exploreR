###################
# server.R
# 
# For all your server needs 
###################
server <- function(input, output) {
  # Data generation for demonstration
  set.seed(122)
  histdata <- rnorm(500)
  
  output$plot1 <- renderPlot({
    data <- histdata[seq_len(input$slider)]
    hist(data, main = "Histogram of Data", col = 'darkgray', border = 'white')
  })
  
  # Display session info
  output$sessionInfo <- renderPrint({
    sessionInfo()
  })
  
  # Example data table
  output$dataTable <- DT::renderDataTable({
    head(mtcars)
  })
}