#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

#should make a functions script and have limma there and then get input from the user and 


library(shiny)
library(DT)
library(readr)
options(shiny.maxRequestSize = 15 * 1024^2)

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("ExploreR - Proteomics and Metabolomics Analysis"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
          fileInput("exprData", "Choose Expression Data File", accept = ".csv"),
          actionButton("runAnalysis", "Run Analysis"),
          conditionalPanel(
            checkboxGroupInput("show_vars", "Columns of covariates to show:",
                               names(c('Sex','DVT','PE','D_dimer_positive','Age','BMI','D_dimer_FEU','CRP','WBC','Platelets','Hb','EVF','Creatinine'), 
                               selected = names('Sex'))))
        ),

        # Show a plot of the generated distribution
        mainPanel(
          DTOutput("resultsTable")
          
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  observeEvent(input$runAnalysis, {
    req(input$exprData)
   # req(input$designMatrix)
    
    # Read in the expression data and design matrix
    exprData <- read_csv(input$exprData$datapath)
   # designMatrix <- read_csv(input$designMatrix$datapath, row.names = 1)
    
    # Ensure the design matrix is a matrix
    #designMatrix <- as.matrix(designMatrix)
    
    # Run the limma analysis
    #fit <- lmFit(exprData, designMatrix)
    #fit <- eBayes(fit)
    
    #results <- topTable(fit, adjust.method = "BH", number = Inf)
    output$resultsTable <- renderDT({
          datatable(exprData, options = list(pageLength = 10))
    })
  })
}


# Run the application 
shinyApp(ui = ui, server = server)
