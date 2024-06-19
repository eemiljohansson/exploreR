
library(limma)
library(dplyr)

run_limma_analysis <- function(expr_data, covariates) {
  # Run the limma analysis
  #fit <- lmFit(expr_data, design_matrix)
  fit <- expr_data %>%
    select(where(is.numeric)) %>% 
    t() %>% 
    lmFit(design=designMatrix)
  fit <- eBayes(fit)
  results <- topTable(fit, adjust.method = "BH", number = Inf)
  
  return(results)
}


validateDesignMatrix <- function(designMatrix) {
  # Check if designMatrix is a data frame
  if (!is.data.frame(designMatrix)) {
    stop("Input must be a data frame.")
  }
  
  # Convert non-numeric columns to factors
  designMatrix <- designMatrix %>%
    mutate(across(everything(), ~ ifelse(is.numeric(.), ., as.factor(.))))
  
  return(designMatrix)
}
