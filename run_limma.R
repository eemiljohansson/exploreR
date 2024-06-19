# run_limma.R

library(limma)

run_limma_analysis <- function(expr_data, covariates) {
  # Ensure the design matrix is a matrix
  design_matrix <- model.matrix(~ ., data = covariates)
  
  # Run the limma analysis
  fit <- lmFit(expr_data, design_matrix)
  fit <- eBayes(fit)
  results <- topTable(fit, adjust.method = "BH", number = Inf)
  
  return(results)
}
