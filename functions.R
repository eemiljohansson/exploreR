
library(limma)
library(dplyr)
library(purrr)

checkData <- function(data){
  rows_with_na <- sum(apply(is.na(data), 1, any))
  
  # Print number of rows with NA values
  cat("Number of rows with NA values:", rows_with_na, "\n","Not good! Dropping all rows with NA values", "\n")
  
  # Filter out rows with NA values
  revisedData <- data %>%
    filter_all(all_vars(!is.na(.)))
  
  # Print the filtered tibble
  print(revisedData)
} 
  
validateDesignMatrix <- function(designMatrix) {
  # Check if designMatrix is a data frame
  if (!is.data.frame(designMatrix)) {
    stop("Input must be a data frame.")
  }
  
  # Check each column for non-numeric values
  non_numeric_columns <- map_lgl(designMatrix, function(col) {
    any(!is.numeric(col))
  })
  
  # Print warning for columns with non-numeric values
  if (any(non_numeric_columns)) {
    warning(paste("Columns", names(non_numeric_columns)[non_numeric_columns], "contain non-numeric values. Converting to numeric...",'\n'))
  }
  
  # Check for NA values in each column
  na_columns <- purrr::map_lgl(designMatrix, function(col) {
    any(is.na(col))
  })
  
  # Print message if any columns have NA values and drop them
  if (any(na_columns)) {
    warning("Columns with NA values have been dropped.")
    designMatrix <- designMatrix %>%
      select_if(~ !anyNA(.))
  }
  
  reviseddesignMatrix <- designMatrix |> 
    mutate(Sex = case_when(Sex == 'Female' ~ 1,
                           Sex == 'Male' ~ 0),
           VTE = case_when(VTE == 'Yes' ~ 1,
                           VTE == 'No' ~ 0),
           D_dimer_positive = case_when(D_dimer_positive == 'Yes' ~ 1,
                                        D_dimer_positive == 'No' ~ 0))
  
  return(reviseddesignMatrix)
}

run_limma_analysis <- function(expr_data, covariates, design) {
  # Run the limma analysis
  #fit <- lmFit(expr_data, design_matrix)
  fit <- expr_data %>%
    select(where(is.numeric)) %>% 
    t() %>% 
    lmFit(design=design)
  fit <- eBayes(fit)
  results <- topTable(fit, adjust.method = "BH", number = Inf)
  
  return(results)
}

