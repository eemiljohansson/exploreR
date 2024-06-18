# Function to create age distribution boxplots
visualize_age_distribution <- function(df = input$file, age_column = "Age", cohort_column = "VTE", fill_color = "VTE", 
                                       title = "Age distribution of the cohorts used in this study", 
                                       subtitle = "Cardiovascular, metabolic and healthy cohort(s)") {
  # Check for necessary columns
  if (!all(c(age_column, cohort_column) %in% names(df))) {
    stop("Required columns are missing from the dataframe")
  }
  
  # Filter for non-missing age values
  filtered_data <- df %>%
    filter(!is.na(!!sym(age_column)))
  
  # Create the boxplot
  plot <- filtered_data %>%
    ggplot(aes_string(x = cohort_column, y = age_column, fill = fill_color)) +
    geom_boxplot() +
    labs(x = "Cohort", y = "Age at sampling/diagnosis", title = title, subtitle = subtitle) +
    theme_light() +
    guides(fill = FALSE) +
    coord_flip()
  
  return(plot)
}
