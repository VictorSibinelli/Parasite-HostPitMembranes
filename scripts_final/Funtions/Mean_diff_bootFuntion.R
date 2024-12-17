mean_diff_boot <- function(x, cols, cat) {
  shuffled <- x
  
  # Shuffle the categorical column with or without replacement
  shuffled[[cat]] <- sample(shuffled[[cat]], size = length(shuffled[[cat]]), replace = T)
  
  # Calculate the mean for each level of the categorical variable
  results <- diff(tapply(shuffled[[cols]], shuffled[[cat]], mean, na.rm = TRUE))
  
  # Ensure the results are numeric and return as a matrix
  results <- as.matrix(results)
  
  # Remove NA values, if needed
  results <- results[!is.na(results)]
  
  return(results)
}