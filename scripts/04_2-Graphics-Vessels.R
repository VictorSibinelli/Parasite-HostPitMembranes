

######################################################################
#
# Victor Sibinelli (victor.sibinelli@usp.br / sibinelli95@gmail.com)
# 13/07/2024
# Script 03.2 - Graphics - Vessels
#################################################################
library(here)
source(file=here("scripts","03_2-DataAnalysis-Vessels.R"))

## Create a vector of column names
column_names <- colnames(boot_diffs)
# Loop through each column in the boot_diffs data frame using lapply
sapply(seq_along(column_names), function(i) {
  x <- boot_diffs[, i]  # Get the current column data
  col_name <- column_names[i]  # Get the current column name
  
  # Create a density plot
  plot(density(x, na.rm = TRUE), main = "", xlab = "Differences", ylab = "Density")
  
  # Add quantile lines
  abline(v = quantile(x, c(0.025, 0.975), na.rm = TRUE), lwd = 2, col = "black")
  
  # Add a vertical line for the observed value
  abline(v = x[1], col = "red", lwd = 2)
  
  # Calculate p-value and display it on the plot
  p_value <- sum(abs(x) >= abs(x[1]), na.rm = TRUE) / length(x)
  text(x = mean(x, na.rm = TRUE), y = max(density(x)$y, na.rm = TRUE) * 0.9, paste("p-value =",p_value))
  
  # Title for the plot using the current column name
  title(main = col_name)
})