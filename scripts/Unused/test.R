
# Perform t-tests for vessel wall thickness
for (pair in species_pairs) {
  # Extract data for each species in the pair
  parasite_data <- species_data[[pair[1]]]
  host_data <- species_data[[pair[2]]]
  
  # Check if there is data for both species
  if (nrow(parasite_data) == 0 || nrow(host_data) == 0) {
    cat("No data available for species pair:", pair[1], "and", pair[2], "\n")
    next # Skip this pair if data is missing
  }
  
  # Use tryCatch to handle errors during the t-test
  tryCatch({
    # Perform t-test for vessel wall thickness
    ttest <- t.test(parasite_data$wthickness, host_data$wthickness, var.equal = TRUE)
    
    # Append t-test results to the data frame
    VWall_results <- rbind(VWall_results, data.frame(
      Parasite = pair[1],
      Host = pair[2],
      ParasiteMean = ttest$estimate[1],
      HostMean = ttest$estimate[2],
      MeanDifference = diff(ttest$estimate),
      pvalue = ttest$p.value,
      stringsAsFactors = FALSE
    ))
    row.names(VWall_results) <- NULL
  }, error = function(e) {
    # Handle the error, print a message and continue the loop
    cat("Error in t-test for species pair:", pair[1], "and", pair[2], "\n")
    print(e)
  })
}

# Display the results
print(VWall_results)

# Write results to a CSV file
fwrite(VWall_results, file = here("outputs", "tables", "VWall_results.csv"))
