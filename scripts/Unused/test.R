


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




# Perform resampling to create null distribution --> If it is the first time runing this script,
#uncoment lines 54-59
WT_resample <- t(replicate(iterations,
                           shuffle_means(x = WT_EV, cols = "wthickness", cat = "parasitism")))


WT_resample[1, ] <- WT_obs  # First row: observed values
fwrite(WT_resample, file = here("data", "processed", "ressampled", "WallThickness_ressampled.csv")) 

WT_resample <- as.matrix(read.csv(file = here(
  "data", "processed", "ressampled", "WallThickness_ressampled.csv")))

# Calculate differences between groups
WT_diff <- WT_resample[, 1] - WT_resample[, 2]  

# --------------------------------------------------------
# Calculate p-value from resampling distribution
# --------------------------------------------------------
WT_pvalue <- sum(abs(WT_diff) >= abs(WT_diff[1])) / length(WT_diff)

# Plot resampling density and mark the observed difference
plot(density(WT_diff), main = "Resampling Density Plot")
abline(v = WT_diff[1], col = "red")  # Observed difference line
abline(v = quantile(WT_diff, c(0.025, 0.975)), col = "black", lty = 2)  # 95% CI
text(x = 0.45, y = 0.5,cex=0.9, paste("p-value =", round(WT_pvalue,digits = 3)))
