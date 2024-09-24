###Functions

shuffle_calculation <- function(x, cols, cat) {
  # Copy the data frame
  shuffled <- x
  
  # Shuffle numeric columns (cols) using apply
  shuffled[, cols] <- apply(shuffled[, cols], 2, sample)
  
  # Shuffle the categorical column
  shuffled[, cat] <- sample(shuffled[[cat]])
  
  # Initialize a numeric vector to store the results
  results <- numeric(length(cols))
  
  # Loop through each numeric column to calculate differences in means
  for (i in seq_along(cols)) {
    # Calculate the difference in means for the current numeric column
    results[i] <- diff(tapply(shuffled[[cols[i]]], shuffled[[cat]], mean))
  }
  
  # Assign names to the results vector based on the column names
  names(results) <- cols
  return(results)
}

#functions for shufling data than taking averages by cat
shuffle_means <- function(x, cols, cat, rcat = FALSE, rcol = FALSE) {
  # Copy the data frame
  shuffled <- x
  
  # Shuffle numeric columns (cols) with or without replacement
  shuffled[, cols] <- apply(shuffled[, cols, drop = FALSE], 2, function(col) {
    sample(col, size = length(col), replace = rcol)
  })
  
  # Shuffle the categorical column with or without replacement
  shuffled[[cat]] <- sample(shuffled[[cat]], size = length(shuffled[[cat]]), replace = rcat)
  
  # Calculate the mean for each level of the categorical variable
  results <- tapply(shuffled[[cols]], shuffled[[cat]], mean, na.rm = TRUE)
  results <- results[!is.na(results)]
  if (!is.numeric(results)) {
    results <- as.numeric(results) }
  
  return(results)
}




####Same as shuffle_means 
shuffle_subset <- function(x, cols, index, pairs, shuffle = TRUE) {
  shuffled <- x
  
  # Shuffle the specified numeric column (cols) only if shuffle is TRUE
  if (shuffle) {
    shuffled[, cols] <- apply(shuffled[, cols, drop = FALSE], 2, sample)
  }
  
  ssp_diff <- numeric(length(pairs))  # Initialize result storage
  
  for (i in seq_along(pairs)) {
    pair <- pairs[[i]]  # Get the current pair of species
    
    # Subset the data for the species pair
    pair_data <- shuffled[shuffled[[index]] %in% pair, ]
    
    # Check if both species are present
    if (length(unique(pair_data[[index]])) < 2) {
      cat("No data for species pair:", paste(pair, collapse = " vs "), "\n")
      ssp_diff[i] <- NA_real_  # Assign NA if species are missing
      next  # Skip to the next pair if one species is missing
    }
    
    # Calculate the difference in means for 'cols' between the two species
    means <- tapply(pair_data[[cols]], pair_data[[index]], mean, na.rm = TRUE)
    means <- means[!is.na(means)]
    
    # Ensure there are exactly two means
    if (length(means) == 2) {
      ssp_diff[i] <- diff(means)  # parasite - host
    } else {
      ssp_diff[i] <- NA_real_  # Assign NA if there are not exactly 2 means
    }
  }
  
  # Assign names to the result vector for clarity
  names(ssp_diff) <- sapply(pairs, function(pair) paste(pair, collapse = "_vs_"))
  
  return(ssp_diff)  # Return the result
}

