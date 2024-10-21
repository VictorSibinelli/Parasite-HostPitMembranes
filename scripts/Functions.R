######################################################################
#
# Victor Sibinelli (victor.sibinelli@usp.br / sibinelli95@gmail.com)
# 13/07/2024
# Script - Functions
######################################################################

# ------------------------------------------------------
# Function to plot density before and after outlier removal on the same graph
# ------------------------------------------------------
plot_density_comparison <- function(original_density, updated_density, species, variable) {
  plot(original_density, 
       main = paste(species, "-", variable), 
       xlab = variable, 
       ylab = "Density", 
       lwd = 2, 
       col = "red", 
       xlim = range(c(original_density$x, updated_density$x), na.rm = TRUE), 
       cex.main = 1.2, 
       cex.lab = 1.0, 
       cex.axis = 0.8)  # Adjusted text sizes
  
  lines(updated_density, lwd = 2, col = "blue")
  legend("topright", legend = c("Before", "After"), col = c("red", "blue"), lwd = 2, cex = 0.8)  # Smaller legend text
}
relevel_factors <- function(dataframes) {
  for (df_name in dataframes) {
    df <- get(df_name)  # Get the object by name
    
    if (!is.data.frame(df)) next  # Skip non-data frames
    
    # Check for the 'ssp' or 'Grouping' column and reorder factors
    if ("ssp" %in% colnames(df)) {  
      df$ssp <- factor(df$ssp, levels = c(
        "Psittacanthus robustus", "Vochysia thyrsoidea",
        "Phoradendron perrotettii", "Tapirira guianensis",
        "Struthanthus rhynchophyllus", "Tipuana tipu",
        "Viscum album", "Populus nigra"
      ))
    }
    if ("Grouping" %in% colnames(df)){
      df$Grouping <- factor(df$Grouping, levels = c(
        "Parasite","Host",
        "Psittacanthus robustus", "Vochysia thyrsoidea",
        "Phoradendron perrotettii", "Tapirira guianensis",
        "Struthanthus rhynchophyllus", "Tipuana tipu",
        "Viscum album", "Populus nigra"
      ))
    }
    # Check for 'parasitism' column and reorder factors
    if ("parasitism" %in% colnames(df)) {  
      df$parasitism <- factor(df$parasitism, levels = c("Parasite", "Host"))
    }
    
    assign(df_name, df, envir = .GlobalEnv)  # Save the modified data frame
  }
}

# Now pass the list of data frame names to the function


# ------------------------------------------------------
# Function to replace outliers with NA and test for normality
# ------------------------------------------------------
replace_outliers_test_normality <- function(data, variable, species, max_iter = 10) {
  outlier_count <- 0
  original_data <- data[[variable]][data$ssp == species]
  shapiro_p <- NA
  
  for (i in 1:max_iter) {
    outliers <- boxplot.stats(data[[variable]][data$ssp == species])$out
    if (length(outliers) == 0) break
    
    outlier_count <- length(outliers)  # Update outlier count
    data[[variable]][data[[variable]] %in% outliers & data$ssp == species] <- NA
    
    shapiro_test <- shapiro.test(data[[variable]][data$ssp == species])
    shapiro_p <- shapiro_test$p.value
    if (shapiro_p > 0.05) break
  }
  
  original_density <- density(original_data, na.rm = TRUE)
  updated_density <- density(data[[variable]][data$ssp == species], na.rm = TRUE)
  
  return(list(data = data, 
              outlier_count = outlier_count, 
              shapiro_p = shapiro_p, 
              original_density = original_density, 
              updated_density = updated_density))
}

# ------------------------------------------------------
# Function to perform normality test across species and variables
# ------------------------------------------------------
perform_normality_test <- function(data, species_list, variables) {
  results <- data.frame(species = character(), stringsAsFactors = FALSE)
  
  for (ssp in species_list) {
    shapiro_results <- sapply(variables, function(var) {
      if (sum(!is.na(data[[var]][data$ssp == ssp])) > 2) {
        return(shapiro.test(data[[var]][data$ssp == ssp])$p.value)
      } else {
        return(NA)
      }
    })
    results <- rbind(results, c(ssp, shapiro_results))
  }
  
  colnames(results) <- c("species", variables)
  return(results)
}

# ------------------------------------------------------
# Function for shuffling data and calculating differences in means
# ------------------------------------------------------
shuffle_calculation <- function(x, cols, cat) {
  shuffled <- x
  
  # Shuffle numeric columns (cols)
  shuffled[, cols] <- apply(shuffled[, cols], 2, sample)
  
  # Shuffle the categorical column
  shuffled[, cat] <- sample(shuffled[[cat]])
  
  # Calculate the differences in means for each numeric column
  results <- numeric(length(cols))
  
  for (i in seq_along(cols)) {
    results[i] <- diff(tapply(shuffled[[cols[i]]], shuffled[[cat]], mean))
  }
  
  names(results) <- cols
  return(results)
}

shuffle_means <- function(x, cols, cat, rcat = FALSE, rcol = FALSE) {
  shuffled <- x
  
  # Shuffle numeric columns (cols) with or without replacement
  shuffled[, cols] <- apply(shuffled[, cols, drop = FALSE], 2, function(col) {
    sample(col, size = length(col), replace = rcol)
  })
  
  # Shuffle the categorical column with or without replacement
  shuffled[[cat]] <- sample(shuffled[[cat]], size = length(shuffled[[cat]]), replace = rcat)
  
  # Calculate the mean for each level of the categorical variable
  results <- tapply(shuffled[[cols]], shuffled[[cat]], mean, na.rm = TRUE)
  
  # Ensure the results are numeric and return as a matrix
  results <- as.matrix(results)
  
  # Remove NA values, if needed
  results <- results[!is.na(results)]
  
  return(results)
}

# ------------------------------------------------------
# Function for shuffling data subsets and calculating species pair differences
# ------------------------------------------------------
shuffle_subset <- function(x, cols, index, pairs, shuffle = TRUE) {
  shuffled <- x
  
  if (shuffle) {
    shuffled[, cols] <- apply(shuffled[, cols, drop = FALSE], 2, sample)
  }
  
  ssp_diff <- numeric(length(pairs))  # Initialize result storage
  
  for (i in seq_along(pairs)) {
    pair <- pairs[[i]]  # Get the current pair of species
    
    # Subset the data for the species pair
    pair_data <- shuffled[shuffled[[index]] %in% pair, ]
    
    if (length(unique(pair_data[[index]])) < 2) {
      cat("No data for species pair:", paste(pair, collapse = " vs "), "\n")
      ssp_diff[i] <- NA_real_  # Assign NA if species are missing
      next
    }
    
    # Calculate the difference in means for 'cols' between the two species
    means <- tapply(pair_data[[cols]], pair_data[[index]], mean, na.rm = TRUE)
    means <- means[!is.na(means)]
    
    if (length(means) == 2) {
      ssp_diff[i] <- diff(means)  # Parasite - host
    } else {
      ssp_diff[i] <- NA_real_  # Assign NA if there are not exactly 2 means
    }
  }
  
  names(ssp_diff) <- sapply(pairs, function(pair) paste(pair, collapse = "_vs_"))
  return(ssp_diff)  # Return the result
}
