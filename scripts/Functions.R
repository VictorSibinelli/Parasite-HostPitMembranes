######################################################################
#
# Victor Sibinelli (victor.sibinelli@usp.br / sibinelli95@gmail.com)
# 13/07/2024
# Script - Functions
######################################################################



mean_diff_boot <- function(x, cols, cat) {
  shuffled <- x
  
  # Shuffle the categorical column with or without replacement
  shuffled[[cols]] <- sample(shuffled[[cols]], size = length(shuffled[[cols]]), replace = T)
  
  # Calculate the mean for each level of the categorical variable
  results <- diff(tapply(shuffled[[cols]], shuffled[[cat]], mean, na.rm = TRUE))
  
  # Ensure the results are numeric and return as a matrix
  results <- as.matrix(results)
  
  # Remove NA values, if needed
  results <- results[!is.na(results)]
  
  return(results)
}

bootstrap_summary <- function(data, group_var, variables) {
  data %>%
    group_by(across(all_of(group_var))) %>%
    summarise(across(
      all_of(variables),
      ~ mean(sample(.x, replace = TRUE), na.rm=T),
      .names = "{col}"
    ), .groups = "drop")
}

shuffle_means <- function(x, cols, cat) {
  # Shuffle numeric column (cols) with replacement
  x[[cols]] <- sample(x[[cols]], size = nrow(x), replace = TRUE)
  
  # Calculate the mean for each level of the categorical variable
  results <- tapply(x[[cols]], x[[cat]], mean, na.rm = TRUE)
  
  # Convert to a named vector or data frame
  results <- data.frame(Category = names(results), Mean = as.numeric(results))
  
  return(results)
}

relevel_factors <- function(dataframes) {
  for (df_name in dataframes) {
    # Retrieve the data frame by its name
    df <- get(df_name, envir = .GlobalEnv)
    
    if (!is.data.frame(df)) next  # Skip non-data frame objects
    
    # Relevel the 'ssp' column
    if ("ssp" %in% colnames(df)) {  
      df$ssp <- factor(df$ssp, levels = c(
        "Psittacanthus robustus", "Vochysia thyrsoidea",
        "Phoradendron perrotettii", "Tapirira guianensis",
        "Struthanthus rhynchophyllus", "Tipuana tipu",
        "Viscum album", "Populus nigra"
      ))
    }
    
    # Relevel the 'Grouping' column
    if ("Grouping" %in% colnames(df)) {
      df$Grouping <- factor(df$Grouping, levels = c(
        "Parasite", "Host",
        "Psittacanthus robustus", "Vochysia thyrsoidea",
        "Phoradendron perrotettii", "Tapirira guianensis",
        "Struthanthus rhynchophyllus", "Tipuana tipu",
        "Viscum album", "Populus nigra"
      ))
    }
    
    # Relevel the 'parasitism' column
    if ("parasitism" %in% colnames(df)) {  
      df$parasitism <- factor(df$parasitism, levels = c("Parasite", "Host"))
    }
    
    # Save the modified data frame back to the global environment
    assign(df_name, df, envir = .GlobalEnv)
  }
}
perform_t_tests <- function(data, trait_col, indiv_test) {
  # Ensure 'indiv' column exists
  if (!"indiv" %in% colnames(data)) {
    stop(paste("Error: 'indiv' column not found in dataset for trait:", trait_col))
  }
  
  results_df <- data.frame(
    indiv_tested = character(),
    difference_observed = numeric(),
    CI95 = character(),
    pvalue = numeric(),
    Rdiff = numeric(),  # New column
    stringsAsFactors = FALSE
  )
  
  for (ind in indiv_test) {
    filtered_data <- data %>% filter(indiv %in% ind)
    
    if (length(unique(filtered_data$indiv)) == 2) {
      group1 <- filtered_data[[trait_col]][filtered_data$indiv == ind[1]]
      group2 <- filtered_data[[trait_col]][filtered_data$indiv == ind[2]]
      
      if (length(group1) > 1 & length(group2) > 1) {
        result <- t.test(group1, group2, var.equal = FALSE)
        indiv_tested <- paste(ind[1], "vs", ind[2])
        difference_observed <- result$estimate[1] - result$estimate[2]
        CI95 <- paste0("(", round(result$conf.int[1], 2), ", ", round(result$conf.int[2], 2), ")")
        pvalue <- result$p.value
        mean_group2 <- mean(group2, na.rm = TRUE)  # Compute mean of group 2
        print(result)
        # Compute Rdiff safely (avoid division by zero)
        Rdiff <- ifelse(mean_group2 != 0, (difference_observed / mean_group2)*100, NA)
        
        results_df <- rbind(
          results_df,
          data.frame(
            indiv_tested = indiv_tested,
            difference_observed = round(difference_observed, 3),
            CI95 = CI95,
            pvalue = round(pvalue, 5),
            Rdiff = round(Rdiff, 3),  # Store Rdiff in the results
            stringsAsFactors = FALSE
          )
        )
      } else {
        cat("\nSkipping pair:", ind, "- Not enough data for valid t-test.\n")
      }
    } else {
      cat("\nSkipping pair:", ind, "- Insufficient or invalid data.\n")
    }
  }
  
  return(results_df)
}


