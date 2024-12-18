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


