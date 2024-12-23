
rm(list = ls())  # Remove all objects from memory

# Load required libraries and custom functions
library(here)
library(tidyverse)  # Load packages
source(here("scripts", "Functions.R"))    # Load custom functions

Mean_data <- read.csv(here("data","processed","Mean_data.csv"))
Median_data <- read.csv(here("data","processed","Median_data.csv"))

# Factor level reordering across all data frames
relevel_factors(ls())

# Define species pairs for comparison
species_pairs <- list(
  c("Psittacanthus robustus", "Vochysia thyrsoidea"),
  c("Phoradendron perrotettii", "Tapirira guianensis"),
  c("Struthanthus rhynchophyllus", "Tipuana tipu"),
  c("Viscum album", "Populus nigra")
  )



 
data_list <- list(Means = Mean_data, Medians = Median_data)
Bootstraped_data <- list()
Variables <- colnames(Mean_data)[sapply(Mean_data, is.numeric)]
Iterations <- 250000

for (d in seq_along(data_list)){
  
  df_name <- names(data_list)[d]
  m <- matrix(NA, nrow = Iterations, ncol = length(Variables))
  colnames(m) <- Variables  # Set column names for the matrix
  x <- data_list[[d]]
  for (v in Variables) {
    # Generate bootstrap values for the current variable
    boot <- replicate(Iterations, mean_diff_boot(x, cols = v, cat = "parasitism"))
    
    # Assign the generated bootstrap values to the corresponding column in the matrix
    m[, v] <- boot
  }
  
  # Store the matrix in the Bootstraped_data list with appropriate names
  Bootstraped_data[[paste(df_name, "parasitism", sep = "_")]] <- m
  
  rm(m)
}


Iterations <- 100000
 

for (d in seq_along(data_list)) {
  df_name <- names(data_list)[d]
  x <- data_list[[d]]
  
  for (pair in species_pairs) {
    # Subset data for the species pair
    subset_data <- subset(x, ssp %in% pair)
    pair_id <- paste(pair, collapse = "X")
    
    # Initialize a matrix for bootstrap values
    m <- matrix(NA, nrow = Iterations, ncol = length(Variables))
    colnames(m) <- Variables  # Set column names for the matrix
    
    for (v in Variables) {
      # Generate bootstrap values for the current variable
      boot <- replicate(Iterations, mean_diff_boot(subset_data, cols = v, cat = "ssp"))
      
      # Assign the generated bootstrap values to the corresponding column in the matrix
      m[, v] <- boot
      
    }
    
    # Store the matrix in the Bootstraped_data list with appropriate names
    Bootstraped_data[[paste(df_name, pair_id, sep = "_")]] <- m
    
  }
  rm(x,m)
}




###################################################################
Iterations <- 10000
CI95 <- list()

for (d in seq_along(data_list)) {
  df_name <- names(data_list)[d]
  x <- data_list[[d]]
  
  for (p in levels(x$parasitism)) {
    subset_data <- subset(x, parasitism == p)
    para_type <- paste(p, "CI95", sep = "_")
    
    # Generate bootstrap values for all variables in one matrix
    m <- sapply(Variables, function(v) {
      replicate(Iterations, mean(sample(subset_data[[v]], replace = TRUE), na.rm=T))
    })
    
    # Add the result matrix to the list
    CI95[[paste(df_name, para_type, sep = "_")]] <- m
  }
}




for (d in seq_along(data_list)) {
  df_name <- names(data_list)[d]
  x <- data_list[[d]]
  
  for (p in levels(x$ssp)) {
    subset_data <- subset(x, ssp == p)
    para_type <- paste(p, "CI95", sep = "_")
    
    # Generate bootstrap values for all variables in one matrix
    m <- sapply(Variables, function(v) {
      replicate(Iterations, mean(sample(subset_data[[v]], replace = TRUE),na.rm=T))
    })
    
    # Add the result matrix to the list
    CI95[[paste(df_name, para_type, sep = "_")]] <- m
  }
}

Boots <- list(Bootstraped_data,CI95)


lapply(seq_along(Boots), function(i) {
  # Extract the sublist and its name
  sublist <- Boots[[i]]
  
  # Iterate over the data frames in the sublist
  lapply(seq_along(sublist), function(j) {
    # Extract the data frame and its name
    df <- sublist[[j]]
    df_name <- names(sublist)[j]
    
    # Create a file name using sublist and data frame names
    table_name <- paste0(df_name, ".csv")

    # Save the data frame
    write.csv(df, file = here("data", "processed", "ressampled", table_name))
  })
})


