
rm(list = ls())  # Remove all objects from memory

# Load required libraries and custom functions
library(here)
source(here("scripts", "00-library.R"))  # Load packages
source(here("scripts", "Functions.R"))    # Load custom functions

Mean_data <- read.csv(here("data","processed","Mean_data.csv"))
Median_data <- read.csv(here("data","processed","Median_data.csv"))

# Factor level reordering across all data frames
relevel_factors(ls())

# Define species pairs for comparison
species_pairs <- list(
  c("Psittacanthus robustus", "Vochysia thyrsoidea"),
  c("Phoradendron perrotettii", "Tapirira guianensis"),
  c("Struthanthus rhynchophyllus", "Tipuana tipu")#,
  # c("Viscum album", "Populus nigra")
  )


# Calculate means grouped by 'parasitism'
 Obs_means <-bind_rows(Mean_data %>%
                         group_by(parasitism) %>%
                         summarise(across(everything()[-c(1, 2)], ~ mean(.x, na.rm = TRUE))) %>%
                         rename(Group = parasitism),
                       
                       Mean_data %>%
                         group_by(ssp) %>%
                         summarise(across(everything()[-c(1,11)], ~ mean(.x, na.rm = TRUE)))%>%
                         rename(Group = ssp)
 )
   
 Obs_medians <-bind_rows(Median_data %>%
                         group_by(parasitism) %>%
                         summarise(across(everything()[-c(1, 2)], ~ mean(.x, na.rm = TRUE))) %>%
                         rename(Group = parasitism),
                       
                       Mean_data %>%
                         group_by(ssp) %>%
                         summarise(across(everything()[-c(1,11)], ~ mean(.x, na.rm = TRUE)))%>%
                         rename(Group = ssp)
 )
 
 
data_list <- list(Means = Mean_data, medians = Median_data)
Bootstraped_data <- list()
Variables <- colnames(Mean_data)[sapply(Mean_data, is.numeric)]
Iterations <- 10

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
CI95 <- list()


