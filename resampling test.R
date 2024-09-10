
# List of species pairs
species_pairs <- list(
  c("Psittacanthus robustus", "Vochysia thyrsoidea"),
  c("Phoradendron perrotettii", "Tapirira guianensis"),
  c("Struthanthus rhynchophyllus", "Tipuana tipu")  #,
  #c("Viscum album", "Populus nigra")
)

for (pair in species_pairs) {
  # Subset data for the current species pair
  subset_data <- HydraulicData[HydraulicData$ssp %in% pair, ]
  
  
  
  
}


# Function to permute data by individuals while preserving the species association
permute_data_by_individual_and_species <- function(data, indiv_col, ssp_col) {
  # Create a unique identifier for each individual-species pair
  data$indiv_ssp <- paste(data[[indiv_col]], data[[ssp_col]], sep = "_")
  
  # Identify unique individual-species combinations
  unique_combinations <- unique(data$indiv_ssp)
  
  # Shuffle the unique individual-species combinations
  shuffled_combinations <- sample(unique_combinations)
  
  # Create a copy of the data to permute
  permuted_data <- data
  
  # Reassign the permuted individual-species combinations
  for (i in seq_along(unique_combinations)) {
    # Extract the original and shuffled identifiers
    original_comb <- unique_combinations[i]
    shuffled_comb <- shuffled_combinations[i]
    
    # Get the new indiv and ssp from the shuffled combination
    new_indiv <- strsplit(shuffled_comb, "_")[[1]][1]
    new_ssp <- strsplit(shuffled_comb, "_")[[1]][2]
    
    # Apply the shuffled indiv and ssp to the corresponding rows
    permuted_data[[indiv_col]][permuted_data$indiv_ssp == original_comb] <- new_indiv
    permuted_data[[ssp_col]][permuted_data$indiv_ssp == original_comb] <- new_ssp
  }
  
  # Drop the temporary column
  permuted_data$indiv_ssp <- NULL
  
  return(permuted_data)
}



grouped_means <- HydraulicData %>%
  group_by(indiv) %>%
  summarise(
    ssp = first(ssp),  # Or use any other summary function if needed
    HydraulicDiameter = mean(HydraulicDiameter, na.rm = TRUE),
    vdensity = mean(vdensity, na.rm = TRUE),
    Kmax = mean(Kmax, na.rm = TRUE),
    parasitism = first(parasitism),  # This assumes the same for each individual
    .groups = "drop"
  )


# List of species pairs
species_pairs <- list(
  c("Psittacanthus robustus", "Vochysia thyrsoidea"),
  c("Phoradendron perrotettii", "Tapirira guianensis"),
  c("Struthanthus rhynchophyllus", "Tipuana tipu")  #,
  #c("Viscum album", "Populus nigra")
  
  
diff(tapply(grouped_means$Kmax,grouped_means$parasitism,mean))
Itres <- matrix(data=NA,nrow=10000,ncol=1)  
Itres[1] <- diff(tapply(grouped_means$Kmax,grouped_means$parasitism,mean))  


iterations <- 10000  
  for (i in 2:iterations){
    grouped_means$ssp <- sample(grouped_means$ssp)
    grouped_means$parasitism <- sample(grouped_means$parasitism)
    Itres[i] <- diff(tapply(grouped_means$Kmax,grouped_means$parasitism,mean))

  }
hist(Itres)
abline(v = Itres[1], col = "red")
abline(v = Itres[1]*-1, col = "red")

