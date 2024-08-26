######################################################################
#
# Victor Sibinelli (victor.sibinelli@usp.br / sibinelli95@gmail.com)
# 13/07/2024
# Script 03.3 - Data Analysis - Pit membranes
#################################################################
library(here)
source(here("scripts", "02_3-TestAssumptions-Pit.R"))
rm(list = ls())

# load data
pitdata_clean <- read.csv(here("data", "processed", "pitdata_clean.csv"))
pitOdata_clean <- read.csv(here("data", "processed", "pitO_clean.csv"))

# List of data frame names
dataframes <- ls()
# Relevel the factors for each data frame
for (df_name in dataframes) {
  df <- get(df_name) # Get the data frame by name

  if ("ssp" %in% colnames(df)) { # Check if 'ssp' column exists
    df$ssp <- factor(df$ssp, levels = c(
      "Psittacanthus robustus", "Vochysia thyrsoidea",
      "Phoradendron perrotettii", "Tapirira guianensis",
      "Struthanthus rhynchophyllus", "Tipuana tipu",
      "Viscum album", "Populus nigra"
    ))
    assign(df_name, df) # Assign the modified data frame back to its name
  }
  rm(df, df_name) # remove duplicated dataframe
}

# Initialize an empty data frame to store the results
results <- data.frame(ssp = character(),
                      mean_pcavg = numeric(),
                      mean_peavg = numeric(),
                      means_difference = numeric(),
                      p_value = numeric(),
                      stringsAsFactors = FALSE)

# Loop through each unique species in the dataset
for (ssp in unique(pitdata_clean$ssp)) {
  # Subset data for the current species
  data <- pitdata_clean[pitdata_clean$ssp == ssp, ]
  
  # Perform Welch's t-test comparing pcavg and peavg
  test_result <- t.test(data$pcavg, data$peavg, var.equal = FALSE)
  
  # Extract the values
  p_value <- test_result$p.value
  mean_pcavg <- test_result$estimate[1]*1000
  mean_peavg <- test_result$estimate[2]*1000
  means_diff <- mean_pcavg - mean_peavg
  
  # Store the results in the data frame
  results <- rbind(results, data.frame(ssp = ssp,
                                       mean_pcavg = mean_pcavg,
                                       mean_peavg = mean_peavg,
                                       means_difference = means_diff,  # Corrected here
                                       p_value = p_value))
}

# View the results
print(results)
fwrite(results, file = here("outputs", "tables", "pcavgXpeavg.csv"))
##### All samples showed slightly thicker centers

# Data subsets
species_data <- list(
  "P. robustus" = pitdata_clean[pitdata_clean$ssp == "Psittacanthus robustus", ],
  "V. thyrsoidea" = pitdata_clean[pitdata_clean$ssp == "Vochysia thyrsoidea", ],
  "P. perrotettii" = pitdata_clean[pitdata_clean$ssp == "Phoradendron perrotettii", ],
  "T. guianensis" = pitdata_clean[pitdata_clean$ssp == "Tapirira guianensis", ],
  "S. rhynchophyllus" = pitdata_clean[pitdata_clean$ssp == "Struthanthus rhynchophyllus", ],
  "T. tipu" = pitdata_clean[pitdata_clean$ssp == "Tipuana tipu", ],
  "V. album" = pitdata_clean[pitdata_clean$ssp == "Viscum album", ],
  "P. nigra" = pitdata_clean[pitdata_clean$ssp == "Populus nigra", ]
)

species_data2 <- list(
  "P. robustus" = pitOdata_clean[pitOdata_clean$ssp == "Psittacanthus robustus", ],
  "V. thyrsoidea" = pitOdata_clean[pitOdata_clean$ssp == "Vochysia thyrsoidea", ],
  "P. perrotettii" = pitOdata_clean[pitOdata_clean$ssp == "Phoradendron perrotettii", ],
  "T. guianensis" = pitOdata_clean[pitOdata_clean$ssp == "Tapirira guianensis", ],
  "S. rhynchophyllus" = pitOdata_clean[pitOdata_clean$ssp == "Struthanthus rhynchophyllus", ],
  "T. tipu" = pitOdata_clean[pitOdata_clean$ssp == "Tipuana tipu", ],
  "V. album" = pitOdata_clean[pitOdata_clean$ssp == "Viscum album", ],
  "P. nigra" = pitOdata_clean[pitOdata_clean$ssp == "Populus nigra", ]
)

# Initialize result data frames
pit_mean_diff <- data.frame(Parasite = character(),
                            Host = character(),
                            ParasiteMean = numeric(),
                            HostMean = numeric(),
                            MeanDifference = numeric(),
                            pvalue = numeric(),
                            stringsAsFactors = FALSE)

pcd_results <- data.frame(Parasite = character(),
                          Host = character(),
                          ParasiteMean = numeric(),
                          HostMean = numeric(),
                          MeanDifference = numeric(),
                          pvalue = numeric(),
                          stringsAsFactors = FALSE)

pitDiameter_results <- data.frame(Parasite = character(),
                                  Host = character(),
                                  ParasitePitDiameter = numeric(),
                                  HostPitDiameter = numeric(),
                                  MeanDifference = numeric(),
                                  pvalue = numeric(),
                                  stringsAsFactors = FALSE)

pitOpening_results <- data.frame(Parasite = character(),
                                  Host = character(),
                                  ParasitePitOpeningDiameter = numeric(),
                                  HostPitOpeningDiameter = numeric(),
                                  MeanDifference = numeric(),
                                  pvalue = numeric(),
                                  stringsAsFactors = FALSE)

# Perform t-tests for Pit membrane thickness
for (pair in list(c("P. robustus", "V. thyrsoidea"), 
                  c("P. perrotettii", "T. guianensis"), 
                  c("S. rhynchophyllus", "T. tipu"), 
                  c("V. album", "P. nigra"))) {
  
  parasite_data <- species_data[[pair[1]]]
  host_data <- species_data[[pair[2]]]
  
  # Perform t-test for pit membrane thickness
  ttest <- t.test(parasite_data$pitavg, host_data$pitavg, var.equal = FALSE)
  pit_mean_diff <- rbind(pit_mean_diff, data.frame(
    Parasite = pair[1],
    ParasiteMean = ttest$estimate[1]*1000,
    Host = pair[2],
    HostMean = ttest$estimate[2]*1000,
    MeanDifference = diff(ttest$estimate)*1000,
    pvalue = ttest$p.value,
    stringsAsFactors = FALSE
  ))
}

# Save pit membrane thickness results to CSV
pit_mean_diff
fwrite(pit_mean_diff, file = here("outputs", "tables", "pit_membrane_diff.csv"))

# Perform t-tests for Pit chamber depth
for (pair in list(c("P. robustus", "V. thyrsoidea"), 
                  c("P. perrotettii", "T. guianensis"), 
                  c("S. rhynchophyllus", "T. tipu"), 
                  c("V. album", "P. nigra"))) {
  
  parasite_data <- species_data[[pair[1]]]
  host_data <- species_data[[pair[2]]]
  
  # Perform t-test for pit chamber depth
  ttest <- t.test(parasite_data$pcd, host_data$pcd, var.equal = FALSE)
  pcd_results <- rbind(pcd_results, data.frame(
    Parasite = pair[1],
    ParasiteMean = ttest$estimate[1]*1000,
    Host = pair[2],
    HostMean = ttest$estimate[2]*1000,
    MeanDifference = diff(ttest$estimate)*1000,
    pvalue = ttest$p.value,
    stringsAsFactors = FALSE
  ))
}

# Save pit chamber depth results to CSV
pcd_results
fwrite(pcd_results, file = here("outputs", "tables", "pcd_results.csv"))


# Perform t-tests for Pit Diameter
for (pair in list(c("P. robustus", "V. thyrsoidea"), 
                  c("P. perrotettii", "T. guianensis"), 
                  c("S. rhynchophyllus", "T. tipu"), 
                  c("V. album", "P. nigra"))) {
  
  parasite_data <- species_data2[[pair[1]]]
  host_data <- species_data2[[pair[2]]]
  
  # Perform t-test for pit diameter
  ttest <- t.test(parasite_data$PitDiameter, host_data$PitDiameter, var.equal = FALSE)
  pitDiameter_results <- rbind(pitDiameter_results, data.frame(
    Parasite = pair[1],
    Host = pair[2],
    ParasiteMean = ttest$estimate[1],
    HostMean = ttest$estimate[2],
    MeanDifference = diff(ttest$estimate),
    pvalue = ttest$p.value,
    stringsAsFactors = FALSE
  ))
}



# Perform t-tests for Pit Opening Diameter
for (pair in list(c("P. robustus", "V. thyrsoidea"), 
                  c("P. perrotettii", "T. guianensis"), 
                  c("S. rhynchophyllus", "T. tipu"), 
                  c("V. album", "P. nigra"))) {
  
  parasite_data <- species_data2[[pair[1]]]
  host_data <- species_data2[[pair[2]]]
  
  # Perform t-test for pit Opening diameter
  ttest <- t.test(parasite_data$PitOpening, host_data$PitOpening, var.equal = FALSE)
  pitOpening_results <- rbind(pitOpening_results, data.frame(
    Parasite = pair[1],
    Host = pair[2],
    ParasiteMean = ttest$estimate[1],
    HostMean = ttest$estimate[2],
    MeanDifference = diff(ttest$estimate),
    pvalue = ttest$p.value,
    stringsAsFactors = FALSE
  ))
}
fwrite(pitDiameter_results, file = here("outputs", "tables", "PitDiameter_results.csv"))
fwrite(pitOpening_results, file = here("outputs", "tables", "PitOpening_results.csv"))