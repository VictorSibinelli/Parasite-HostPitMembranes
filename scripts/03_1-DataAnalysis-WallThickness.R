######################################################################
#
# Victor Sibinelli (victor.sibinelli@usp.br / sibinelli95@gmail.com)
# 13/07/2024
# Script 03.1 - Data Analysis - Vessel walls
#################################################################
library(here)
source(here("scripts", "02_1-TestAssumptions-WallThickness.R"))
rm(list=ls())

wdata_clean <- read.csv(here("data", "processed", "wdata_clean.csv"))



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
# Data subsets
species_data <- list(
  "P. robustus" = wdata_clean[wdata_clean$ssp == "Psittacanthus robustus", ],
  "V. thyrsoidea" = wdata_clean[wdata_clean$ssp == "Vochysia thyrsoidea", ],
  "P. perrotettii" = wdata_clean[wdata_clean$ssp == "Phoradendron perrotettii", ],
  "T. guianensis" = wdata_clean[wdata_clean$ssp == "Tapirira guianensis", ],
  "S. rhynchophyllus" = wdata_clean[wdata_clean$ssp == "Struthanthus rhynchophyllus", ],
  "T. tipu" = wdata_clean[wdata_clean$ssp == "Tipuana tipu", ],
  "V. album" = wdata_clean[wdata_clean$ssp == "Viscum album", ],
  "P. nigra" = wdata_clean[wdata_clean$ssp == "Populus nigra", ])

VWall_results <- data.frame(Parasite = character(),
                            Host = character(),
                            ParasiteMean = numeric(),
                            HostMean = numeric(),
                            MeanDifference = numeric(),
                            pvalue = numeric(),
                            stringsAsFactors = FALSE)

# Perform t-tests for Pit membrane thickness
for (pair in list(c("P. robustus", "V. thyrsoidea"), 
                  c("P. perrotettii", "T. guianensis"), 
                  c("S. rhynchophyllus", "T. tipu"), 
                  c("V. album", "P. nigra"))){
  
  parasite_data <- species_data[[pair[1]]]
  host_data <- species_data[[pair[2]]]
  
  # Perform t-test for vessel wall thickness
  ttest <- t.test(parasite_data$wthickness, host_data$wthickness, var.equal = T)
  VWall_results <- rbind(VWall_results, data.frame(
    Parasite = pair[1],
    Host = pair[2],
    ParasiteMean = ttest$estimate[1],
    HostMean = ttest$estimate[2],
    MeanDifference = diff(ttest$estimate),
    pvalue= ttest$p.value,
    stringsAsFactors = FALSE
  ))
}

fwrite(VWall_results, file = here("outputs", "tables", "VWall_results.csv"))

fwrite

  