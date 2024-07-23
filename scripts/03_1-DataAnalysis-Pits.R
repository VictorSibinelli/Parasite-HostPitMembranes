######################################################################
#
#Victor Sibinelli (victor.sibinelli@usp.br / sibinelli95@gmail.com)
#13/07/2024
# Script 03.1 - Data Analysis - Pit membranes
#################################################################
source(here("scripts","02-TestAssumptions.R"))
rm(list=ls())

#load data
pitdata_clean <- read.csv(here("data","processed","pitdata_clean.csv"))


# List of data frame names
dataframes <- ls()

# Relevel the factors for each data frame
for (df_name in dataframes) {
  df <- get(df_name)  # Get the data frame by name
  
  if ("ssp" %in% colnames(df)) {  # Check if 'ssp' column exists
    df$ssp <- factor(df$ssp, levels = c("Psittacanthus robustus", "Vochysia thyrsoidea",
                                        "Phoradendeon perrotettii", "Tapirira guianensis",
                                        "Struthanthus rhynchophyllus", "Tipuana tipu",
                                        "Viscum album", "Populus nigra"))
    assign(df_name, df)  # Assign the modified data frame back to its name
  }
  rm(df,df_name) #remove duplicated dataframe
}


###testing difference between pe and pc

# Loop through each unique species in the dataset
for (ssp in unique(pitdata_clean$ssp)) {
  # Subset data for the current species
  data <- pitdata_clean[pitdata_clean$ssp == ssp, ]
  
  # Perform Welch's t-test comparing pcavg and peavg
  test_result <- t.test(data$pcavg, data$peavg, var.equal = FALSE)
  
  # Print the species and the test result
  cat("\nSpecies:", ssp, "\n")
  print(test_result)
}

#####All samples showed slightly thicker centers


##sub seting species
psittacanthus_data <-pitdata_clean[pitdata_clean$ssp == "Psittacanthus robustus", ]
vochysia_data <-pitdata_clean[pitdata_clean$ssp == "Vochysia thyrsoidea", ]
phoradendron_data <-pitdata_clean[pitdata_clean$ssp == "Phoradendeon perrotettii", ]
tapirira_data <-pitdata_clean[pitdata_clean$ssp == "Tapirira guianensis", ]
struthanthus_data <-pitdata_clean[pitdata_clean == "Struthanthus rhynchophyllus", ]
tipuana_data <-pitdata_clean[pitdata_clean$ssp == "Tipuana tipu", ]
viscum_data <-pitdata_clean[pitdata_clean == "Viscum album", ]
populus_data <-pitdata_clean[pitdata_clean == "Populus nigra", ]


#Welch's test for Pit membrane thickness
t.test(pitdata_clean$pitavg[pitdata_clean$parasitism=="p"], pitdata_clean$pitavg[pitdata_clean$parasitism=="h"], var.equal = FALSE)

# Perform t-tests and store the results in a data frame
pit_mean_diff <- data.frame(Parasite = character(),
                            ParasiteMean = numeric(),
                            Host = character(),
                            HostMean = numeric(),
                            MeanDifference = numeric(),
                            stringsAsFactors = FALSE)

# Perform t-tests and extract the difference of means
ttest1 <- t.test(psittacanthus_data$pitavg, vochysia_data$pitavg, var.equal = FALSE)
mean_diff1 <- c(ttest1$estimate, abs(diff(ttest1$estimate)))

ttest2 <- t.test(phoradendron_data$pitavg, tapirira_data$pitavg, var.equal = FALSE)
mean_diff2 <- c(ttest2$estimate, abs(diff(ttest2$estimate)))

ttest3 <- t.test(struthanthus_data$pitavg, tipuana_data$pitavg, var.equal = FALSE)
mean_diff3 <- c(ttest3$estimate, abs(diff(ttest3$estimate)))

ttest4 <- t.test(viscum_data$pitavg, populus_data$pitavg, var.equal = FALSE)
mean_diff4 <- c(ttest4$estimate, abs(diff(ttest4$estimate)))

# Add the results to the data frame
pit_mean_diff <- rbind(pit_mean_diff, data.frame(Parasite = "P. robustus", ParasiteMean = mean_diff1[1], Host = "V. thyrsoidea", HostMean = mean_diff1[2], MeanDifference = mean_diff1[3]))
pit_mean_diff <- rbind(pit_mean_diff, data.frame(Parasite = "P. perrotettii", ParasiteMean = mean_diff2[1], Host = "T. guianensis", HostMean = mean_diff2[2], MeanDifference = mean_diff2[3]))
pit_mean_diff <- rbind(pit_mean_diff, data.frame(Parasite = "S. rhynchophyllus", ParasiteMean = mean_diff3[1], Host = "T. tipu", HostMean = mean_diff3[2], MeanDifference = mean_diff3[3]))
pit_mean_diff <- rbind(pit_mean_diff, data.frame(Parasite = "V. album", ParasiteMean = mean_diff4[1], Host = "P. nigra", HostMean = mean_diff4[2], MeanDifference = mean_diff4[3]))

pit_means <- tapply(pitdata_clean$pitavg, pitdata_clean$ssp, mean, na.rm=T)
pcd_measn <- tapply(pitdata_clean$pcd, pitdata_clean$ssp, mean, na.rm=T)
# Print the results

print(pit_mean_diff)
fwrite(pit_mean_diff,file=here("outputs", "tables","pit_membrane_diff.csv"))


#Welch's test for Pit chamber depth

t.test(psittacanthus_data$pcd, vochysia_data$pcd, var.equal = FALSE) #highly significant 

t.test(phoradendron_data$pcd, tapirira_data$pcd, var.equal = FALSE) #highly significant 

t.test(struthanthus_data$pcd, tipuana_data$pcd, var.equal = FALSE) #highly significant 

t.test(viscum_data$pcd, populus_data$pcd, var.equal = FALSE) #highly significant 

# Perform t-tests and store the results in a data frame
pcd_results <- data.frame(Parasite = character(),
                            ParasiteMean = numeric(),
                            Host = character(),
                            HostMean = numeric(),
                            MeanDifference = numeric(),
                            stringsAsFactors = FALSE)


# Perform t-tests and extract the difference of means
ttest1 <- t.test(psittacanthus_data$pcd, vochysia_data$pcd, var.equal = FALSE)
mean_diff1 <- c(ttest1$estimate, abs(diff(ttest1$estimate)))

ttest2 <- t.test(phoradendron_data$pcd, tapirira_data$pcd, var.equal = FALSE)
mean_diff2 <- c(ttest2$estimate, abs(diff(ttest2$estimate)))

ttest3 <- t.test(struthanthus_data$pcd, tipuana_data$pcd, var.equal = FALSE)
mean_diff3 <- c(ttest3$estimate, abs(diff(ttest3$estimate)))

ttest4 <- t.test(viscum_data$pcd, populus_data$pcd, var.equal = FALSE)
mean_diff4 <- c(ttest4$estimate, abs(diff(ttest4$estimate)))

# Add the results to the data frame
pcd_results <- rbind(pcd_results, data.frame(Parasite = "P. robustus", ParasiteMean = mean_diff1[1], Host = "V. thyrsoidea", HostMean = mean_diff1[2], MeanDifference = mean_diff1[3]))
pcd_results <- rbind(pcd_results, data.frame(Parasite = "P. perrotettii", ParasiteMean = mean_diff2[1], Host = "T. guianensis", HostMean = mean_diff2[2], MeanDifference = mean_diff2[3]))
pcd_results <- rbind(pcd_results, data.frame(Parasite = "S. rhynchophyllus", ParasiteMean = mean_diff3[1], Host = "T. tipu", HostMean = mean_diff3[2], MeanDifference = mean_diff3[3]))
pcd_results <- rbind(pcd_results, data.frame(Parasite = "V. album", ParasiteMean = mean_diff4[1], Host = "P. nigra", HostMean = mean_diff4[2], MeanDifference = mean_diff4[3]))


# Print the results
print(pcd_results)
fwrite(pcd_results,file=here("outputs", "tables","pcd_results.csv"))
