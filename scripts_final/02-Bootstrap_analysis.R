library(here)
library(tidyverse)
source(here("scripts", "Functions.R")) 
Median_data<- read.csv(here("data", "processed", "Median_data.csv"))
TEM_data <- read.csv(here("data", "processed", "PitMembrane_data.csv"))
relevel_factors(ls())

# Define species pairs for comparison
species_pairs <- list(
  c("Psittacanthus robustus", "Vochysia thyrsoidea"),
  c("Phoradendron perrotettii", "Tapirira guianensis"),
  c("Struthanthus rhynchophyllus", "Tipuana tipu"),
  c("Viscum album", "Populus nigra")
)

Obs_medians <-bind_rows(Median_data %>%
                          group_by(parasitism) %>%
                          summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE))) %>%
                          rename(Group = parasitism),
                        
                        Median_data %>%
                          group_by(ssp) %>%
                          summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)))%>%
                          rename(Group = ssp)) %>% 
  cbind( bind_rows(TEM_data %>% 
              group_by(parasitism) %>% 
                summarise(across(where(is.numeric), ~ median(.x, na.rm = TRUE))) %>%
              rename(Group = parasitism),
            
            TEM_data %>% 
              group_by(ssp) %>% 
              summarise(across(where(is.numeric), ~ median(.x, na.rm = TRUE))) %>%
              rename(Group = ssp))[,-1]
)

Var <- colnames(Obs_medians)[sapply(Obs_medians, is.numeric)]
##Parasite x Host
#####################################################################################

Median_obs_diff <- Obs_medians[Obs_medians$Group == "Parasite", ][-1]- Obs_medians[Obs_medians$Group == "Host", ][-1]

median_boot <- data.table::fread(here("data", "processed", "ressampled","Medians_parasitism.csv"))[,-1]

median_boot[1,] <- Median_obs_diff






PxH_median_pvalues <- numeric(length(Var)) # Initialize a named numeric vector for p-values
names(PxH_median_pvalues) <- Var

PxH_median_pvalues <- sapply(Var, function(v) {
  data <- median_boot[[v]]
  plot(density(data), main = paste("Median", v, "PxH", sep = " "))
  abline(v =data[1], col = "red", lwd = 2)
  abline(v = quantile(data, c(0.025, 0.975)), col =  "black", lwd = 2)
  bootp <- sum(abs(data)>=abs(data[1]))/length(data)
  text(x = mean(data), y = max(density(data)$y) * 0.9, labels = paste("p =", bootp))
  bootp
})


median_effect <- Median_obs_diff/Obs_medians[Obs_medians$Group == "Host", ][-1]
# Add a column with row names and convert to a tibble
PxH_results <- rbind( 
  Median_obs_diff, 
  PxH_median_pvalues, 
  median_effect * 100
) %>% 
  tibble() %>% 
  mutate(Label = c("Median diff", "p-value", "Effect size")) %>% 
  select(Label, everything())  # Reorder to make 'Label' the first column






##################################################
##Psittacanthus x Vochysuia

Median_obs_diff <- Obs_medians[Obs_medians$Group == "Psittacanthus robustus", ][-1]- Obs_medians[Obs_medians$Group == "Vochysia thyrsoidea", ][-1]

median_boot <- data.table::fread(here("data", "processed", "ressampled","Medians_Psittacanthus robustusXVochysia thyrsoidea.csv"))[,-1]

median_boot[1,] <- Median_obs_diff




PrVt_median_pvalues <- numeric(length(Var))# Initialize a named numeric vector for p-values
names(PrVt_median_pvalues) <- Var

PrVt_median_pvalues <- sapply(Var, function(v) {
  data <- median_boot[[v]]
  plot(density(data), main = paste("Median", v, "PrxVt", sep = " "))
  abline(v =data[1], col = "red", lwd = 2)
  abline(v = quantile(data, c(0.025, 0.975)), col =  "black", lwd = 2)
  bootp <- sum(abs(data)>=abs(data[1]))/length(data)
  text(x = mean(data), y = max(density(data)$y) * 0.9, labels = paste("p =", round(bootp, 3)))
  bootp
})


median_effect <- Median_obs_diff/Obs_medians[Obs_medians$Group == "Vochysia thyrsoidea", ][-1]
# Add a column with row names and convert to a tibble
PrxVt_results <- rbind(
  Median_obs_diff, 
  PrVt_median_pvalues, 
  median_effect * 100
) %>% 
  tibble() %>% 
  mutate(Label = c("Median diff", "p-value", "Effect size")) %>% 
  select(Label, everything())  # Reorder to make 'Label' the first column



###
##Phoradendron

Median_obs_diff <- Obs_medians[Obs_medians$Group == "Phoradendron perrotettii", ][-1]- Obs_medians[Obs_medians$Group == "Tapirira guianensis", ][-1]

median_boot <- data.table::fread(here("data", "processed", "ressampled","Medians_Phoradendron perrotettiiXTapirira guianensis.csv"))[,-1]

median_boot[1,] <- Median_obs_diff

PpTg_median_pvalues <- numeric(length(Var)) %>% as.tibble()  # Initialize a named numeric vector for p-values
names(PpTg_median_pvalues) <- Var

PpTg_median_pvalues <- sapply(Var, function(v) {
  data <- median_boot[[v]]
  plot(density(data), main = paste("Median", v, "PpxTg", sep = " "))
  abline(v =data[1], col = "red", lwd = 2)
  abline(v = quantile(data, c(0.025, 0.975)), col =  "black", lwd = 2)
  bootp <- sum(abs(data)>=abs(data[1]))/length(data)
  text(x = mean(data), y = max(density(data)$y) * 0.9, labels = paste("p =", round(bootp, 3)))
  bootp
})


median_effect <- Median_obs_diff/Obs_medians[Obs_medians$Group == "Tapirira guianensis", ][-1]
# Add a column with row names and convert to a tibble
PpxTg_results <- rbind(
  Median_obs_diff, 
  PpTg_median_pvalues, 
  median_effect * 100
) %>% 
  tibble() %>% 
  mutate(Label = c("Median diff", "p-value", "Effect size")) %>% 
  select(Label, everything())  # Reorder to make 'Label' the first column


###########################################
#struthanthus



Median_obs_diff <- Obs_medians[Obs_medians$Group == "Struthanthus rhynchophyllus", ][-1]- Obs_medians[Obs_medians$Group == "Tipuana tipu", ][-1]

median_boot <- data.table::fread(here("data", "processed", "ressampled","Medians_Struthanthus rhynchophyllusXTipuana tipu.csv"))[,-1]

median_boot[1,] <- Median_obs_diff



SrTt_median_pvalues <- numeric(length(Var)) # Initialize a named numeric vector for p-values
names(SrTt_median_pvalues) <- Var

SrTt_median_pvalues <- sapply(Var, function(v) {
  data <- median_boot[[v]]
  plot(density(data), main =paste("Median", v, "SrxTt", sep = " "))
  abline(v =data[1], col = "red", lwd = 2)
  abline(v = quantile(data, c(0.025, 0.975)), col =  "black", lwd = 2)
  bootp <- sum(abs(data)>=abs(data[1]))/length(data)

  text(x = mean(data), y = max(density(data)$y) * 0.9, labels = paste("p =", round(bootp, 3)))
  bootp
})

median_effect <- Median_obs_diff/Obs_medians[Obs_medians$Group == "Tipuana tipu", ][-1]
# Add a column with row names and convert to a tibble
SrxTt_results <- rbind(
  Median_obs_diff, 
  SrTt_median_pvalues, 
  median_effect * 100
) %>% 
  tibble() %>% 
  mutate(Label = c("Median diff", "p-value", "Effect size")) %>% 
  select(Label, everything())  # Reorder to make 'Label' the first column

#######################################
##viscum


Median_obs_diff <- Obs_medians[Obs_medians$Group == "Viscum album", ][-1]- Obs_medians[Obs_medians$Group == "Populus nigra", ][-1]

median_boot <- data.table::fread(here("data", "processed", "ressampled","Medians_Viscum albumXPopulus nigra.csv"))[,-1]

median_boot[1,] <- Median_obs_diff



VaPn_median_pvalues <- numeric(length(Var))  # Initialize a named numeric vector for p-values
names(PpTg_median_pvalues) <- Var

VaPn_median_pvalues <- sapply(Var, function(v) {
  data <- median_boot[[v]]
  plot(density(data), main =paste("Median", v, "VaxPn", sep = " "))
  abline(v =data[1], col = "red", lwd = 2)
  abline(v = quantile(data, c(0.025, 0.975)), col =  "black", lwd = 2)
  bootp <- sum(abs(data)>=abs(data[1]))/length(data)
  text(x = mean(data), y = max(density(data)$y) * 0.9, labels = paste("p =", round(bootp, 3)))
  bootp
})


median_effect <- Median_obs_diff/Obs_medians[Obs_medians$Group == "Populus nigra", ][-1]
# Add a column with row names and convert to a tibble
VaxPn_results <- rbind(
  Median_obs_diff, 
  VaPn_median_pvalues, 
  median_effect * 100
) %>% 
  tibble() %>% 
  mutate(Label = c("Median diff", "p-value", "Effect size")) %>% 
  select(Label, everything())  # Reorder to make 'Label' the first column













# Get the names of all objects containing "_results"
results_names <- grep(pattern = "_results", ls(), value = TRUE)

# Loop over the names and save each as a CSV file
for (name in results_names) {
  # Get the object dynamically using `get`
  obj <- get(name)
  
  # Ensure the object is a data frame before saving
  if (is.data.frame(obj)) {
    # Create a file name based on the object name
    file_name <- paste0(name, ".csv")
    
    # Save the data frame as a CSV file
    write.csv(obj, file = here("outputs","tables","Bootstrap",paste0("Bootstrap",file_name,collapse = "_")), row.names = FALSE)
  }
  print(paste(file_name, "saved", sep=" "))
}

