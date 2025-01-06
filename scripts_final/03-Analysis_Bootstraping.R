library(here)
library(tidyverse)
source(here("scripts", "Functions.R")) 

Median_data<- read.csv(here("data", "processed", "Median_data.csv"))
Mean_data<- read.csv(here("data", "processed", "Mean_data.csv"))
# Factor level reordering across all data frames
relevel_factors(ls())

# Define species pairs for comparison
species_pairs <- list(
  c("Psittacanthus robustus", "Vochysia thyrsoidea"),
  c("Phoradendron perrotettii", "Tapirira guianensis"),
  c("Struthanthus rhynchophyllus", "Tipuana tipu"),
  c("Viscum album", "Populus nigra")
)
Var <- colnames(Mean_data)[sapply(Mean_data, is.numeric)]
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
                        
                        Median_data %>%
                          group_by(ssp) %>%
                          summarise(across(everything()[-c(1,11)], ~ mean(.x, na.rm = TRUE)))%>%
                          rename(Group = ssp)
)


######################################################################################
##Parasite x Host
#####################################################################################
Mean_obs_diff <- Obs_means[Obs_means$Group == "Parasite", ][-1]- Obs_means[Obs_means$Group == "Host", ][-1]
Median_obs_diff <- Obs_medians[Obs_medians$Group == "Parasite", ][-1]- Obs_medians[Obs_medians$Group == "Host", ][-1]
mean_boot <- data.table::fread(here("data", "processed", "ressampled","Means_parasitism.csv"))[,2:10]
median_boot <- data.table::fread(here("data", "processed", "ressampled","Medians_parasitism.csv"))[,2:10]
mean_boot[1,] <- Mean_obs_diff
median_boot[1,] <- Median_obs_diff


PxH_mean_pvalues <- numeric(length(Var))  # Initialize a named numeric vector for p-values
names(PxH_mean_pvalues) <- Var


PxH_mean_pvalues <- sapply(Var, function(v) {
  data <- mean_boot[[v]]
  plot(density(data), main =paste("Mean", v, "PxH", sep = " "))
  abline(v =data[1], col = "red", lwd = 2)
  abline(v = quantile(data, c(0.025, 0.975)), col =  "black", lwd = 2)
  bootp <- sum(abs(data)>=abs(data[1]))/length(data)
  text(x = mean(data), y = max(density(data)$y) * 0.9, labels = paste("p =", bootp))
  bootp
})


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

mean_effect <- Mean_obs_diff/Obs_means[Obs_means$Group == "Parasite", ][-1]
median_effect <- Median_obs_diff/Obs_medians[Obs_medians$Group == "Parasite", ][-1]
# Add a column with row names and convert to a tibble
PxH_results <- rbind(
  Mean_obs_diff, 
  PxH_mean_pvalues, 
  mean_effect * 100, 
  Median_obs_diff, 
  PxH_median_pvalues, 
  median_effect * 100
) %>% 
  tibble() %>% 
  mutate(Label = c("Mean diff", "p-value", "Effect size", "Median diff", "p-value", "Effect size")) %>% 
  select(Label, everything())  # Reorder to make 'Label' the first column


##################################################
##Psittacanthus x Vochysuia
Mean_obs_diff <- Obs_means[Obs_means$Group == "Psittacanthus robustus", ][-1]- Obs_means[Obs_means$Group == "Vochysia thyrsoidea", ][-1]
Median_obs_diff <- Obs_medians[Obs_medians$Group == "Psittacanthus robustus", ][-1]- Obs_medians[Obs_medians$Group == "Vochysia thyrsoidea", ][-1]
mean_boot <- data.table::fread(here("data", "processed", "ressampled","Means_Psittacanthus robustusXVochysia thyrsoidea.csv"))[,2:10]
median_boot <- data.table::fread(here("data", "processed", "ressampled","Medians_Psittacanthus robustusXVochysia thyrsoidea.csv"))[,2:10]
mean_boot[1,] <- Mean_obs_diff
median_boot[1,] <- Median_obs_diff


PrVt_mean_pvalues <- numeric(length(Var))# Initialize a named numeric vector for p-values
names(PrVt_mean_pvalues) <- Var

PrVt_mean_pvalues <- sapply(Var, function(v) {
  data <- mean_boot[[v]]
  plot(density(data), main = paste("Mean", v, "PrxVt", sep = " "))
  abline(v =data[1], col = "red", lwd = 2)
  abline(v = quantile(data, c(0.025, 0.975)), col =  "black", lwd = 2)
  bootp <- sum(abs(data)>=abs(data[1]))/length(data)
  text(x = mean(data), y = max(density(data)$y) * 0.9, labels = paste("p =", round(bootp, 3)))
  bootp
})

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

mean_effect <- Mean_obs_diff/Obs_means[Obs_means$Group == "Psittacanthus robustus", ][-1]
median_effect <- Median_obs_diff/Obs_medians[Obs_medians$Group == "Psittacanthus robustus", ][-1]
# Add a column with row names and convert to a tibble
PrxVt_results <- rbind(
  Mean_obs_diff, 
  PrVt_mean_pvalues, 
  mean_effect * 100, 
  Median_obs_diff, 
  PrVt_median_pvalues, 
  median_effect * 100
) %>% 
  tibble() %>% 
  mutate(Label = c("Mean diff", "p-value", "Effect size", "Median diff", "p-value", "Effect size")) %>% 
  select(Label, everything())  # Reorder to make 'Label' the first column



###
##Phoradendron
Mean_obs_diff <- Obs_means[Obs_means$Group == "Phoradendron perrotettii", ][-1]- Obs_means[Obs_means$Group == "Tapirira guianensis", ][-1]
Median_obs_diff <- Obs_medians[Obs_medians$Group == "Phoradendron perrotettii", ][-1]- Obs_medians[Obs_medians$Group == "Tapirira guianensis", ][-1]
mean_boot <- data.table::fread(here("data", "processed", "ressampled","Means_Phoradendron perrotettiiXTapirira guianensis.csv"))[,2:10]
median_boot <- data.table::fread(here("data", "processed", "ressampled","Medians_Phoradendron perrotettiiXTapirira guianensis.csv"))[,2:10]
mean_boot[1,] <- Mean_obs_diff
median_boot[1,] <- Median_obs_diff


PpTg_mean_pvalues <- numeric(length(Var)) # Initialize a named numeric vector for p-values
names(PpTg_mean_pvalues) <- Var

PpTg_mean_pvalues <- sapply(Var, function(v) {
  data <- mean_boot[[v]]
  plot(density(data), main =paste("Mean", v, "PpxTg", sep = " "))
  abline(v =data[1], col = "red", lwd = 2)
  abline(v = quantile(data, c(0.025, 0.975)), col =  "black", lwd = 2)
  bootp <- sum(abs(data)>=abs(data[1]))/length(data)
  text(x = mean(data), y = max(density(data)$y) * 0.9, labels = paste("p =", round(bootp, 3)))
  bootp
})

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

mean_effect <- Mean_obs_diff/Obs_means[Obs_means$Group == "Phoradendron perrotettii", ][-1]
median_effect <- Median_obs_diff/Obs_medians[Obs_medians$Group == "Phoradendron perrotettii", ][-1]
# Add a column with row names and convert to a tibble
PpxTg_results <- rbind(
  Mean_obs_diff, 
  PpTg_median_pvalues, 
  mean_effect * 100, 
  Median_obs_diff, 
  PpTg_median_pvalues, 
  median_effect * 100
) %>% 
  tibble() %>% 
  mutate(Label = c("Mean diff", "p-value", "Effect size", "Median diff", "p-value", "Effect size")) %>% 
  select(Label, everything())  # Reorder to make 'Label' the first column


###########################################
#struthanthus


Mean_obs_diff <- Obs_means[Obs_means$Group == "Struthanthus rhynchophyllus", ][-1]- Obs_means[Obs_means$Group == "Tipuana tipu", ][-1]
Median_obs_diff <- Obs_medians[Obs_medians$Group == "Struthanthus rhynchophyllus", ][-1]- Obs_medians[Obs_medians$Group == "Tipuana tipu", ][-1]
mean_boot <- data.table::fread(here("data", "processed", "ressampled","Means_Struthanthus rhynchophyllusXTipuana tipu.csv"))[,2:10]
median_boot <- data.table::fread(here("data", "processed", "ressampled","Medians_Struthanthus rhynchophyllusXTipuana tipu.csv"))[,2:10]
mean_boot[1,] <- Mean_obs_diff
median_boot[1,] <- Median_obs_diff


SrTt_mean_pvalues <- numeric(length(Var)) # Initialize a named numeric vector for p-values
names(SrTt_mean_pvalues) <- Var

SrTt_mean_pvalues <- sapply(Var, function(v) {
  data <- mean_boot[[v]]
  plot(density(data), main = paste("Mean", v, "SrxTt", sep = " "))
  abline(v =data[1], col = "red", lwd = 3)
  abline(v = quantile(data, c(0.025, 0.975)), col =  "black", lwd = 2)
  bootp <- sum(abs(data)>=abs(data[1]))/length(data)
  text(x = mean(data), y = max(density(data)$y) * 0.9, labels = paste("p =", round(bootp, 3)))
  bootp
})

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

mean_effect <- Mean_obs_diff/Obs_means[Obs_means$Group == "Struthanthus rhynchophyllus", ][-1]
median_effect <- Median_obs_diff/Obs_medians[Obs_medians$Group == "Struthanthus rhynchophyllus", ][-1]
# Add a column with row names and convert to a tibble
SrxTt_results <- rbind(
  Mean_obs_diff, 
  SrTt_mean_pvalues, 
  mean_effect * 100, 
  Median_obs_diff, 
  SrTt_median_pvalues, 
  median_effect * 100
) %>% 
  tibble() %>% 
  mutate(Label = c("Mean diff", "p-value", "Effect size", "Median diff", "p-value", "Effect size")) %>% 
  select(Label, everything())  # Reorder to make 'Label' the first column

#######################################
##viscum

Mean_obs_diff <- Obs_means[Obs_means$Group == "Viscum album", ][-1]- Obs_means[Obs_means$Group == "Populus nigra", ][-1]
Median_obs_diff <- Obs_medians[Obs_medians$Group == "Viscum album", ][-1]- Obs_medians[Obs_medians$Group == "Populus nigra", ][-1]
mean_boot <- data.table::fread(here("data", "processed", "ressampled","Means_Viscum albumXPopulus nigra.csv"))[,2:10]
median_boot <- data.table::fread(here("data", "processed", "ressampled","Medians_Viscum albumXPopulus nigra.csv"))[,2:10]
mean_boot[1,] <- Mean_obs_diff
median_boot[1,] <- Median_obs_diff


VaPn_mean_pvalues <- numeric(length(Var)) # Initialize a named numeric vector for p-values
names(VaPn_mean_pvalues) <- Var

VaPn_mean_pvalues <- sapply(Var, function(v) {
  data <- mean_boot[[v]]
  plot(density(data), main = paste("Mean", v, "VaxPn", sep = " "))
  abline(v =data[1], col = "red", lwd = 2)
  abline(v = quantile(data, c(0.025, 0.975)), col =  "black", lwd = 2)
  bootp <- sum(abs(data)>=abs(data[1]))/length(data)
  text(x = mean(data), y = max(density(data)$y) * 0.9, labels = paste("p =", round(bootp, 3)))
  bootp
})

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

mean_effect <- Mean_obs_diff/Obs_means[Obs_means$Group == "Viscum album", ][-1]
median_effect <- Median_obs_diff/Obs_medians[Obs_medians$Group == "Viscum album", ][-1]
# Add a column with row names and convert to a tibble
VaxPn_results <- rbind(
  Mean_obs_diff, 
  VaPn_mean_pvalues, 
  mean_effect * 100, 
  Median_obs_diff, 
  VaPn_median_pvalues, 
  median_effect * 100
) %>% 
  tibble() %>% 
  mutate(Label = c("Mean diff", "p-value", "Effect size", "Median diff", "p-value", "Effect size")) %>% 
  select(Label, everything())  # Reorder to make 'Label' the first column

