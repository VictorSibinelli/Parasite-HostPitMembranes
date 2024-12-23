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
                        
                        Mean_data %>%
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


PxH_mean_pvalues <- numeric(length(Var)) %>% as.tibble()  # Initialize a named numeric vector for p-values
names(PxH_mean_pvalues) <- Var

PxH_mean_pvalues <- sapply(Var, function(v) {
  data <- mean_boot[[v]]
  plot(density(data), main = v)
  abline(v =data[1], col = "red", lwd = 2)
  abline(v = quantile(data, c(0.025, 0.975)), col =  "black", lwd = 2)
  bootp <- mean(data >= abs(data[1]))
  text(x = mean(data), y = max(density(data)$y) * 0.9, labels = paste("p =", round(bootp, 3)))
  bootp
})


PxH_median_pvalues <- numeric(length(Var)) %>% as.tibble()  # Initialize a named numeric vector for p-values
names(PxH_median_pvalues) <- Var

PxH_mean_pvalues <- sapply(Var, function(v) {
  data <- median_boot[[v]]
  plot(density(data), main = v)
  abline(v =data[1], col = "red", lwd = 2)
  abline(v = quantile(data, c(0.025, 0.975)), col =  "black", lwd = 2)
  bootp <- mean(data >= abs(data[1]))
  text(x = mean(data), y = max(density(data)$y) * 0.9, labels = paste("p =", round(bootp, 3)))
  bootp
})


##################################################
##Psittacanthus x Vochysuia
Mean_obs_diff <- Obs_means[Obs_means$Group == "Psittacanthus robustus", ][-1]- Obs_means[Obs_means$Group == "Vochysia thyrsoidea", ][-1]
Median_obs_diff <- Obs_medians[Obs_medians$Group == "Psittacanthus robustus", ][-1]- Obs_medians[Obs_medians$Group == "Vochysia thyrsoidea", ][-1]
mean_boot <- data.table::fread(here("data", "processed", "ressampled","Means_Psittacanthus robustusXVochysia thyrsoidea.csv"))[,2:10]
median_boot <- data.table::fread(here("data", "processed", "ressampled","Medians_Psittacanthus robustusXVochysia thyrsoidea.csv"))[,2:10]
mean_boot[1,] <- Mean_obs_diff
median_boot[1,] <- Median_obs_diff


PrVt_mean_pvalues <- numeric(length(Var)) %>% as.tibble()  # Initialize a named numeric vector for p-values
names(PrVt_mean_pvalues) <- Var

PrVt_mean_pvalues <- sapply(Var, function(v) {
  data <- mean_boot[[v]]
  plot(density(data), main = v)
  abline(v =data[1], col = "red", lwd = 2)
  abline(v = quantile(data, c(0.025, 0.975)), col =  "black", lwd = 2)
  bootp <- mean(data >= abs(data[1]))
  text(x = mean(data), y = max(density(data)$y) * 0.9, labels = paste("p =", round(bootp, 3)))
  bootp
})

PrVt_median_pvalues <- numeric(length(Var)) %>% as.tibble()  # Initialize a named numeric vector for p-values
names(PrVt_median_pvalues) <- Var

PrVt_median_pvalues <- sapply(Var, function(v) {
  data <- median_boot[[v]]
  plot(density(data), main = v)
  abline(v =data[1], col = "red", lwd = 2)
  abline(v = quantile(data, c(0.025, 0.975)), col =  "black", lwd = 2)
  bootp <- mean(data >= abs(data[1]))
  text(x = mean(data), y = max(density(data)$y) * 0.9, labels = paste("p =", round(bootp, 3)))
  bootp
})


###
##Phoradendron
Mean_obs_diff <- Obs_means[Obs_means$Group == "Phoradendron perrotettii", ][-1]- Obs_means[Obs_means$Group == "Tapirira guianensis", ][-1]
Median_obs_diff <- Obs_medians[Obs_medians$Group == "Phoradendron perrotettii", ][-1]- Obs_medians[Obs_medians$Group == "Tapirira guianensis", ][-1]
mean_boot <- data.table::fread(here("data", "processed", "ressampled","Means_Phoradendron perrotettiiXTapirira guianensis.csv"))[,2:10]
median_boot <- data.table::fread(here("data", "processed", "ressampled","Medians_Phoradendron perrotettiiXTapirira guianensis.csv"))[,2:10]
mean_boot[1,] <- Mean_obs_diff
median_boot[1,] <- Median_obs_diff


PpTg_mean_pvalues <- numeric(length(Var)) %>% as.tibble()  # Initialize a named numeric vector for p-values
names(PpTg_mean_pvalues) <- Var

PpTg_mean_pvalues <- sapply(Var, function(v) {
  data <- mean_boot[[v]]
  plot(density(data), main = v)
  abline(v =data[1], col = "red", lwd = 2)
  abline(v = quantile(data, c(0.025, 0.975)), col =  "black", lwd = 2)
  bootp <- mean(data >= abs(data[1]))
  text(x = mean(data), y = max(density(data)$y) * 0.9, labels = paste("p =", round(bootp, 3)))
  bootp
})

PpTg_median_pvalues <- numeric(length(Var)) %>% as.tibble()  # Initialize a named numeric vector for p-values
names(PpTg_median_pvalues) <- Var

PpTg_median_pvalues <- sapply(Var, function(v) {
  data <- median_boot[[v]]
  plot(density(data), main = v)
  abline(v =data[1], col = "red", lwd = 2)
  abline(v = quantile(data, c(0.025, 0.975)), col =  "black", lwd = 2)
  bootp <- mean(data >= abs(data[1]))
  text(x = mean(data), y = max(density(data)$y) * 0.9, labels = paste("p =", round(bootp, 3)))
  bootp
})

