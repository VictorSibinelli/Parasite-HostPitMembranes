#########
#######Tpm and Pcd

library(tidyverse)
library(here)
source(here("scripts","00-library.R"))
source(here("scripts", "Functions.R"))
PitMembrane_data<- read.csv(here("data", "processed", "PitMembrane_data.csv"))
species_pairs <- list(
  c("Psittacanthus robustus", "Vochysia thyrsoidea"),
  c("Phoradendron perrotettii", "Tapirira guianensis"),
  c("Struthanthus rhynchophyllus", "Tipuana tipu"),
  c("Viscum album", "Populus nigra")
)

relevel_factors(ls())

# Conduct a t-test for parasitism and store the p-value
Pcd_PxH <- t.test(pcd ~ parasitism, data = PitMembrane_data, var.equal = FALSE, na.rm = TRUE)
Pcd_p_values <- data.frame(Model = "PxH", p_value = Pcd_PxH$p.value)

# Prepare initial results dataframe
Pcd_results <- PitMembrane_data %>% 
  group_by(Group = parasitism) %>% 
  summarise(
    mean = mean(pcd, na.rm = TRUE),
    sd = sd(pcd, na.rm = TRUE),
    n = n(),
    ci_lower = mean - qt(0.975, df = n() - 1) * (sd / sqrt(n())),
    ci_upper = mean + qt(0.975, df = n() - 1) * (sd / sqrt(n()))
  )

# Loop through species pairs and perform pairwise t-tests
for (pair in species_pairs) {
  # Subset data for the current pair
  subset_data <- subset(PitMembrane_data, ssp %in% pair)
  
  # Perform t-test
  t <- t.test(pcd ~ ssp, data = subset_data, var.equal = FALSE, na.rm = TRUE)
  
  # Print t-test results
  print(pair)
  print(t)
  
  # Calculate group-level statistics and confidence intervals
  res <- subset_data %>% 
    group_by(Group = ssp) %>% 
    summarise(
      mean = mean(pcd, na.rm = TRUE),
      sd = sd(pcd, na.rm = TRUE),
      n = n(),
      ci_lower = mean - qt(0.975, df = n() - 1) * (sd / sqrt(n())),
      ci_upper = mean + qt(0.975, df = n() - 1) * (sd / sqrt(n()))
    )
  
  # Combine results into the main results dataframe
  Pcd_results <- bind_rows(Pcd_results, res)
  
  # Append p-value to the p-values dataframe
  Pcd_p_values <- bind_rows(
    Pcd_p_values,
    data.frame(Model = paste(pair, collapse = " x "), p_value = t$p.value)
  )
}

# Display results
print("Confidence Interval and Group Summary Results:")
print(Pcd_results)

print("P-Values:")
print(Pcd_p_values)

# Conduct a t-test for parasitism and store the p-value
Tpm_PxH <- t.test(Tpm ~ parasitism, data = PitMembrane_data, var.equal = FALSE, na.rm = TRUE)
Tpm_p_values <- data.frame(Model = "PxH", p_value = Tpm_PxH$p.value)

# Prepare initial results dataframe
Tpm_results <- PitMembrane_data %>% 
  group_by(Group = parasitism) %>% 
  summarise(
    mean = mean(Tpm, na.rm = TRUE),
    sd = sd(Tpm, na.rm = TRUE),
    n = n(),
    ci_lower = mean - qt(0.975, df = n() - 1) * (sd / sqrt(n())),
    ci_upper = mean + qt(0.975, df = n() - 1) * (sd / sqrt(n()))
  )

# Loop through species pairs and perform pairwise t-tests
for (pair in species_pairs) {
  # Subset data for the current pair
  subset_data <- subset(PitMembrane_data, ssp %in% pair)
  
  # Perform t-test
  t <- t.test(Tpm ~ ssp, data = subset_data, var.equal = FALSE, na.rm = TRUE)
  
  # Print t-test results
  print(pair)
  print(t)
  
  # Calculate group-level statistics and confidence intervals
  res <- subset_data %>% 
    group_by(Group = ssp) %>% 
    summarise(
      mean = mean(Tpm, na.rm = TRUE),
      sd = sd(Tpm, na.rm = TRUE),
      n = n(),
      ci_lower = mean - qt(0.975, df = n() - 1) * (sd / sqrt(n())),
      ci_upper = mean + qt(0.975, df = n() - 1) * (sd / sqrt(n()))
    )
  
  # Combine results into the main results dataframe
  Tpm_results <- bind_rows(Tpm_results, res)
  
  # Append p-value to the p-values dataframe
  Tpm_p_values <- bind_rows(
    Tpm_p_values,
    data.frame(Model = paste(pair, collapse = " x "), p_value = t$p.value)
  )
}

# Display results
print("Confidence Interval and Group Summary Results:")
print(Tpm_results)

print("P-Values:")
print(Tpm_p_values)

Pcd_table <- data.frame(
  Pair=character(),
  Parasite_mean=numeric(),
  Host_mean=numeric(),
  Relative_diff=numeric(),
  pvalue=numeric()
)

