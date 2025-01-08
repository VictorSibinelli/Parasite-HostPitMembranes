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
# Initialize the results table
Pcd_results <- data.frame(
  Model = character(),
  Parasite = character(),
  Host = character(),
  Relative_Difference = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

# 1. Add the Parasite vs Host comparison
Pcd_PxH <- t.test(pcd ~ parasitism, data = PitMembrane_data, var.equal = FALSE, na.rm = TRUE)

# Calculate group-level statistics for Parasite and Host
group_stats <- PitMembrane_data %>% 
  group_by(Group = parasitism) %>% 
  summarise(
    mean = mean(pcd, na.rm = TRUE),
    ci_lower = mean - qt(0.975, df = n() - 1) * (sd(pcd, na.rm = TRUE) / sqrt(n())),
    ci_upper = mean + qt(0.975, df = n() - 1) * (sd(pcd, na.rm = TRUE) / sqrt(n()))
  )

# Extract values for Parasite and Host
parasite <- group_stats[group_stats$Group == "Parasite", ]
host <- group_stats[group_stats$Group == "Host", ]

# Calculate Relative Difference
relative_diff <- (host$mean - parasite$mean) / parasite$mean

# Format the Parasite and Host columns
parasite_formatted <- paste0(round(parasite$ci_lower, 2), " - ", round(parasite$mean, 2), " - ", round(parasite$ci_upper, 2))
host_formatted <- paste0(round(host$ci_lower, 2), " - ", round(host$mean, 2), " - ", round(host$ci_upper, 2))

# Create a row for Parasite vs Host comparison
comparison_row <- data.frame(
  Model = "Parasite x Host",
  Parasite = parasite_formatted,
  Host = host_formatted,
  Relative_Difference = round(relative_diff, 3),
  p_value = signif(Pcd_PxH$p.value, 3),
  stringsAsFactors = FALSE
)

# Append Parasite vs Host comparison to the results table
Pcd_results <- rbind(Pcd_results, comparison_row)



# 2. Loop through species pairs
for (pair in species_pairs) {
  # Subset data for the current pair
  subset_data <- subset(PitMembrane_data, ssp %in% pair)
  
  # Perform t-test
  t <- t.test(pcd ~ ssp, data = subset_data, var.equal = FALSE, na.rm = TRUE)
  print(t)
  
  # Calculate group-level statistics
  res <- subset_data %>% 
    group_by(Group = ssp) %>% 
    summarise(
      mean = mean(pcd, na.rm = TRUE),
      ci_lower = mean - qt(0.975, df = n() - 1) * (sd(pcd, na.rm = TRUE) / sqrt(n())),
      ci_upper = mean + qt(0.975, df = n() - 1) * (sd(pcd, na.rm = TRUE) / sqrt(n()))
    )
  
  # Extract values for both groups
  group1 <- res[1, ]
  group2 <- res[2, ]
  
  # Calculate Relative Difference
  relative_diff <- (group2$mean - group1$mean) / group1$mean
  print(-(group2$mean - group1$mean))
  # Format the columns
  group1_formatted <- paste0(round(group1$ci_lower, 2), " - ", round(group1$mean, 2), " - ", round(group1$ci_upper, 2))
  group2_formatted <- paste0(round(group2$ci_lower, 2), " - ", round(group2$mean, 2), " - ", round(group2$ci_upper, 2))
  
  # Append species pair comparison to the results table
  Pcd_results <- rbind(
    Pcd_results,
    data.frame(
      Model = paste(pair, collapse = " x "),
      Parasite = group1_formatted,
      Host = group2_formatted,
      Relative_Difference = round(relative_diff, 3),
      p_value = signif(t$p.value, 3)
    )
  )
}

# Print the results table
print(Pcd_results)



Tpm_results <- data.frame(
  Model = character(),
  Parasite = character(),
  Host = character(),
  Relative_Difference = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

# 1. Add the Parasite vs Host comparison
Tpm_PxH <- t.test(Tpm ~ parasitism, data = PitMembrane_data, var.equal = FALSE, na.rm = TRUE)
diff(Tpm_PxH$estimate)
# Calculate group-level statistics for Parasite and Host
group_stats <- PitMembrane_data %>% 
  group_by(Group = parasitism) %>% 
  summarise(
    mean = mean(Tpm, na.rm = TRUE),
    ci_lower = mean - qt(0.975, df = n() - 1) * (sd(Tpm, na.rm = TRUE) / sqrt(n())),
    ci_upper = mean + qt(0.975, df = n() - 1) * (sd(Tpm, na.rm = TRUE) / sqrt(n()))
  )

# Extract values for Parasite and Host
parasite <- group_stats[group_stats$Group == "Parasite", ]
host <- group_stats[group_stats$Group == "Host", ]

# Calculate Relative Difference
relative_diff <- (host$mean - parasite$mean) / parasite$mean

# Format the Parasite and Host columns
parasite_formatted <- paste0(round(parasite$ci_lower, 2), " - ", round(parasite$mean, 2), " - ", round(parasite$ci_upper, 2))
host_formatted <- paste0(round(host$ci_lower, 2), " - ", round(host$mean, 2), " - ", round(host$ci_upper, 2))

# Create a row for Parasite vs Host comparison
comparison_row <- data.frame(
  Model = "Parasite x Host",
  Parasite = parasite_formatted,
  Host = host_formatted,
  Relative_Difference = round(relative_diff, 3),
  p_value = signif(Tpm_PxH$p.value, 3),
  stringsAsFactors = FALSE
)

# Append Parasite vs Host comparison to the results table
Tpm_results <- rbind(Tpm_results, comparison_row)



# 2. Loop through species pairs
for (pair in species_pairs) {
  # Subset data for the current pair
  subset_data <- subset(PitMembrane_data, ssp %in% pair)
  
  # Perform t-test
  t <- t.test(Tpm ~ ssp, data = subset_data, var.equal = FALSE, na.rm = TRUE)
  print(t)
  # Calculate group-level statistics
  res <- subset_data %>% 
    group_by(Group = ssp) %>% 
    summarise(
      mean = mean(Tpm, na.rm = TRUE),
      ci_lower = mean - qt(0.975, df = n() - 1) * (sd(Tpm, na.rm = TRUE) / sqrt(n())),
      ci_upper = mean + qt(0.975, df = n() - 1) * (sd(Tpm, na.rm = TRUE) / sqrt(n()))
    )
  
  # Extract values for both groups
  group1 <- res[1, ]
  group2 <- res[2, ]
  
  # Calculate Relative Difference
  relative_diff <- (group2$mean - group1$mean) / group1$mean
  print(-(group2$mean - group1$mean))
  # Format the columns
  group1_formatted <- paste0(round(group1$ci_lower, 2), " - ", round(group1$mean, 2), " - ", round(group1$ci_upper, 2))
  group2_formatted <- paste0(round(group2$ci_lower, 2), " - ", round(group2$mean, 2), " - ", round(group2$ci_upper, 2))
  
  # Append species pair comparison to the results table
  Tpm_results <- rbind(
    Tpm_results,
    data.frame(
      Model = paste(pair, collapse = " x "),
      Parasite = group1_formatted,
      Host = group2_formatted,
      Relative_Difference = round(relative_diff, 3),
      p_value = signif(t$p.value, 3)
    )
  )
}

# Print the results table
print(Tpm_results)
