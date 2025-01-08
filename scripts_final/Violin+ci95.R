# Load necessary libraries
library(dplyr)
library(ggplot2)
library(here)

# Load data
Wall_data <- read.csv(here("data", "processed", "Wall_data.csv"))
VesselDiameter_data <- read.csv(here("data", "processed", "VesselDiameter_data.csv")) %>%
  group_by(ssp, indiv) %>%
  filter(VesselDiameter >= quantile(VesselDiameter, 0.9)) %>%
  ungroup()
Hydraulic_data <- read.csv(here("data", "processed", "HydraulicData.csv"))
PitFraction_data <- read.csv(here("data", "processed", "PitFraction_data.csv"))
PitDiOp_data <- read.csv(here("data", "processed", "PitDiOp_data.csv"))
PitMembrane_data <- read.csv(here("data", "processed", "PitMembrane_data.csv"))

# Load bootstrap confidence interval data
directory <- here("data", "processed", "ressampled")
file_list <- list.files(path = directory, pattern = "Medians_.*CI95.*\\.csv$", full.names = TRUE)

CI95_boot <- lapply(file_list, function(file) {
  read.csv(file) %>% select(-1)
})

names(CI95_boot) <- gsub("Medians_|_CI95\\.csv", "", basename(file_list))

# Load functions
source(here("scripts", "Functions.R"))

# Define a function to format confidence intervals
format_ci <- function(t_test_results) {
  data.frame(
    Group = names(t_test_results),
    CI = sapply(t_test_results, function(x) {
      paste0(
        round(x$conf.int[1], 2), " - ",
        round(x$estimate, 2), " - ",
        round(x$conf.int[2], 2)
      )
    }),
    stringsAsFactors = FALSE
  )
}

# List of species pairs for comparison
species_pairs <- list(
  c("Psittacanthus robustus", "Vochysia thyrsoidea"),
  c("Phoradendron perrotettii", "Tapirira guianensis"),
  c("Struthanthus rhynchophyllus", "Tipuana tipu"),
  c("Viscum album", "Populus nigra")
)

# Define short names for species
short_names <- c(
  "Psittacanthus robustus" = expression(italic("P. robustus")),
  "Vochysia thyrsoidea" = expression(italic("V. thyrsoidea")),
  "Phoradendron perrotettii" = expression(italic("P. perrotettii")),
  "Tapirira guianensis" = expression(italic("T. guianensis")),
  "Struthanthus rhynchophyllus" = expression(italic("S. rhynchophyllus")),
  "Tipuana tipu" = expression(italic("T. tipu")),
  "Viscum album" = expression(italic("V. album")),
  "Populus nigra" = expression(italic("P. nigra"))
)

# Compute confidence intervals for traits
CI95 <- lapply(seq_len(ncol(CI95_boot[[1]])), function(v) {
  ci_data <- data.frame(Grouping = character(0), Mean = numeric(0), Lower_CI = numeric(0), Upper_CI = numeric(0))
  
  for (i in seq_along(CI95_boot)) {
    data <- CI95_boot[[i]]
    Grouping <- names(CI95_boot)[i]
    mean_v <- mean(data[, v], na.rm = TRUE)
    ci_v <- quantile(data[, v], c(0.025, 0.975), na.rm = TRUE)
    
    ci_data <- rbind(ci_data, data.frame(Grouping = Grouping, Mean = mean_v, Lower_CI = ci_v[1], Upper_CI = ci_v[2]))
  }
  ci_data
})

names(CI95) <- colnames(CI95_boot[[1]])

# Generate and save plots
for (e in seq_along(CI95)) {
  data <- CI95[[e]]
  data$Grouping <- factor(data$Grouping, levels = c(
    "Parasite", "Host",
    "Psittacanthus robustus", "Vochysia thyrsoidea",
    "Phoradendron perrotettii", "Tapirira guianensis",
    "Struthanthus rhynchophyllus", "Tipuana tipu",
    "Viscum album", "Populus nigra"
  ))
  
  g <- data %>%
    ggplot(aes(Grouping, Mean)) +
    geom_point(size = 4, aes(color = Grouping)) +
    geom_errorbar(aes(ymin = Lower_CI, ymax = Upper_CI), width = 0.2) +
    coord_flip() +
    labs(
      title = paste0("Bootstrap WT 95% Confidence Intervals ", names(CI95)[e]),
      x = "Effect", 
      y = "Estimate"
    ) +
    scale_color_manual(values = rep(c("firebrick", "black"), length.out = nlevels(data$Grouping))) +
    theme_minimal() +
    theme(legend.position = "none")
  
  print(g)
}

# Format and combine confidence interval results
pcd_CI95 <- tapply(PitMembrane_data$pcd, PitMembrane_data$parasitism, t.test)
Tpm_CI95 <- tapply(PitMembrane_data$Tpm, PitMembrane_data$parasitism, t.test)
pcd_CI95_ssp <- tapply(PitMembrane_data$pcd, PitMembrane_data$ssp, t.test)
Tpm_CI95_ssp <- tapply(PitMembrane_data$Tpm, PitMembrane_data$ssp, t.test)

Pcd_CI_df <- format_ci(pcd_CI95)
Tpm_CI_df <- format_ci(Tpm_CI95)
Pcd_CI_ssp_df <- format_ci(pcd_CI95_ssp)
Tpm_CI_ssp_df <- format_ci(Tpm_CI95_ssp)

Pcd_CI_combined <- rbind(Pcd_CI_df, Pcd_CI_ssp_df)
Tpm_CI_combined <- rbind(Tpm_CI_df, Tpm_CI_ssp_df)

formatted_columns <- lapply(CI95, function(data) {
  paste0(round(data$Lower_CI, 2), " - ", round(data$Mean, 2), " - ", round(data$Upper_CI, 2))
})

result_df <- do.call(cbind, formatted_columns)
rownames(result_df) <- CI95[[1]]$Grouping
result_df <- as.data.frame(result_df)
result_df$Pcd <- Pcd_CI_combined$CI[match(rownames(result_df), Pcd_CI_combined$Group)]
result_df$Tpm <- Tpm_CI_combined$CI[match(rownames(result_df), Tpm_CI_combined$Group)]

print(result_df)


Wall_data %>% ggplot()


