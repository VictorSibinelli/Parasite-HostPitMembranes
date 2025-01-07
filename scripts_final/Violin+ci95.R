####trait violin plots


# Load data
Wall_data <- read.csv(here("data", "processed", "Wall_data.csv"))
VesselDiameter_data<- read.csv(here("data", "processed", "VesselDiameter_data.csv")) %>% group_by(ssp,indiv) %>%
  filter(VesselDiameter >= quantile(VesselDiameter, 0.9)) %>%
  ungroup()
Hydraulic_data<- read.csv(here("data", "processed", "HydraulicData.csv"))
PitFraction_data<- read.csv(here("data", "processed", "PitFraction_data.csv"))
PitDiOp_data<- read.csv(here("data", "processed", "PitDiOp_data.csv"))
PitMembrane_data<- read.csv(here("data", "processed", "PitMembrane_data.csv"))
directory <- here("data", "processed", "ressampled")  # Replace with your directory path
file_list <- list.files(path = directory, pattern = "Medians_.*CI95.*\\.csv$", full.names = TRUE)

# Create a named list where each element is the content of the corresponding .csv file
CI95_boot <- lapply(file_list, function(file) {
  read.csv(file) %>% select(-1)
})

names(CI95_boot) <- basename(file_list)
# Remove "Medians_" and "_CI95.csv" from the names
names(CI95_boot) <- gsub("Medians_|_CI95\\.csv", "", basename(file_list))


source(here("scripts", "Functions.R"))


# List of species pairs for comparison
species_pairs <- list(
  c("Psittacanthus robustus", "Vochysia thyrsoidea"),
  c("Phoradendron perrotettii", "Tapirira guianensis"),
  c("Struthanthus rhynchophyllus", "Tipuana tipu"),
  c("Viscum album", "Populus nigra")
)

relevel_factors(ls())

levels(Wall_data$ssp)
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

CI95 <- list()

# Iterate over each column index in the first bootstrapped data frame
for (v in seq_len(ncol(CI95_boot[[1]]))) {
  
  # Initialize a data frame to store results for the current column
  ci_data <- data.frame(
    Grouping = character(0),  # Species name or Groupingentifier
    Mean = numeric(0),  # Mean values for CI
    Lower_CI = numeric(0),  # Lower bound of CI
    Upper_CI = numeric(0)   # Upper bound of CI
  )
  
  # Iterate over each bootstrapped data set
  for (i in seq_along(CI95_boot)) {
    data <- CI95_boot[[i]]
    Grouping <- names(CI95_boot)[i]
    
    # Calculate the mean and 95% CI for the current column
    mean_v <- mean(data[, v], na.rm = TRUE)
    ci_v <- quantile(data[, v], c(0.025, 0.975), na.rm = TRUE)
    
    # Append the results to the ci_data data frame
    ci_data <- rbind(
      ci_data,
      data.frame(
        Grouping = Grouping,
        Mean = mean_v,
        Lower_CI = ci_v[1],
        Upper_CI = ci_v[2]
      )
    )
  }
  
  # Store the results for the current column in the CI95 list
  CI95[[v]] <- ci_data
}

names(CI95) <- colnames(CI95_boot[[1]])



for (e in seq_along(CI95)) {
  data <- CI95[[e]]
  
  # Ensure `Grouping` is a factor for consistent color mapping
  data$Grouping <- as.factor(data$Grouping)
  data$Grouping <- factor(data$Grouping, levels = c(
    "Parasite", "Host",
    "Psittacanthus robustus", "Vochysia thyrsoidea",
    "Phoradendron perrotettii", "Tapirira guianensis",
    "Struthanthus rhynchophyllus", "Tipuana tipu",
    "Viscum album", "Populus nigra"
  ))
  v <- names(CI95)[e]
  # Generate the plot
  g <- data %>%
    ggplot(aes(Grouping, Mean)) +
    geom_point(size = 4, aes(color = Grouping)) +
    geom_errorbar(aes(ymin = Lower_CI, ymax = Upper_CI), width = 0.2) +
    coord_flip() +
    labs(
      title = paste0("Bootstrap WT 95% Confidence Intervals ", v),
      x = "Effect", 
      y = "Estimate"
    ) +
    scale_color_manual(
      values = rep(c("firebrick", "black"), length.out = nlevels(data$Grouping))  # Alternate colors
    ) +
    theme_minimal() +
    theme(legend.position = "none")  # Remove legend if not needed
  
  # Print the plot
  print(g)
  
  # Save the plot to a file
 
}


