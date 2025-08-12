library(tidyverse)
library(here)
library(patchwork)


# Load bootstrap confidence interval data
directory <- here("data", "processed", "ressampled")
file_list <- list.files(path = directory, pattern = "Medians_.*CI95.*\\.csv$", full.names = TRUE)

CI95_boot <- lapply(file_list, function(file) {
  read.csv(file) %>% select(-1)
})

names(CI95_boot) <- gsub("Medians_|_CI95\\.csv", "", basename(file_list))

# Load functions
source(here("scripts", "Functions.R"))
relevel_factors(ls())
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


# Define the trait names and their abbreviations
trait_names <- c(
  "VesselDiameter" = "D",
  "TopVesselDiameter" = "Dtop",
  "HydraulicDiameter" = "Dh",
  "VesselDensity" = "VD",
  "VesselFraction" = "Fv",
  "Kmax" = "Kmax",
  "Wthickness" = "Tvw",
  "PitOpening" = "Dpa",
  "PitDiameter" = "Dpit",
  "PitFraction" = "Fp",
  "pcd" = "Hpit",
  "Tpm" = "Tpm"
)



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

# Order the traits alphabetically by their abbreviation
ordered_trait_names <- trait_names

plots <- list()  # Create an empty list to store plots

for (e in seq_along(ordered_trait_names)) {
  # Get the full name of the trait based on the abbreviation
  trait_full_name <- names(ordered_trait_names)[e]
  data <- CI95[[trait_full_name]]
  
  # Reverse the order of 'ssp' levels
  data$Grouping <- fct_rev(factor(data$Grouping, levels = c(
    "Parasite", "Host",
    "Psittacanthus robustus", "Vochysia thyrsoidea",
    "Phoradendron perrotettii", "Tapirira guianensis",
    "Struthanthus rhynchophyllus", "Tipuana tipu",
    "Viscum album", "Populus nigra"
  )))
  
  # Create the plot
  g <- data %>%
    ggplot(aes(Grouping, Mean)) +
    geom_point(size = 5, aes(color = Grouping)) +
    geom_errorbar(aes(ymin = Lower_CI, ymax = Upper_CI), width = 1) +
    coord_flip() +
    labs(
      title = paste0("95% CI ", ordered_trait_names[trait_full_name]),  # Use abbreviation for title
      y = "Estimate"
    ) +
    scale_color_manual(values = rep(c("black", "firebrick"), length.out = nlevels(data$Grouping))) +
    theme_minimal() +
    scale_x_discrete(labels = short_names) +
    theme(
      legend.position = "none",
      axis.text.y = element_text(face = "italic", size = 15),
      axis.text.x = element_text(size = 15),  # Italicize group names
      plot.title = element_text(hjust = 0.5, size = 15),
      axis.title.x = element_blank(),
      axis.title.y = element_blank()  # Center align the title
    )
  
  plots[[e]] <- g  # Store each plot in the list
}

# Combine the plots
combined_plot1 <- wrap_plots(plots, ncol = 2)  # Adjust `ncol` to set the layout

# Save the combined plot at 600 DPI
ggsave(
  filename = here("outputs", "figs", "CI95.png"),  # Specify the file name
  plot = combined_plot1,  # The plot to save
  width = 10,  # Width in inches
  height = 14,  # Height in inches
  dpi = 600  # DPI for high-resolution
)
