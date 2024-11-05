

######################################################################
#
# Victor Sibinelli (victor.sibinelli@usp.br / sibinelli95@gmail.com)
# 13/07/2024
# Script 03.2 - Graphics - Vessels
#################################################################
library(here)
rm(list=ls())
source(here("scripts","Functions.R"))
HydraulicData <- read.csv(here("data", "processed", "HydraulicData.csv"))
# List of file names to load into a list
files <- list.files(path = here("outputs", "tables"), pattern = "Vessels_MonteCarlo_CI95_.*\\.csv", full.names = TRUE)
# Load each file as a data frame in a list
CI95_MT <- lapply(files, fread)
names(CI95_MT) <- gsub(".*CI95_(.*)\\.csv", "\\1", basename(files))  # Assign names based on file names
Vessels_Pvalues <-read.csv(here("outputs","tables","Vessels_MonteCarlo_Pvalues.csv"))
Vessels_AIC <- list(
  HydraulicDiameter = fread(here("outputs", "tables", "HydraulicDiameter_results.R")),
  VesselDensity = fread(here("outputs", "tables", "vdensity_results.R")),
  Kmax = fread(here("outputs", "tables", "Kmax_results.R"))
) %>% lapply(as_tibble)

# Set the desired order for the groups
desired_order <- c("Parasite", "Host",
                   "Psittacanthus robustus", "Vochysia thyrsoidea",
                   "Phoradendron perrotettii", "Tapirira guianensis",
                   "Struthanthus rhynchophyllus", "Tipuana tipu",
                   "Viscum album", "Populus nigra")


# Define short names
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
relevel_factors(ls())


##############################Ressampling
###Hydraulic Diameter
Vessels_AIC[[1]]

g <- HydraulicData %>%  
  ggplot(aes_string(x = "ssp", y = " HydraulicDiameter", fill = "parasitism")) +
  geom_jitter(aes(color = pic),
              size = 1, alpha = 0.4, position = position_jitter(width = 0.3)) +
  geom_boxplot(color = "black", alpha = 1, position = position_dodge(width = 0.8)) +
  scale_fill_manual(
    values = c("Parasite" = "firebrick", "Host" = "grey"),
    name = "Parasitism"
  ) +
  scale_color_viridis_d(option = "D") +
  geom_vline(xintercept = c(2.5, 4.5, 6.5)) +
  theme_classic() +
  scale_x_discrete(labels = short_names) +
  labs(title = "Hydraulic Diameter",
       x = "Species",
       y = "Hydraulic Diameter (µm)") +  # Use dot for multiplication
  annotate("text", x = seq_along(unique(HydraulicData$ssp)),
           y = max(HydraulicData$HydraulicDiameter)*1.1, label = c("A","B","A","A","A","B"), size = 6) +
  annotate("text", x = c(1.5,3.5,5.5,7.5),
           y = max(HydraulicData$HydraulicDiameter)*1.1, 
           label = "*", 
           size = 6) +
  theme(legend.position = "right",
        axis.text.x = element_text(size = 12),        # X-axis tick labels size
        axis.text.y = element_text(size = 12)) +      # Y-axis tick labels size
  guides(color = "none")  # Remove the legend for `ssp`

print(g)



###Vessel density

Vessels_AIC[2]

g <- HydraulicData %>%  
  ggplot(aes_string(x = "ssp", y = "vdensity", fill = "parasitism")) +
  geom_jitter(aes(color = pic),
              size = 1, alpha = 0.4, position = position_jitter(width = 0.3)) +
  geom_boxplot(color = "black", alpha = 1, position = position_dodge(width = 0.8)) +
  scale_fill_manual(
    values = c("Parasite" = "firebrick", "Host" = "grey"),
    name = "Parasitism"
  ) +
  scale_color_viridis_d(option = "D") +
  geom_vline(xintercept = c(2.5, 4.5, 6.5)) +
  theme_classic() +
  scale_x_discrete(labels = short_names) +
  labs(title = "Vessel Density",
       x = "Species",
       y = "Vessel Density (Vessels/mm²)") +  # Use dot for multiplication
  annotate("text", x = seq_along(unique(HydraulicData$ssp)),
           y = max(HydraulicData$vdensity)*1.1, label = c("A","A","A","A","A","B"), size = 6) +
  annotate("text", x = c(1.5,3.5,5.5,7.5),
           y = max(HydraulicData$vdensity)*1.1, 
           label = "*", 
           size = 6) +
  theme(legend.position = "right",
        axis.text.x = element_text(size = 12),        # X-axis tick labels size
        axis.text.y = element_text(size = 12)) +      # Y-axis tick labels size
  guides(color = "none")  # Remove the legend for `ssp`

print(g)


#Kmax

Vessels_AIC[3]
Vessels_Pvalues
  g <- HydraulicData %>%  
    ggplot(aes_string(x = "ssp", y = "Kmax", fill = "parasitism")) +
    geom_jitter(aes(color = pic),
                size = 1, alpha = 0.4, position = position_jitter(width = 0.3)) +
    geom_boxplot(color = "black", alpha = 1, position = position_dodge(width = 0.8)) +
    scale_fill_manual(
      values = c("Parasite" = "firebrick", "Host" = "grey"),
      name = "Parasitism"
    ) +
    scale_color_viridis_d(option = "D") +
    geom_vline(xintercept = c(2.5, 4.5, 6.5)) +
    theme_classic() +
    scale_x_discrete(labels = short_names) +
    labs(title = "Kmax",
         x = "Species",
         y = "Max Conductivity (kg·s·MPa⁻¹·m⁻²)") +  # Use dot for multiplication
    annotate("text", x = seq_along(unique(HydraulicData$ssp)),
             y = max(HydraulicData$Kmax)*1.1, label = c("A","B","A","B","A","B"), size = 6) +
    annotate("text", x = c(1.5,3.5,5.5,7.5),
             y = max(HydraulicData$Kmax)*1.1, 
             label = "*", 
             size = 6) +
    theme(legend.position = "right",
          axis.text.x = element_text(size = 12),        # X-axis tick labels size
          axis.text.y = element_text(size = 8)) +      # Y-axis tick labels size
    guides(color = "none")  # Remove the legend for `ssp`
  
  print(g)

### graph of ssp ci95
  
  for (i in seq_along(CI95_MT)) {
    # Convert the data.table to a data.frame for compatibility with ggplot
    current_data <- as.data.frame(CI95_MT[[i]]) 
    
    # Set the factor levels for 'ssp' based on desired_order
    current_data$ssp <- factor(current_data$ssp, levels = rev(desired_order))
    
    # Now, get unique levels from the updated factor
    unique_levels <- length(unique(current_data$ssp))  
    
    g <- ggplot(current_data, aes(ssp, Mean)) +
      geom_point(size = 4, aes(color = ssp)) +
      geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2) +
      coord_flip() +
      labs(title = paste("95% Confidence Intervals for", names(CI95_MT)[i]),
           x = "Effect", y = "Estimate") +
      scale_color_manual(
        values = rep(c("black", "firebrick"), length.out = unique_levels)  # Adjust colors for unique levels
      ) +
      theme_minimal() +
      theme(legend.position = "none")  # Remove legend if not needed
    
    print(g)
  }
  
 ###############################################################################
  ##LME
  g <- HydraulicData %>%  
    ggplot(aes_string(x = "ssp", y = " HydraulicDiameter", fill = "parasitism")) +
    geom_jitter(aes(color = pic),
                size = 1, alpha = 0.4, position = position_jitter(width = 0.3)) +
    geom_boxplot(color = "black", alpha = 1, position = position_dodge(width = 0.8)) +
    scale_fill_manual(
      values = c("Parasite" = "firebrick", "Host" = "grey"),
      name = "Parasitism"
    ) +
    scale_color_viridis_d(option = "D") +
    geom_vline(xintercept = c(2.5, 4.5, 6.5)) +
    theme_classic() +
    scale_x_discrete(labels = short_names) +
    labs(title = "Hydraulic Diameter",
         x = "Species",
         y = "Hydraulic Diameter (µm)") +  # Use dot for multiplication
    annotate("text", x = seq_along(unique(HydraulicData$ssp)),
             y = max(HydraulicData$HydraulicDiameter)*1.1, label = c("A","B","A","B","A","B"), size = 6) +
    theme(legend.position = "right",
          axis.text.x = element_text(size = 12),        # X-axis tick labels size
          axis.text.y = element_text(size = 12)) +      # Y-axis tick labels size
    guides(color = "none")  # Remove the legend for `ssp`
  
  print(g)
  
  
  
  ###Vessel density
  
  g <- HydraulicData %>%  
    ggplot(aes_string(x = "ssp", y = "vdensity", fill = "parasitism")) +
    geom_jitter(aes(color = pic),
                size = 1, alpha = 0.4, position = position_jitter(width = 0.3)) +
    geom_boxplot(color = "black", alpha = 1, position = position_dodge(width = 0.8)) +
    scale_fill_manual(
      values = c("Parasite" = "firebrick", "Host" = "grey"),
      name = "Parasitism"
    ) +
    scale_color_viridis_d(option = "D") +
    geom_vline(xintercept = c(2.5, 4.5, 6.5)) +
    theme_classic() +
    scale_x_discrete(labels = short_names) +
    labs(title = "Vessel Density",
         x = "Species",
         y = "Vessel Density (Vessels/mm²)") +  # Use dot for multiplication
    annotate("text", x = seq_along(unique(HydraulicData$ssp)),
             y = max(HydraulicData$vdensity)*1.1, label = c("A","B","A","B","A","B"), size = 6) +
    theme(legend.position = "right",
          axis.text.x = element_text(size = 12),        # X-axis tick labels size
          axis.text.y = element_text(size = 12)) +      # Y-axis tick labels size
    guides(color = "none")  # Remove the legend for `ssp`
  
  print(g)
  
  
  #Kmax
  g <- HydraulicData %>%  
    ggplot(aes_string(x = "ssp", y = "Kmax", fill = "parasitism")) +
    geom_jitter(aes(color = pic),
                size = 1, alpha = 0.4, position = position_jitter(width = 0.3)) +
    geom_boxplot(color = "black", alpha = 1, position = position_dodge(width = 0.8)) +
    scale_fill_manual(
      values = c("Parasite" = "firebrick", "Host" = "grey"),
      name = "Parasitism"
    ) +
    scale_color_viridis_d(option = "D") +
    geom_vline(xintercept = c(2.5, 4.5, 6.5)) +
    theme_classic() +
    scale_x_discrete(labels = short_names) +
    labs(title = "Kmax",
         x = "Species",
         y = "Max Conductivity (kg·s·MPa⁻¹·m⁻²)") +  # Use dot for multiplication
    annotate("text", x = seq_along(unique(HydraulicData$ssp)),
             y = max(HydraulicData$Kmax)*1.1, label = c("A","B","A","B","A","B"), size = 6) +
    theme(legend.position = "right",
          axis.text.x = element_text(size = 12),        # X-axis tick labels size
          axis.text.y = element_text(size = 8)) +      # Y-axis tick labels size
    guides(color = "none")  # Remove the legend for `ssp`
  
  print(g)
  
  