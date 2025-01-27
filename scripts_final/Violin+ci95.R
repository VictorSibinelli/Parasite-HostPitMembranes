# Load necessary libraries
library(dplyr)
library(ggplot2)
library(here)
library(patchwork)

# Load data
Wall_data <- read.csv(here("data", "processed", "Wall_data.csv"))
VesselDiameter_data <- read.csv(here("data", "processed", "VesselDiameter_data.csv"))
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
relevel_factors(ls())
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

plots <- list()  # Create an empty list to store plots

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
    geom_point(size = 5, aes(color = Grouping)) +
    geom_errorbar(aes(ymin = Lower_CI, ymax = Upper_CI), width = 1) +
    coord_flip() +
    labs(
      title = paste0("95% CI ", names(CI95)[e]), 
      y = "Estimate"
    ) +
    scale_color_manual(values = rep(c("firebrick", "black"), length.out = nlevels(data$Grouping))) +
    theme_minimal() +
    scale_x_discrete(labels = short_names)+
    theme(
      legend.position = "none",
      axis.text.y = element_text(face = "italic",size = 15),
      axis.text.x = element_text(size = 15),# Italicize group names
      plot.title = element_text(hjust = 0.5,size = 15),
      axis.title.x = element_blank(),
      axis.title.y = element_blank()# Center align the title
    )
  
  plots[[e]] <- g  # Store each plot in the list
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

Pcd_CI_combined <- Pcd_CI_combined %>%
  separate(CI, into = c("Lower_CI", "Mean", "Upper_CI"), sep = " - ", convert = TRUE)
Pcd_CI_combined$Group <- factor(Pcd_CI_combined$Group, levels = c(
  "Parasite", "Host",
  "Psittacanthus robustus", "Vochysia thyrsoidea",
  "Phoradendron perrotettii", "Tapirira guianensis",
  "Struthanthus rhynchophyllus", "Tipuana tipu",
  "Viscum album", "Populus nigra"
))
Tpm_CI_combined <- Tpm_CI_combined %>%
  separate(CI, into = c("Lower_CI", "Mean", "Upper_CI"), sep = " - ", convert = TRUE)
Tpm_CI_combined$Group <- factor(Tpm_CI_combined$Group, levels = c(
  "Parasite", "Host",
  "Psittacanthus robustus", "Vochysia thyrsoidea",
  "Phoradendron perrotettii", "Tapirira guianensis",
  "Struthanthus rhynchophyllus", "Tipuana tipu",
  "Viscum album", "Populus nigra"
))
# Plot for Pcd
Pcd_plot <- Pcd_CI_combined %>%
  ggplot(aes(x = Group, y = Mean)) +
  geom_point(size = 5, aes(color = Group)) +
  geom_errorbar(aes(ymin = Lower_CI, ymax = Upper_CI), width = 0.2) +
  coord_flip() +
  labs(
    title = "95% CI for Pcd", 
    y = "Pcd Estimate"
  ) +
  scale_color_manual(values = rep(c("firebrick", "black"),times=5)) +
  theme_minimal() +
  scale_x_discrete(labels = short_names)+
  theme(
    legend.position = "none",
    axis.text.y = element_text(face = "italic", size = 15),
    axis.text.x = element_text(size = 15),
    plot.title = element_text(hjust = 0.5, size = 15),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )

# Plot for Tpm
Tpm_plot <- Tpm_CI_combined %>%
  ggplot(aes(x = Group, y = Mean)) +
  geom_point(size = 5, aes(color = Group)) +
  geom_errorbar(aes(ymin = Lower_CI, ymax = Upper_CI), width = 0.2) +
  coord_flip() +
  labs(
    title = "95% CI for Tpm", 
    y = "Tpm Estimate"
  ) +
  scale_x_discrete(labels = short_names)+
  scale_color_manual(values = rep(c("firebrick", "black"),times=5)) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.y = element_text(face = "italic", size = 15),
    axis.text.x = element_text(size = 15),
    plot.title = element_text(hjust = 0.5, size = 15),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )

# Print the plots
print(Pcd_plot)
print(Tpm_plot)
plots[[11]] <- Pcd_plot
plots[[12]] <- Tpm_plot
# Combine all plots into one figure
combined_plot1 <- wrap_plots(plots, ncol = 2)  # Adjust `ncol` to set the layout
print(combined_plot1)
x11(width = 10,height = 14)






w_plot <- Wall_data %>% 
  ggplot(aes(x = ssp, y = WallThickness, fill = parasitism)) +
  geom_jitter(
              size = 1, alpha=0.5, position = position_jitter(width = 0.3)) +
  geom_violin(color = "black", alpha = 0.5, position = position_dodge(width = 0.9)) +
  scale_fill_manual(
    values = c("Parasite" = "firebrick", "Host" = "grey"),
    name = "Parasitism"
  ) +

  geom_vline(xintercept = c(2.5, 4.5, 6.5)) +
  theme_classic() +
  scale_x_discrete(labels = short_names, guide = guide_axis(angle = 45)) +
  labs(
       x = "Species",
       y = "Tvw (µm)") +
  annotate("text", x = c(1.5, 3.5, 5.5, 7.5),
           y = max(Wall_data$WallThickness, na.rm = TRUE)+1, label = c("ns"), size = 5) +
  theme(
    legend.position = "none",  # Removes the legend
    axis.title.x = element_blank(),  # Removes x-axis title
    axis.text.x = element_blank(),  # Removes x-axis tick labels
    axis.ticks.x = element_blank()  # Removes x-axis tick marks
  ) +
  guides(fill = "none", color = "none")+
  scale_y_continuous(limits = c(0,max(Wall_data$WallThickness)*1.1),
    expand = expansion(mult = c(0, 0.1))  # Adds extra space above the highest y value
  )  # Ensures no guides for `fill` or `color`
# Remove the legend for `ssp`
w_plot

d_plot <- VesselDiameter_data %>% ggplot( aes(x = ssp, y = VesselDiameter, fill = parasitism)) +
  geom_jitter(
              size = 0.5, alpha=0.3, position = position_jitter(width = 0.3)) +
  geom_violin(color = "black", alpha = 0.7, position = position_dodge(width = 0.9),adjust=2) +
  scale_fill_manual(
    values = c("Parasite" = "firebrick", "Host" = "grey"),
    name = "Parasitism"
  )+

  geom_vline(xintercept = c(2.5,4.5,6.5)) +
  theme_classic() +
  scale_x_discrete(labels = short_names,guide = guide_axis(angle = 45))  +
  labs(
       x = "Species",
       y = "D (µm)") +
  annotate("text", x = c(1.5,3.5,5.5,7.5),
           y = max(VesselDiameter_data$VesselDiameter) + 5, label = c("*"), size = 8)+
  annotate("text", x = c(1.5,3.5,5.5,7.5),
           y = max(VesselDiameter_data$VesselDiameter) + 5,
           label = "",
           size = 6)+
  theme(
    legend.position = "none",  # Removes the legend
    axis.title.x = element_blank(),  # Removes x-axis title
    axis.text.x = element_blank(),  # Removes x-axis tick labels
    axis.ticks.x = element_blank()  # Removes x-axis tick marks
  ) +
  guides(fill = "none", color = "none")+
  scale_y_continuous(limits = c(0,max(VesselDiameter_data$VesselDiameter)*1.15),
                     expand = expansion(mult = c(0, 0.1)))
                     
d_plot


dtop_plot <- VesselDiameter_data %>% group_by(indiv,ssp) %>% 
  filter(VesselDiameter>=quantile(VesselDiameter,0.9)) %>% ungroup() %>%
  ggplot( aes(x = ssp, y = VesselDiameter, fill = parasitism)) +
  geom_jitter(
    size = 0.5, alpha=0.3, position = position_jitter(width = 0.3)) +
  geom_violin(color = "black", alpha = 0.7, position = position_dodge(width = 0.9),adjust=2) +
  scale_fill_manual(
    values = c("Parasite" = "firebrick", "Host" = "grey"),
    name = "Parasitism"
  )+
  
  geom_vline(xintercept = c(2.5,4.5,6.5)) +
  theme_classic() +
  scale_x_discrete(labels = short_names,guide = guide_axis(angle = 45))  +
  labs(
    x = "Species",
    y = "Dtop (µm)") +
  annotate("text", x = c(1.5,3.5,5.5,7.5),
           y = max(VesselDiameter_data$VesselDiameter) + 5, label = c("*"), size = 8)+
  annotate("text", x = c(1.5,3.5,5.5,7.5),
           y = max(VesselDiameter_data$VesselDiameter) + 5,
           label = "",
           size = 6)+
  theme(
    legend.position = "none",  # Removes the legend
    axis.title.x = element_blank(),  # Removes x-axis title
    axis.text.x = element_blank(),  # Removes x-axis tick labels
    axis.ticks.x = element_blank()  # Removes x-axis tick marks
  ) +
  guides(fill = "none", color = "none")+
  scale_y_continuous(limits = c(0,max(VesselDiameter_data$VesselDiameter)*1.15),
                     expand = expansion(mult = c(0, 0.1)))
dtop_plot

hd_plot <- Hydraulic_data %>% 
  ggplot(aes(x = ssp, y = HydraulicDiameter, fill = parasitism)) +
  geom_jitter(
              size = 1, alpha=0.5, position = position_jitter(width = 0.3)) +
  geom_violin(color = "black", alpha = 0.6, position = position_dodge(width = 1),adjust=2) +
  scale_fill_manual(
    values = c("Parasite" = "firebrick", "Host" = "grey"),
    name = "Parasitism"
  ) +

  geom_vline(xintercept = c(2.5, 4.5, 6.5)) +
  theme_classic() +
  scale_x_discrete(labels = short_names, guide = guide_axis(angle = 45)) +
  scale_y_log10() +  # Transform y-axis to log scale
  labs(
       x = "Species",
       y = "Dh (µm)") +
  annotate("text", x = c(1.5, 3.5, 5.5, 7.5),
           y = max(Hydraulic_data$HydraulicDiameter)*0.9,  # Adjust for log scale
           label = c("*","ns","*","*"), size = c(8,5,8,8)) +
  annotate("text", x = c(1.5, 3.5, 5.5, 7.5),
           y = max(Hydraulic_data$HydraulicDiameter),
           label = "",
           size = 6) +
  theme(
    legend.position = "none",  # Removes the legend
    axis.title.x = element_blank(),  # Removes x-axis title
    axis.text.x = element_blank(),  # Removes x-axis tick labels
    axis.ticks.x = element_blank()  # Removes x-axis tick marks
  ) +
  guides(fill = "none", color = "none")
hd_plot

vd_plot <- Hydraulic_data %>% 
  ggplot(aes(x = ssp, y = VesselDensity, fill = parasitism)) +
  geom_jitter(
              size = 1, alpha = 0.7, position = position_jitter(width = 0.3)) +
  geom_violin(color = "black", alpha = 0.5, position = position_dodge(width = 0.9)) +
  scale_fill_manual(
    values = c("Parasite" = "firebrick", "Host" = "grey"),
    name = "Parasitism"
  ) +
  # scale_color_viridis_d(option = "C") +
  geom_vline(xintercept = c(2.5, 4.5, 6.5)) +
  theme_classic() +
  scale_x_discrete(labels = short_names, guide = guide_axis(angle = 45)) +
  scale_y_log10() +  # Transform y-axis to log scale
  labs(
       x = "Species",
       y = "VD (vessels·mm⁻¹)") +
  annotate("text", x = c(1.5, 3.5, 5.5, 7.5),
           y = max(Hydraulic_data$VesselDensity) * 1.1,  # Adjust for log scale
           label = c("ns"), size = c(5)) +
  annotate("text", x = c(1.5, 3.5, 5.5, 7.5),
           y = max(Hydraulic_data$VesselDensity) * 1.1,
           label = "",
           size = 6) +
  theme(
    legend.position = "none",  # Removes the legend
    axis.title.x = element_blank(),  # Removes x-axis title
    axis.text.x = element_blank(),  # Removes x-axis tick labels
    axis.ticks.x = element_blank()  # Removes x-axis tick marks
  ) +
  guides(fill = "none", color = "none")+
  scale_y_continuous(limits=c(0,310),
                       expand = expansion(mult = c(0, 0.1))  # Adds extra space above the highest y value
  )

vd_plot

fv_plot <-  Hydraulic_data %>% 
  ggplot(aes(x = ssp, y = VesselFraction, fill = parasitism)) +
  geom_jitter(
              size = 1, alpha=0.5, position = position_jitter(width = 0.3)) +
  geom_violin(color = "black", alpha = 0.5, position = position_dodge(width = 0.9)) +
  scale_fill_manual(
    values = c("Parasite" = "firebrick", "Host" = "grey"),
    name = "Parasitism"
  ) +

  geom_vline(xintercept = c(2.5, 4.5, 6.5)) +
  theme_classic() +
  scale_x_discrete(labels = short_names, guide = guide_axis(angle = 45)) +
  scale_y_log10() +  # Transform y-axis to log scale
  labs(
       x = "Species",
       y = "Fv (%)") +
  annotate("text", x = c(1.5, 3.5, 5.5, 7.5),
           y = max(Hydraulic_data$VesselFraction) * 1.1,  # Adjust for log scale
           label = c("*","*","ns","*"), size = c(8,8,5,8)) +
  annotate("text", x = c(1.5, 3.5, 5.5, 7.5),
           y = max(Hydraulic_data$VesselFraction) * 1.1,
           label = "",
           size = 6) +
  theme(
    legend.position = "none",  # Removes the legend
    axis.title.x = element_blank(),  # Removes x-axis title
    axis.text.x = element_blank(),  # Removes x-axis tick labels
    axis.ticks.x = element_blank()  # Removes x-axis tick marks
  ) +
  guides(fill = "none", color = "none")+
  scale_y_continuous(limits=c(0,35),
                     expand = expansion(mult = c(0, 0.1))  # Adds extra space above the highest y value
  )
fv_plot


kmax_plot <- Hydraulic_data %>% 
  ggplot(aes(x = ssp, y = Kmax, fill = parasitism)) +
  geom_jitter(
              size = 1, alpha=0.5, position = position_jitter(width = 0.3)) +
  geom_violin(color = "black", alpha = 0.5, position = position_dodge(width = 0.9)) +
  scale_fill_manual(
    values = c("Parasite" = "firebrick", "Host" = "grey"),
    name = "Parasitism"
  ) +

  geom_vline(xintercept = c(2.5, 4.5, 6.5)) +
  theme_classic() +
  scale_x_discrete(labels = short_names, guide = guide_axis(angle = 45)) +
  scale_y_log10(expand=c(0,0.2)) +  # Transform y-axis to log scale
  labs(
       x = "Species",
       y = "log Kmax (kg·s·MPa⁻¹·m⁻ µm)") +
  annotate("text", x = c(1.5, 3.5, 5.5, 7.5),
           y = max(Hydraulic_data$Kmax) * 1.1,  # Adjust for log scale
           label = c("*","*","ns","*"), size = c(8,8,5,8)) +
  annotate("text", x = c(1.5, 3.5, 5.5, 7.5),
           y = max(Hydraulic_data$Kmax) * 1.1,
           label = "",
           size = 6) +
  theme(
    legend.position = "none",  # Removes the legend
    axis.title.x = element_blank(),  # Removes x-axis title
    axis.text.x = element_blank(),  # Removes x-axis tick labels
    axis.ticks.x = element_blank()  # Removes x-axis tick marks
  ) +
  guides(fill = "none", color = "none")
  
kmax_plot


pd_plot <- PitDiOp_data %>% 
  ggplot(aes(x = ssp, y = PitDiameter, fill = parasitism)) +
  geom_jitter(
              size = 1, alpha=0.5, position = position_jitter(width = 0.3)) +
  geom_violin(color = "black", alpha = 0.5, position = position_dodge(width = 0.9)) +
  scale_fill_manual(
    values = c("Parasite" = "firebrick", "Host" = "grey"),
    name = "Parasitism"
  ) +

  geom_vline(xintercept = c(2.5, 4.5, 6.5)) +
  theme_classic() +
  scale_x_discrete(labels = short_names, guide = guide_axis(angle = 45)) +
  labs(
       x = "Species",
       y = "Dpit (µm)") +
  annotate("text", x = c(1.5, 3.5, 5.5, 7.5),
           y = max(PitDiOp_data$PitDiameter,na.rm = T) * 1.1,  # Adjust for log scale
           label = c("*"), size = c(8)) +
  annotate("text", x = c(1.5, 3.5, 5.5, 7.5),
           y = max(PitDiOp_data$PitDiameter) * 1.1,
           label = "",
           size = 6) +
  theme(
    legend.position = "none",  # Removes the legend
    axis.title.x = element_blank(),  # Removes x-axis title
    axis.text.x = element_blank(),  # Removes x-axis tick labels
    axis.ticks.x = element_blank()  # Removes x-axis tick marks
  ) +
  guides(fill = "none", color = "none")+
  scale_y_continuous(limits=c(0,22),
                     expand = expansion(mult = c(0, 0.1))  # Adds extra space above the highest y value
  )
pd_plot


po_plot <- PitDiOp_data %>% 
  ggplot(aes(x = ssp, y = PitOpening, fill = parasitism)) +
  geom_jitter(
              size = 1, alpha=0.5, position = position_jitter(width = 0.3)) +
  geom_violin(color = "black", alpha = 0.5, position = position_dodge(width = 0.9)) +
  scale_fill_manual(
    values = c("Parasite" = "firebrick", "Host" = "grey"),
    name = "Parasitism"
  ) +

  geom_vline(xintercept = c(2.5, 4.5, 6.5)) +
  theme_classic() +
  scale_x_discrete(labels = short_names, guide = guide_axis(angle = 45)) +
  labs(
       x = "Species",
       y = "Dop (µm)") +
  annotate("text", x = c(1.5, 3.5, 5.5, 7.5),
           y = max(PitDiOp_data$PitOpening,na.rm = T) * 1.1,  # Adjust for log scale
           label = c("ns","ns","*","*"), size = c(5,5,8,8)) +
  annotate("text", x = c(1.5, 3.5, 5.5, 7.5),
           y = max(PitDiOp_data$PitOpening) * 1.1,
           label = "",
           size = 6) +
  theme(
    legend.position = "none",  # Removes the legend
    axis.title.x = element_blank(), 
    axis.text.x = element_blank()# Removes x-axis title
 # Removes x-axis tick marks
  ) +
  guides(fill = "none", color = "none")+
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))  # Adds extra space above the highest y value
  )
  
po_plot



fp_plot <- PitFraction_data %>% 
  ggplot(aes(x = ssp, y = PitFraction, fill = parasitism)) +
  geom_jitter(
    size = 1, alpha=0.5, position = position_jitter(width = 0.3)) +
  geom_violin(color = "black", alpha = 0.5, position = position_dodge(width = 0.9)) +
  scale_fill_manual(
    values = c("Parasite" = "firebrick", "Host" = "grey"),
    name = "Parasitism"
  ) +
  
  geom_vline(xintercept = c(2.5, 4.5, 6.5)) +
  theme_classic() +
  scale_x_discrete(labels = short_names, guide = guide_axis(angle = 45)) +
  labs(
    x = "Species",
    y = "Fp (%)") +
  annotate("text", x = c(1.5, 3.5, 5.5, 7.5),
           y = max(PitFraction_data$PitFraction) * 1.1,  # Adjust for log scale
           label = c("ns","ns","*","*"), size = c(5,5,8,8)) +
  annotate("text", x = c(1.5, 3.5, 5.5, 7.5),
           y = max(PitFraction_data$PitFraction) * 1.1,
           label = "",
           size = 6) +
  theme(
    legend.position = "none"
  ) +
  guides(fill = "none", color = "none")+
  scale_y_continuous(limits=c(0,110),
                     expand = expansion(mult = c(0, 0.1))  # Adds extra space above the highest y value
  )
fp_plot

pcd_plot <- PitMembrane_data %>% 
  ggplot(aes(x = ssp, y = pcd, fill = parasitism)) +
  geom_jitter(
              size = 1, alpha=0.5, position = position_jitter(width = 0.3)) +
  geom_violin(color = "black", alpha = 0.5, position = position_dodge(width = 0.9)) +
  scale_fill_manual(
    values = c("Parasite" = "firebrick", "Host" = "grey"),
    name = "Parasitism"
  ) +

  geom_vline(xintercept = c(2.5, 4.5, 6.5)) +
  theme_classic() +
  scale_x_discrete(labels = short_names, guide = guide_axis(angle = 45)) +
  labs(
       x = "Species",
       y = "Hpit (nm)") +
  annotate("text", x = c(1.5, 3.5, 5.5, 7.5),
           y = max(PitMembrane_data$pcd,na.rm = T) * 1.1,  # Adjust for log scale
           label = c("*"), size = c(8)) +
  annotate("text", x = c(1.5, 3.5, 5.5, 7.5),
           y = max(PitMembrane_data$pcd) * 1.1,
           label = "",
           size = 6) +
  theme(legend.position = "right",
        axis.text.x = element_text(size = 10),  # X-axis tick labels size
        axis.text.y = element_text(size = 10)   # Y-axis tick labels size
  ) +
  guides(color = "none",fill="none")+
  scale_y_continuous(limits=(c(0,2100)),
    expand = expansion(mult = c(0, 0.1))  # Adds extra space above the highest y value
  )
pcd_plot


Tpm_plot <- PitMembrane_data %>% 
  ggplot(aes(x = ssp, y = Tpm, fill = parasitism)) +
  geom_jitter(
              size = 1, alpha=0.5, position = position_jitter(width = 0.3)) +
  geom_violin(color = "black", alpha = 0.5, position = position_dodge(width = 0.9)) +
  scale_fill_manual(
    values = c("Parasite" = "firebrick", "Host" = "grey"),
    name = "Parasitism"
  ) +

  geom_vline(xintercept = c(2.5, 4.5, 6.5)) +
  theme_classic() +
  scale_x_discrete(labels = short_names, guide = guide_axis(angle = 45)) +
  labs(,
       x = "Species",
       y = "Tpm (nm)") +
  annotate("text", x = c(1.5, 3.5, 5.5, 7.5),
           y = max(PitMembrane_data$Tpm,na.rm = T) * 1.1,  # Adjust for log scale
           label = c("*"), size = c(8)) +
  scale_y_continuous(limits = c(0,1200),
    expand = expansion(mult = c(0, 0.1))  # Adds extra space above the highest y value
  )+
  annotate("text", x = c(1.5, 3.5, 5.5, 7.5),
           y = max(PitMembrane_data$Tpm) * 1.1,
           label = "",
           size = 6) +
  theme(legend.position = "right",
        axis.text.x = element_text(size = 10),  # X-axis tick labels size
        axis.text.y = element_text(size = 10),
        plot.margin = margin(1,1,1,1)
  ) +
  guides(color = "none")

Tpm_plot


combined_plot2 <- 
  (d_plot + hd_plot + dtop_plot) /
  (fv_plot + kmax_plot + w_plot) /
  ( vd_plot+ pd_plot + po_plot) /
  (pcd_plot + Tpm_plot +fp_plot) +
  plot_layout(guides = "collect") &  # Collect legends and place them
  theme(legend.position = "bottom")            # Place the legend at the bottom

x11(width = 12,height = 12)
ggsave(file=here("outputs","figs","trait_violin_plot.png"),
       combined_plot2,units = "in",dpi = 600,height = 14,width = 10)
ggsave(file=here("outputs","figs","trait_CI95.png"),
       combined_plot1,units = "in",dpi = 600,height = 14,width = 10)
