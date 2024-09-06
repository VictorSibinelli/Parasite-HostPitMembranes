######################################################################
#
# Victor Sibinelli (victor.sibinelli@usp.br / sibinelli95@gmail.com)
# 13/07/2024
# Scrpt 04.3 -Graphics- Pits
#################################################################
library(here)
source(here("scripts", "03_3-DataAnalysis-Pits.R"))
rm(list=ls())

pitdata_clean <- read.csv(here("data", "processed", "pitdata_clean.csv"))
pitOdata <- read.csv(here("data", "processed", "pitOdata.csv"))
PitMembrane_results <- read.csv(here("outputs","tables","pit_membrane_diff.csv"))
pcd_results <- read.csv(here("outputs","tables","pcd_results.csv"))

# List of data frame names
dataframes <- ls()
# Relevel the factors for each data frame
for (df_name in dataframes) {
  df <- get(df_name) # Get the data frame by name
  
  if ("ssp" %in% colnames(df)) { # Check if 'ssp' column exists
    df$ssp <- factor(df$ssp, levels = c(
      "Psittacanthus robustus", "Vochysia thyrsoidea",
      "Phoradendron perrotettii", "Tapirira guianensis",
      "Struthanthus rhynchophyllus", "Tipuana tipu",
      "Viscum album", "Populus nigra"
    ))
    assign(df_name, df) # Assign the modified data frame back to its name
  }
  rm(df, df_name) # remove duplicated dataframe
}

pitOdata <- pitOdata %>%
  mutate(parasitism = case_when(
    ssp %in% c("Psittacanthus robustus", 
               "Phoradendron perrotettii", 
               "Struthanthus rhynchophyllus", 
               "Viscum album") ~ "Parasite",
    TRUE ~ "Host"
  ))

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
### Pit graphs

# # grouped parasite and host pit membrane
# png(here("outputs", "figs", "p-h_pitmembrane.png"), width = 9, height = 7, units = "in", res = 600)
# boxplot(
#   data = pitdata_clean, pitavg * 1000 ~ parasitism,
#   ylab = "Pit Membrane Thickness (nm)", xlab = "", tcl = T, xaxt = "n", at = c(1, 3),
#   col = c("grey", "firebrick"), ylim = c(0, 1000), cex.lab = 1.5, cex.axis = 1.5, whisklwd = 3, staplelwd = 3
# )
# legend("topleft", cex = 2, bty = "n", fill = c("grey", "firebrick"), legend = c("Host", "Parasite"))
# text(x = c(2), y = c(0.3120534, 0.5938832) * 1000, "--", cex = 2, col = "black")
# segments(x0 = 2, y0 = 0.3120534 * 1000, x1 = 2, y1 = 0.5938832 * 1000, lwd = 3)
# text(x = 2.2, y = (0.2818298 / 2 + 0.3120534) * 1000, "0.28", cex = 2)
# text(x = c(1, 3), y = c(2, 1.2) - 0.1, c("A", "B"), cex = 1.5)
# dev.off()

# grouped parasite host pcd

# png(here("outputs", "figs", "p-h_pitchamber.png"), width = 9, height = 7, units = "in", res = 600)
# par(mar = c(5, 6, 4, 2) + 0.1, mgp = c(3, 1, 0))
# boxplot(
#   data = pitdata_clean, pcd * 1000 ~ parasitism, na.rm = T, las = 2, ylim = c(0, 2000), cex.lab = 1.5, at = c(1, 3),
#   ylab = "", xlab = "", tcl = T, xaxt = "n", cex.axis = 1.5,
#   col = c("grey", "firebrick"), whisklwd = 3, staplelwd = 3, boxwex = 0.5
# )
# legend(cex = 2, "topright", bty = "n", fill = c("grey", "firebrick"), legend = c("Host", "Parasite"))
# segments(x0 = 2, y0 = 0.5342017 * 1000, x1 = 2, y1 = 0.8891008 * 1000, lwd = 3)
# text(
#   cex = 2, x = 2,
#   y = 1000 * c(mean(pitdata_clean$pcd[pitdata_clean$parasitism == "h"], na.rm = T), mean(pitdata_clean$pcd[pitdata_clean$parasitism == "p"], na.rm = T)), "--"
# )
# text(x = 2.2, y = (0.3608233 / 2 + 0.5342017) * 1000, 0.36, cex = 2)
# text(x = c(1, 3), y = c(2, 1.2), c("A", "B"), cex = 1.5)
# title(ylab = "Pit Chamber depth (nm)", line = 4, cex.lab = 1.5)
# dev.off()


# pairwise pit membrane thickness and pcd
# Save the plot as a PNG file with 600 dpi resolution
png(here("outputs", "figs", "pairwise_tpm.png"), 
    width = 24, height = 9, units = "in", res = 600)

# Set graphical parameters
par(mar = c(7, 6, 1, 2) + 0.1,bty="l", mgp = c(3, 1, 0), mfrow = c(1, 1))

# First boxplot
boxplot(
  data = pitdata_clean, 
  pitavg * 1000 ~ ssp, 
  na.rm = TRUE, 
  las = 2, 
  cex = 1,
  ylab = "", 
  xlab = NA, 
  tcl = TRUE, 
  xaxt = "n",
  col = c("firebrick", "grey"), 
  at = c(1, 4, 7, 10, 12, 15, 18, 21),
  ylim = c(0, 1100), 
  cex.lab = 1.5, 
  cex.axis = 1.5, 
  whisklwd = 4, 
  staplelwd = 4)

# Add y-axis label and species names
title(ylab = "Pit Membrane Thickness (nm)", line = 4, cex.lab = 1.5)
text(
  x = c(1, 4, 7, 10, 12, 15, 18, 21), 
  y = -150,
  labels = short_names,
  xpd = NA, 
  cex = 1.7, 
  srt = 35, 
  col = "black", 
  adj = 0.50, 
  font = 3
)

# Add x-axis and annotations
axis(1, at = c(1, 4, 7, 10, 12, 15, 18, 21), labels = NA)
text(x = c(1, 4, 7, 10, 12, 15, 18, 21), y = 1000, labels = c("A", "B"), cex = 2.5)
abline(v = c(5.5, 11, 16.5), lty = 2,lwd=2.5)
text(
  x = rep(c(2.5, 8.5, 13.5, 19.5), times = 2), 
  y = c(PitMembrane_results$ParasiteMean, PitMembrane_results$HostMean), 
  labels = "-", 
  cex = 2, 
  col = "black"
)
segments(
  x0 = c(2.5, 8.5, 13.5, 19.5), 
  y0 = PitMembrane_results$ParasiteMean, 
  y1 = PitMembrane_results$HostMean, 
  lwd = 2
)
text(
  x = c(2.5, 8.5, 13.5, 19.5) + 0.5, 
  y = (PitMembrane_results$ParasiteMean + PitMembrane_results$HostMean) / 2,
  labels = round(abs(PitMembrane_results$MeanDifference), 0), 
  cex = 2
)

dev.off()





# Second boxplot
png(here("outputs", "figs", "pairwise_pcd.png"), 
    width = 24, height = 9, units = "in", res = 600)

# Set graphical parameters
par(mar = c(7, 6, 1, 2) + 0.1,bty="l", mgp = c(3, 1, 0), mfrow = c(1, 1))


boxplot(
  data = pitdata_clean, 
  pcd * 1000 ~ ssp, 
  na.rm = TRUE, 
  las = 2, 
  cex = 1,
  ylab = "", 
  xlab = NA, 
  tcl = TRUE, 
  xaxt = "n",
  col = c("firebrick", "grey"), 
  at = c(1, 4, 7, 10, 12, 15, 18, 21),
  ylim = c(0, 2000), 
  cex.lab = 1.5, 
  cex.axis = 1.5, 
  whisklwd = 4, 
  staplelwd = 4
)

# Add y-axis label and species names
title(ylab = "Pit Chamber Depth (nm)", line = 4, cex.lab = 1.5)
text(
  x = c(1, 4, 7, 10, 12, 15, 18, 21), 
  y = -300, 
  labels = short_names,
  xpd = NA, 
  cex = 1.7, 
  srt = 35, 
  col = "black", 
  adj = 0.50, 
  font = 3
)

# Add x-axis and annotations
axis(1, at = c(1, 4, 7, 10, 12, 15, 18, 21), labels = NA)
text(x = c(1, 4, 7, 10, 12, 15, 18, 21), y = 2000, labels = c("A", "B"), cex = 2.5)
abline(v = c(5.5, 11, 16.5), lty = 2, lwd=2.5)
text(
  x = rep(c(2.5, 8.5, 13.5, 19.5), times = 2), 
  y = c(pcd_results$ParasiteMean, pcd_results$HostMean), 
  labels = "-", 
  cex = 2, 
  col = "black"
)
segments(
  x0 = c(2.5, 8.5, 13.5, 19.5), 
  y0 = pcd_results$ParasiteMean, 
  y1 = pcd_results$HostMean, 
  lwd = 2
)
text(
  x = c(2.5, 8.5, 13.5, 19.5) + 0.5, 
  y = (pcd_results$ParasiteMean + pcd_results$HostMean) / 2, 
  labels = round(pcd_results$MeanDifference), 
  cex = 2
)

# Close the graphics device
dev.off()


PD_plot <- ggplot(data=pitOdata, aes(y=PitDiameter, x=ssp, fill=parasitism))+
geom_jitter(aes(color = label), size = 1, alpha = 0.2, position = position_jitter(width = 0.3)) +
  scale_fill_manual(
    values = c("Parasite" = "firebrick", "Host" = "grey"),
    name = "Parasitism"
  ) +
  geom_boxplot(color = "black", alpha = 1, position = position_dodge(width = 0.8))+
  geom_vline(xintercept = c(2.5,4.5,6.5)) +
  theme_classic() +
  scale_x_discrete(labels = short_names) +
  scale_y_continuous(limits=c(0,13), breaks = seq(2,16,by=2))+
  labs(title = "Intervessel Pit Diameter",
       x = "Species",
       y = "Intervessel Pit Diameter (µm)") +
  annotate("text", x = seq_along(unique(pitOdata$ssp)),
           y = max(pitOdata$PitDiameter) + 1,label = rep(c("A", "B"), 
           length.out = length(seq_along(unique(pitOdata$ssp)))), size = 4)+
  theme(legend.position = "right") +  # Legend on the right
  guides(color = "none")  # Remove the legend for `ssp`

ggsave(here("outputs","figs" ,"PitDiameter.png"), plot = PD_plot, dpi = 600, width = 10, height = 7)



PO_plot <- ggplot(data=pitOdata, aes(y=PitOpening, x=ssp, fill=parasitism))+
  geom_jitter(aes(color = label), size = 1, alpha = 0.2, position = position_jitter(width = 0.3)) +
  scale_fill_manual(
    values = c("Parasite" = "firebrick", "Host" = "grey"),
    name = "Parasitism"
  ) +
  geom_boxplot(color = "black", alpha = 1, position = position_dodge(width = 0.8))+
  geom_vline(xintercept = c(2.5,4.5,6.5)) +
  scale_y_continuous(limits=c(0,11))+
  theme_classic() +
  scale_x_discrete(labels = short_names) +
  labs(title = "Intervessel Pit Opening",
       x = "Species",
       y = "Intervessel Pit Opening (µm)") +
  annotate("text", x = seq_along(unique(pitOdata$ssp)),
           y = 10.5, label = c("A","B","A","A","A","B","A"), size = 4)+
  theme(legend.position = "right") +  # Legend on the right
  guides(color = "none")  # Remove the legend for `ssp`
PO_plot
ggsave(here("outputs","figs" ,"PitOpening.png"), plot = PO_plot, dpi = 600, width = 10, height = 7)





# Wave graph base don Kaack 2021.
# Function definition

# Function definition
lucian_model <- function(n_pits, tpm_layer, prob_lp_layer = 0.25) {
  1 - (1 - prob_lp_layer^((tpm_layer + 20) / 40))^n_pits
}

# Create the background plot
number_of_pits <- seq(1e3, 30e3, 1e3)
pit_membrane_thickness <- seq(60, 1180, 40)
probability_of_large_pore_per_layers <- 0.25
background_plot <- outer(number_of_pits, pit_membrane_thickness, lucian_model)

# Define the range of parasitic plants
ptm_range_plant1 <- seq(338, 1021, 1) # Min to max range
plant1_plot <- outer(number_of_pits, ptm_range_plant1, lucian_model)
mean_ptm_plant1 <- mean(pitdata_clean$pitavg[pitdata_clean$parasitism == "p"], na.rm = T) * 1000

# Plot
wave <- plot_ly(
  x = ~pit_membrane_thickness,
  y = ~number_of_pits,
  z = ~background_plot
) %>%
  layout(
    scene = list(
      xaxis = list(title = "Tpm (nm)"),
      yaxis = list(title = "Number of pits in a <br> average vessel"),
      zaxis = list(title = "Probability of encountering <br> a large pore in a vessel")
    )
  ) %>%
  add_surface(colorbar = list(title = "Probability")) %>%
  add_surface(
    x = ~ptm_range_plant1, # add parasitic plant layer
    y = ~number_of_pits,
    z = ~ plant1_plot + 0.011,
    colorscale = list(
      list(0, "white"),
      list(0.05, "pink"),
      list(0.1, "red")
    ),
    colorbar = list(title = "Parasitic Plant Probability")
  ) %>%
  add_trace(
    x = rep(mean_ptm_plant1, length(number_of_pits)),
    y = number_of_pits,
    z = background_plot[, which.min(abs(pit_membrane_thickness - mean_ptm_plant1))] + 0.011,
    type = "scatter3d",
    mode = "lines",
    line = list(color = "black", width = 4),
    showlegend = FALSE
  )
wave
saveWidget(wave, here("outputs", "figs", "interactive_wave_plot.html"))
