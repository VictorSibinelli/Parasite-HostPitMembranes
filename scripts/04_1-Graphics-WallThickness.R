######################################################################
#
# Victor Sibinelli (victor.sibinelli@usp.br / sibinelli95@gmail.com)
# 13/07/2024
# Scrpt 04.3 -Graphics- Pits
#################################################################
library(here)
source(here("scripts", "03_1-DataAnalysis-WallThickness.R"))


x_pos <- c(1, 3, 5, 7, 9, 11, 13, 15)
par(mar = c(7, 6, 1, 2) + 0.1, mgp = c(3, 1, 0), mfrow = c(1, 1))


boxplot(
  data = wdata_clean, wthickness~ ssp, na.rm = T, las = 2, cex = 1,
  ylab = "", xlab = NA, tcl = T, xaxt = "n",
  col = c("firebrick", "grey"), at = c(1, 3, 5, 7, 9, 11, 13, 15),
  ylim = c(0, 8), cex.lab = 1.5, cex.axis = 1.5, whisklwd = 2, staplelwd = 2
)
title(ylab = "Intervessel Wall Thickness (µm)", line = 4, cex.lab = 1.5)
text(
  x = x_pos, y =-2.5, c("P. robustus", "V. thyrsoidea", "P. perrottetti", "T. guianensis", "S. rhynchophyllus", "T. tipu", "V. album", "P. nigra"),
  xpd = NA, cex = 1.5, srt = 35, col = "black", adj = 0.50, font = 3
)
axis(1, at = x_pos, labels = NA)
text(x = x_pos+0.5, y = 6,2, labels = c("A", "B", "A", "A", "A", "B"), cex = 1.5)
abline(v = c(4, 8, 12), lty = 2)

text(x = c(2, 6, 10, 14,2, 6, 10, 14),
     y = c(VWall_results$ParasiteMean,VWall_results$HostMean), "-", cex = 2, col = "black")
segments(x0 = c(2, 6, 10, 14), y0 = VWall_results$ParasiteMean, y1 = VWall_results$HostMean, lwd = 2)
text(x = c(2, 6, 10, 14) + 0.3, y = (VWall_results$ParasiteMean+VWall_results$HostMean)/2, round(VWall_results$MeanDifference, 3), cex = 1.2)




# Update `wdata` to include the `parasitism` variable
wdata <- wdata %>%
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

# Create the plot

g <- ggplot(wdata, aes(x = ssp, y = wthickness, fill = parasitism)) +
  geom_boxplot(color = "black", alpha = 1, position = position_dodge(width = 0.8)) +
  geom_jitter(aes(color = label), size = 1, alpha = 0.4, position = position_jitter(width = 0.3)) +
  scale_fill_manual(
    values = c("Parasite" = "firebrick", "Host" = "grey"),
    name = "Parasitism"
  ) +
 geom_vline(xintercept = c(2.5,4.5,6.5)) +
  theme_classic() +
  scale_x_discrete(labels = short_names) +
  labs(title = "Vessel Wall Thickness",
       x = "Species",
       y = "Vessel Wall Thickness (µm)") +
  annotate("text", x = seq_along(unique(wdata$ssp)),
           y = max(wdata$wthickness) + 5, label = "A", size = 6)+
  theme(legend.position = "right") +  # Legend on the right
  guides(color = "none")  # Remove the legend for `ssp`
g
ggsave(here("outputs","figs","vessel_wall_thickness_plot.png"), plot = g, dpi = 600, width = 10, height = 7)

