######################################################################
#
# Victor Sibinelli (victor.sibinelli@usp.br / sibinelli95@gmail.com)
# 13/07/2024
# Scrpt 04.3 -Graphics- Pits
#################################################################
library(here)
source(here("scripts", "03_1-DataAnalysis-WallThickness.R"))

par(mar = c(7, 6, 1, 2) + 0.1, mgp = c(3, 1, 0), mfrow = c(1, 1))
boxplot(
  data = wdata_clean, wthickness~ ssp, na.rm = T, las = 2, cex = 1,
  ylab = "", xlab = NA, tcl = T, xaxt = "n",
  col = c("firebrick", "grey"), at = c(1, 3, 5, 7, 9, 11, 13, 15),
  ylim = c(0, 8), cex.lab = 1.5, cex.axis = 1.5, whisklwd = 2, staplelwd = 2
)
title(ylab = "Intervessel Wall Thickness (µm)", line = 4, cex.lab = 1.5)
text(
  x = c(1, 3, 5, 7, 9, 11, 13, 15), y =-2.5, c("P. robustus", "V. thyrsoidea", "P. perrottetti", "T. guianensis", "S. rhynchophyllus", "T. tipu", "V. album", "P. nigra"),
  xpd = NA, cex = 1.5, srt = 35, col = "black", adj = 0.50, font = 3
)
axis(1, at = c(1, 3, 5, 7, 9, 11, 13, 15), labels = NA)
text(x = c(1, 3, 5, 7, 9, 11, 13, 15), y = 6,2, labels = c("A", "B", "A", "A", "A", "B"), cex = 1.5)
abline(v = c(4, 8, 12), lty = 2)
text(x = c(2, 2, 6, 6, 10, 10, 14, 14),
     y = c(VWall_results$ParasiteMean,VWall_results$HostMean), "-", cex = 2, col = "black")
segments(x0 = c(2, 6, 10, 14), y0 = VWall_results$ParasiteMean, y1 = VWall_results$HostMean, lwd = 2)
text(x = c(2, 6, 10, 14) + 0.3, y = (VWall_results$ParasiteMean+VWall_results$HostMean)/2, round(VWall_results$MeanDifference, 3), cex = 1.2)


wdata <- wdata %>%
  mutate(parasitism = case_when(
    ssp %in% c("Psittacanthus robustus", 
               "Phoradendron perrotettii", 
               "Struthanthus rhynchophyllus", 
               "Viscum album") ~ "p",
    TRUE ~ "h"
  ))
# Print the updated dataframe to check
print(wdata)


#Define short names
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

ggplot(wdata, aes(x = x_pos, y = wthickness, fill = ssp)) +
  geom_boxplot(aes(color = ssp),color="black", alpha = 2, position = position_dodge(width = 0.8)) +
  geom_jitter(aes(color = label), size = 1, alpha = 0.25, position = position_jitter(width = 0.3))+
  scale_fill_manual(values = c("firebrick", "grey","firebrick", "grey","firebrick", "grey","firebrick", "grey")) +  # Colors for ssp fills
  theme_minimal()+
  scale_x_continuous(breaks = x_pos, labels = short_names) +
  # Labels and theme
  labs(title = "Vessel Wall Thickness",
       x = "Species",
       y = "Vessel Wall Thickness (µm)") +
theme(legend.position = "none",
      axis.title.x = element_text(size = 14),  # Increase x-axis title size
      axis.title.y = element_text(size = 14),  # Increase y-axis title size
      axis.text.x = element_text(size = 12),  # Increase x-axis tick label size
      axis.text.y = element_text(size = 12) 
      )
