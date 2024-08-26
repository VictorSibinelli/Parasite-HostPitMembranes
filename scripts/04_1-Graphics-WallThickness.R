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


