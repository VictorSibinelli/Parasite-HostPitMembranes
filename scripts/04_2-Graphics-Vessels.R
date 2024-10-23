

######################################################################
#
# Victor Sibinelli (victor.sibinelli@usp.br / sibinelli95@gmail.com)
# 13/07/2024
# Script 03.2 - Graphics - Vessels
#################################################################
library(here)
source(file=here("scripts","03_2-DataAnalysis-Vessels.R"))



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

## Create a vector of column names
column_names <- colnames(boot_diffs)
# Loop through each column in the boot_diffs data frame using lapply
sapply(seq_along(column_names), function(i) {
  x <- boot_diffs[, i]  # Get the current column data
  col_name <- column_names[i]  # Get the current column name
  
  # Create a density plot
  plot(density(x, na.rm = TRUE), main = "", xlab = "Differences", ylab = "Density")
  
  # Add quantile lines
  abline(v = quantile(x, c(0.025, 0.975), na.rm = TRUE), lwd = 2, col = "black")
  
  # Add a vertical line for the observed value
  abline(v = x[1], col = "red", lwd = 2)
  
  # Calculate p-value and display it on the plot
  p_value <- sum(abs(x) >= abs(x[1]), na.rm = TRUE) / length(x)
  text(x = mean(x, na.rm = TRUE), y = max(density(x)$y, na.rm = TRUE) * 0.9, paste("p-value =",p_value))
  
  # Title for the plot using the current column name
  title(main = col_name)
})

##graph of ci95



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
          axis.text.y = element_text(size = 12)) +      # Y-axis tick labels size
    guides(color = "none")  # Remove the legend for `ssp`
  
  print(g)






### graph of ssp ci95