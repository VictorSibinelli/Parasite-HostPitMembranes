####covariance
library(viridis)
library(GGally)
library(ggstatsplot)
library(devtools)
library(ggbiplot)
library(here)
Median_data<- read.csv(here("data", "processed", "Median_data.csv"))
head(Median_data)

# Create the ggpairs plot
ggpairs(
  Median_data,
  columns = 3:11,
  aes(colour = parasitism, fill = parasitism),  # Map both colour and fill to parasitism
  lower = list(continuous = "smooth"),         # Smooth lines in lower panels
  upper = list(continuous = wrap("cor", method = "kendall"))  # Correlation in upper panels
) +
  scale_colour_manual(
    values = c("Parasite" = alpha("red",0.7), "Host" =alpha("grey",0.7))  # Customize line and point colors
  ) +
  scale_fill_manual(
    values = c("Parasite" = alpha("red",0.7), "Host" =alpha("grey",0.7)) # Customize fill colors
  ) +
  theme_minimal()  # Apply a minimal theme for better visuals


library(patchwork)

# Create the correlation matrix for "Host"
host_corr <- ggcorrmat(
  data = subset(Median_data, parasitism == "Host"),
  type = "np",
  method = "kendall",
  label = TRUE
) +
  scale_fill_viridis_c(option = "viridis", alpha = 0.8) +
  labs(title = "Host")

# Create the correlation matrix for "Parasite"
parasite_corr <- ggcorrmat(
  data = subset(Median_data, parasitism == "Parasite"),
  type = "np",
  method = "kendall",
  label = TRUE
) +
  scale_fill_viridis_c(option = "viridis", alpha = 0.8) +
  labs(title = "Parasite") +
  theme(legend.position = "none") +  # Remove the legend for the second plot
  guides(fill = "none")             # Remove the color scale for the second plot

# Combine the two plots
combined_plot <- host_corr + parasite_corr + 
  plot_annotation(title = "Kendall's Correlation Matrix")&
  theme(plot.title = element_text(size = 20, hjust = 0.5))  

# Display the combined plot
print(combined_plot)



