# Load libraries
library(plotly)
library(here)
library(tidyverse)
library(matrixStats)
library(ggpubr)
library(ggnewscale)  
library(htmlwidgets)
# Load data
PitMembrane_data <- read.csv(here("data", "processed", "PitMembrane_data.csv"))

lucian_model <- function(n_pits, tpm_layer, prob_lp_layer = 0.25) {
  1 - (1 - prob_lp_layer^((tpm_layer + 20) / 40))^n_pits
}

# Create the background plot
number_of_pits <- seq(0, 30e3, 1e2)
pit_membrane_thickness <- seq(0, 1200, 20)
probability_of_large_pore_per_layers <- 0.25
background_plot <- outer(number_of_pits, pit_membrane_thickness, lucian_model)

# Define the range of parasitic plants
ptm_range_plant1 <- seq(338, 1021, 1) # Min to max range
plant1_plot <- outer(number_of_pits, ptm_range_plant1, lucian_model)
mean_ptm_plant1 <- mean(PitMembrane_data$Tpm[PitMembrane_data$parasitism == "Parasite"], na.rm = T)
median_ptm_plant1 <- median(PitMembrane_data$Tpm[PitMembrane_data$parasitism == "Parasite"], na.rm = T)
# Plot
wave <-
  plot_ly(
    x = ~pit_membrane_thickness,
    y = ~number_of_pits,
    z = ~background_plot
  ) %>%
  layout(
    scene = list(
      xaxis = list(title = "Tpm (nm)"),
      yaxis = list(title = "Number of pits in a <br> average vessel"),
      zaxis = list(title = "Embolism risk")
    )
  ) %>%
  add_surface(colorbar = list(title = "Probability")) %>%
  add_surface(
    x = ~ptm_range_plant1, # add parasitic plant layer
    y = ~number_of_pits,
    z = ~ plant1_plot + 0.001,
    colorscale = list(
      list(0, "white"),
      list(0.05, "pink"),
      list(0.1, "red")
    ),
    colorbar = list(title = "Parasite Tpm Probability")
  ) %>%
  add_trace(
    x = rep(mean_ptm_plant1, length(number_of_pits)),
    y = number_of_pits,
    z = background_plot[, which.min(abs(pit_membrane_thickness - mean_ptm_plant1))] + 0.011,
    type = "scatter3d",
    mode = "lines",
    line = list(color = "black", width = 4),
    showlegend = FALSE
  ) %>%
  add_trace(
    x = rep(median_ptm_plant1, length(number_of_pits)),
    y = number_of_pits,
    z = background_plot[, which.min(abs(pit_membrane_thickness - median_ptm_plant1))] + 0.011,
    type = "scatter3d",
    mode = "lines",
    line = list(color = "blue", width = 4),
    showlegend = FALSE
  )



# Position curves at the BACK of the Y-axis
fixed_y <- max(number_of_pits)

# Loop through species

  data <- PitMembrane_data %>%
    filter(parasitism == "Parasite", !is.na(Tpm))
  

    dens <- density(data$Tpm)
    
    # Normalize density height for Z axis
    dens_z <- dens$y / max(dens$y) * max(background_plot) * 0.3
    dens_y <- rep(fixed_y, length(dens$x))
    
    wave <- wave %>%
      add_trace(
        x = dens$x,
        y = dens_y,
        z = dens_z,
        type = "scatter3d",
        mode = "lines",
        line = list(width = 5, color = "firebrick"),
        name = "Parasite Tpm densiity",
        showlegend = TRUE
      )
wave <- wave %>%
  layout(
    scene = list(
      xaxis = list(
        title = "Tpm (nm)",
        backgroundcolor = "lightgray",
        gridcolor = "white",
        zerolinecolor = "white",
        showbackground = TRUE
      ),
      yaxis = list(
        title = "Number of pits in a <br> average vessel",
        backgroundcolor = "lightgray",
        gridcolor = "white",
        zerolinecolor = "white",
        showbackground = TRUE
      ),
      zaxis = list(
        title = "Embolism risk",
        backgroundcolor = "lightgray",
        gridcolor = "white",
        zerolinecolor = "white",
        showbackground = TRUE
      )
    ),
    paper_bgcolor = "white",
    plot_bgcolor = "white"
  )
wave

# Find the part of the density below 400
below <- dens$x <= 600

# Approximate area under the curve up to 400
area_below <- sum(diff(dens$x)[below[-1]] * dens$y[below][-1])
pct_below <- area_below * 100
pct_below-100


plot(dens)

library(viridis)  # for viridis color palette

# Step 1: Simulate pore model
pores_total <- 12000
layers <- 60
dat <- data.frame()

for (l in 1:10) {
  for (k in c(5, 10, 20)) {
    data <- matrix(floor(rnorm(pores_total * layers, mean = 20, sd = 15)), pores_total, layers)
    data1 <- matrix(0, pores_total, layers)
    
    for (j in 1:pores_total) {
      for (i in 1:(layers - 1)) {
        if (data[j, i] - k >= 0) {
          data1[j, i] <- 1
        }
      }
    }
    
    data2 <- data1
    for (j in 1:pores_total) {
      for (i in 2:(layers - 1)) {
        if (data2[j, i - 1] + data2[j, i] == i) {
          data2[j, i] <- data2[j, i - 1] + data2[j, i]
        } else {
          data2[j, i] <- 0
        }
      }
    }
    
    data3 <- data2
    for (j in 1:pores_total) {
      for (i in 1:(layers - 1)) {
        data3[j, i] <- data3[j, i] / data3[j, i]
      }
    }
    
    data4 <- colCounts(data3, na.rm = TRUE)
    data4 <- data.frame(frequency = data4)
    data4$gold_size <- k
    data4$layer <- seq(1:nrow(data4))
    data4$relative_model <- (data4$frequency / pores_total) * 100
    
    dat <- rbind(dat, data4)
  }
}

# Step 2: Average across repetitions
dat2 <- dat %>%
  group_by(gold_size, layer) %>%
  summarise_all(mean)

###Layers estimation
PitMembrane_data$layer <- ceiling((PitMembrane_data$Tpm+20) / 40)

# Step 4: Base plot with points colored by gold_size
wave2 <- ggplot() +
  geom_point(data = dat2, aes(x = layer, y = frequency, color = factor(gold_size)), size = 4) +
  scale_color_manual(
    name = "Pore size (nm)",
    values = c("5" = "#1b9e77", "10" = "#d95f02", "20" = "#7570b3")
  ) +
  new_scale_color() +
  scale_x_continuous("Layer") +
  scale_y_continuous("Modelled number of penetrable pores") +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    aspect.ratio = 1,
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
    legend.position = "right",
    text = element_text(size = 20),
    axis.title.x = element_text(vjust = -1.0),
    axis.title.y = element_text(vjust = +3.0),
    axis.title.y.right = element_text(
      margin = margin(t = 0, r = 0, b = 0, l = 10),
      angle = 90
    )
  )
wave2

# Filter parasite data and compute density
parasite_data <- PitMembrane_data %>%
  filter(parasitism == "Parasite", !is.na(layer))

# Calculate density only if there are enough values
if (nrow(parasite_data) > 1) {
  dens <- density(parasite_data$layer, from = 0, to = max(dat2$layer))
  df_dens <- data.frame(
    layer = dens$x,
    density_scaled = dens$y / max(dens$y) * max(dat2$frequency)
  )
}

# Compute mean and median
mean_layer <- mean(parasite_data$layer, na.rm = TRUE)
median_layer <- median(parasite_data$layer, na.rm = TRUE)

# Add to the plot
wave2 <- wave2 +
  geom_line(
    data = df_dens,
    aes(x = layer, y = density_scaled),
    color = "firebrick",
    size = 1.2
  ) +
  geom_vline(
    xintercept = median_layer,
    color = "blue",
    size = 1.2,
    linetype = "dashed"
  ) +
  geom_vline(
    xintercept = mean_layer,
    color = "black",
    size = 1.2,
    linetype = "dashed"
  )

# Show plot
print(wave2)


parasite_data <- PitMembrane_data %>%
  filter(parasitism == "Parasite", !is.na(layer))
dens2 <- density(parasite_data$layer)
below <- dens2$x <= 20 
area_below <- sum(diff(dens2$x)[below[-1]] * dens2$y[below][-1])
pct_below <- area_below * 100
pct_below

saveWidget(wave, here("outputs", "figs", "interactive_wave_plot.html"))
ggsave(
  filename = "penetrable_pores.png",
  plot = wave2,
  path = here("outputs", "figs"),
  dpi = 600,
  width = 8,
  height = 6
)

reticulate::py_install("kaleido")
plotly::save_image(
  wave,
  here::here("outputs", "figs", "interactive_wave_plot.png"),
  width = 4000,
  height = 3000,
  scale = 2
)
