######################################################################
#
#Victor Sibinelli (victor.sibinelli@usp.br / sibinelli95@gmail.com)
#13/07/2024
# Scrpt 03 -Graphics
#################################################################








#Wave graph base don Kaack 2021.
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
ptm_range_plant1 <- seq(338, 1021, 1)  # Min to max range
plant1_plot <- outer(number_of_pits, ptm_range_plant1, lucian_model)

# Plot 
plot_ly(x = pit_membrane_thickness,
        y = number_of_pits,
        z = ~background_plot) %>%
  layout(
    scene = list(xaxis=list(title="Tpm (nm)"),
                 yaxis=list(title="Number of pits in a <br> average vessel"),
                 zaxis=list(title="Probability of encountering <br> a large pore in a vessel"))) %>% 
  add_surface() %>%
  add_surface(x = ~ptm_range_plant1,#add parasitic plant layer
              y = ~number_of_pits,
              z = ~plant1_plot + 0.011,
              colorscale = list(
                list(0, "pink"),
                list(0.5, "red"),
                list(1, "pink")
              ))