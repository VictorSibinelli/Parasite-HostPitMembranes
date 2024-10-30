######################################################################
#
# Victor Sibinelli (victor.sibinelli@usp.br / sibinelli95@gmail.com)
# 13/07/2024
# Scrpt 04.3 -Graphics- Pits
#################################################################
library(here)
source(here("scripts", "Functions.R"))
wdata <- read.csv(here("data", "processed", "wdata.csv"))
WT_AIC <- read.csv(here("outputs", "tables", "Wdata_AIC.csv"))
CI95 <- read.csv(here("outputs", "tables", "WTCI95.csv"))
WT_MC <- read.csv(here("outputs", "tables", "WT_MonteCarlo_results.csv"))
CI_95 <- read.csv(here("outputs", "tables", "WT_MonteCarlo_CI95.csv"))

# Set the desired order for the groups
desired_order <- c("Parasite", "Host",
                   "Psittacanthus robustus", "Vochysia thyrsoidea",
                   "Phoradendron perrotettii", "Tapirira guianensis",
                   "Struthanthus rhynchophyllus", "Tipuana tipu",
                   "Viscum album", "Populus nigra")

relevel_factors(ls())
# Reverse the order in 'desired_order' to display top-down
CI95$Group <- factor(CI95$Group, levels = rev(desired_order))
CI_95$Group <- factor(CI_95$Group, levels = rev(desired_order))
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


########for LMEM
# Create the plot
library("viridis")


g <- ggplot(wdata, aes(x = ssp, y = wthickness, fill = parasitism)) +
  geom_jitter(aes(color = label),
              size = 1, alpha = 0.4, position = position_jitter(width = 0.3)) +
  geom_boxplot(color = "black", alpha = 1, position = position_dodge(width = 0.8)) +
  scale_fill_manual(
    values = c("Parasite" = "firebrick", "Host" = "grey"),
    name = "Parasitism"
  ) +
  scale_color_viridis_d(option = "D") +
 geom_vline(xintercept = c(2.5,4.5,6.5)) +
  theme_classic() +
  scale_x_discrete(labels = short_names) +
  labs(title = "Vessel Wall Thickness",
       x = "Species",
       y = "Vessel Wall Thickness (µm)") +
  annotate("text", x = seq_along(unique(wdata$ssp)),
           y = max(wdata$wthickness) + 5, label = c("A","A","A","A","A","B","A"), size = 6)+
  theme(legend.position = "right") +  # Legend on the right
  guides(color = "none")  # Remove the legend for `ssp`
g
ggsave(here("outputs","figs","vessel_wall_thickness_plot.png"), plot = g, dpi = 600, width = 10, height = 7)


# Plot with alternating colors
CI95 %>% 
  ggplot(aes(Group, Estimate)) +
  geom_point(size = 4, aes(color = Group)) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2) +
  coord_flip() +
  labs(title = "Estimates and 95% Confidence Intervals",
       x = "Effect", y = "Estimate") +
  scale_color_manual(
    values = rep(c("black", "firebrick"), length.out = nlevels(CI95$Group))  # Alternate colors based on levels
  ) +
  theme_minimal() +
  theme(legend.position = "none")  # Remove legend if not needed








##################For ressampling
boot_sspWT_long <- boot_sspWT %>%
  pivot_longer(cols = everything(), 
               names_to = "species", 
               values_to = "wthickness") %>%
  mutate(parasitism = case_when(
    species %in% c("Psittacanthus robustus", 
                   "Phoradendron perrotettii", 
                   "Struthanthus rhynchophyllus", 
                   "Viscum album") ~ "Parasite",
    TRUE ~ "Host"
  ))
WT_AIC
WT_MC

g <- ggplot(wdata, aes(x = ssp, y = wthickness, fill = parasitism)) +
  geom_jitter(aes(color = label),
              size = 1, alpha = 0.4, position = position_jitter(width = 0.3)) +
  geom_boxplot(color = "black", alpha = 1, position = position_dodge(width = 0.8)) +
  scale_fill_manual(
    values = c("Parasite" = "firebrick", "Host" = "grey"),
    name = "Parasitism"
  ) +
  scale_color_viridis_d(option = "D") +
  geom_vline(xintercept = c(2.5,4.5,6.5)) +
  theme_classic() +
  scale_x_discrete(labels = short_names) +
  labs(title = "Vessel Wall Thickness",
       x = "Species",
       y = "Vessel Wall Thickness (µm)") +
  annotate("text", x = seq_along(unique(wdata$ssp)),
           y = max(wdata$wthickness) + 5, label = c("A","A","A","A","A","A","A"), size = 6)+
  annotate("text", x = c(1.5,3.5,5.5,7.5),
           y = max(wdata$wthickness) + 5, 
           label = "", 
           size = 6)+
  theme(legend.position = "right",
        axis.text.x = element_text(size = 12),        # X-axis tick labels size
        axis.text.y = element_text(size = 12)         # Y-axis tick labels size
  ) +  # Legend on the right
  guides(color = "none")  # Remove the legend for `ssp`
g



CI_95 %>% 
  ggplot(aes(Group, Estimate)) +
  geom_point(size = 4, aes(color = Group)) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2) +
  coord_flip() +
  labs(title = "Estimates and 95% Confidence Intervals",
       x = "Effect", y = "Estimate") +
  scale_color_manual(
    values = rep(c("firebrick", "black"), length.out = nlevels(CI_95$Group))  # Alternate colors based on levels
  ) +
  theme_minimal() +
  theme(legend.position = "none")  # Remove legend if not needed





