# Load necessary libraries
library(tidyverse)
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

source(here("scripts", "Functions.R"))

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





# Count non-NA numeric observations for each numeric column and group by 'ssp'
wcount <- Wall_data %>% 
  group_by(ssp) %>% 
  summarise(across(where(is.numeric), ~sum(!is.na(.))))

dtopcount <- VesselDiameter_data %>% 
  group_by(indiv, ssp) %>%
  filter(VesselDiameter >= quantile(VesselDiameter, 0.9)) %>%
  ungroup() %>% 
  group_by(ssp) %>% 
  summarise(across(where(is.numeric), ~sum(!is.na(.))))

dcount <- VesselDiameter_data %>%
  group_by(ssp) %>%
  summarise(across(where(is.numeric), ~sum(!is.na(.))))

hdcount <- Hydraulic_data %>%
  group_by(ssp) %>%
  summarise(across(where(is.numeric), ~sum(!is.na(.))))

 pitcount <- PitDiOp_data%>%   group_by(ssp) %>%
   summarise(across(where(is.numeric), ~sum(!is.na(.))))
 
fpitcount <- PitFraction_data %>% group_by(ssp) %>%
  summarise(across(where(is.numeric), ~sum(!is.na(.))))

pitmcount <- PitMembrane_data %>% group_by(ssp) %>% 
  summarise(across(where(is.numeric),~sum(!is.na(.))))

relevel_factors(ls())


w_plot <- Wall_data %>%
  ggplot(aes(x = ssp, y = WallThickness, fill = parasitism)) +
  geom_boxplot(outlier.colour = "red", outlier.shape = 8,
               notch = TRUE, varwidth = TRUE, size=0.8, staplewidth = 0.8, alpha=0.5) +
  geom_violin(color = "black", alpha = 0.01) +
  stat_summary(fun.data = mean_cl_boot, geom = "point", color = "blue", size = 3) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", color = "blue", lwd = 0.8) +
  geom_vline(xintercept = c(2.5, 4.5, 6.5)) +
  scale_fill_manual(values = c("Parasite" = "firebrick", "Host" = "grey")) +
  scale_x_discrete(labels = short_names, guide = guide_axis(angle = 45)) +
  labs(x = "Species", y = "Tvw (µm)") +
  annotate("text", x = c(1.5, 3.5, 5.5, 7.5), y = max(Wall_data$WallThickness) * 1.20,
           label = c("ns","ns","ns","*, #"), size = 8) +
  annotate("text", x = seq(1:8), y = max(Wall_data$WallThickness) * 1.1, 
           label = paste("n=", wcount$WallThickness), size = 8) +
  theme(legend.position = "none", axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  scale_y_continuous(limits = c(0, max(Wall_data$WallThickness) * 1.2), expand = expansion(mult = c(0, 0.1))) +
  guides(fill = "none", color = "none")


d_plot <- VesselDiameter_data %>%
  ggplot(aes(x = ssp, y = VesselDiameter, fill = parasitism)) +
  geom_boxplot(outlier.colour = "red", outlier.shape = 8, 
               notch = TRUE, varwidth = TRUE,size=0.8,staplewidth = 0.5,alpha=0.5) +
  geom_violin(color = "black", alpha = 0.01, adjust = 2) +
  stat_summary(fun.data = mean_cl_boot, geom = "point", color = "blue", size =3) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", color = "blue", lwd = 0.8) +
  geom_vline(xintercept = c(2.5, 4.5, 6.5)) +
  scale_fill_manual(values = c("Parasite" = "firebrick", "Host" = "grey")) +
  scale_x_discrete(labels = short_names, guide = guide_axis(angle = 45)) +
  labs(x = "Species", y = "D (µm)") +
  annotate("text", x = c(1.5, 3.5, 5.5, 7.5), y = max(VesselDiameter_data$VesselDiameter)*1.15, 
           label =c("#","ns","#","*,#"), size =10) +
  annotate("text", x = 1:8, y = max(VesselDiameter_data$VesselDiameter)*1.05, 
           label = paste("n=", dcount$VesselDiameter), size = 6) +
  theme(legend.position = "none", axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  scale_y_continuous(limits = c(0, max(VesselDiameter_data$VesselDiameter) * 1.15), expand = expansion(mult = c(0, 0.1))) +
  guides(fill = "none", color = "none")


dtop_plot <- VesselDiameter_data %>%
  group_by(indiv, ssp) %>%
  filter(VesselDiameter >= quantile(VesselDiameter, 0.9)) %>%
  ungroup() %>%
  ggplot(aes(x = ssp, y = VesselDiameter, fill = parasitism)) +
  geom_boxplot(outlier.colour = "red", outlier.shape = 8,
               notch = TRUE, varwidth = T,size=0.8,staplewidth = 0.5,alpha=0.5) +
  geom_violin(color = "black", alpha = 0.01, adjust = 2) +
  stat_summary(fun.data = mean_cl_boot, geom = "point", color = "blue", size =3) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", color = "blue", lwd = 0.8) +
  geom_vline(xintercept = c(2.5, 4.5, 6.5)) +
  scale_fill_manual(values = c("Parasite" = "firebrick", "Host" = "grey")) +
  scale_x_discrete(labels = short_names, guide = guide_axis(angle = 45)) +
  labs(x = "Species", y = "Dtop (µm)") +
  annotate("text", x = c(1.5, 3.5, 5.5, 7.5), y = max(VesselDiameter_data$VesselDiameter)*1.2,
           label = "*,#", size =10) +
  annotate("text", x = 1:8, y = max(VesselDiameter_data$VesselDiameter)*1.1,
           label = paste("n=", dtopcount$VesselDiameter), size = 6) +
  theme(legend.position = "none", axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  scale_y_continuous(limits = c(0, max(VesselDiameter_data$VesselDiameter) * 1.2), expand = expansion(mult = c(0, 0.1))) +
  guides(fill = "none", color = "none")




# Hydraulic Diameter Plot (hd_plot)
hd_plot <- Hydraulic_data %>%
  ggplot(aes(x = ssp, y = HydraulicDiameter, fill = parasitism)) +
  geom_boxplot(outlier.colour = "red", outlier.shape = 8,
               notch = TRUE, varwidth = TRUE, size = 0.8, staplewidth = 0.5, alpha = 0.5) +
  geom_violin(color = "black", alpha = 0.01, adjust = 2) +
  stat_summary(fun.data = mean_cl_boot, geom = "point", color = "blue", size =3) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", color = "blue", lwd = 0.8) +
  geom_vline(xintercept = c(2.5, 4.5, 6.5)) +
  scale_fill_manual(values = c("Parasite" = "firebrick", "Host" = "grey"), name = "Parasitism") +
  scale_x_discrete(labels = short_names, guide = guide_axis(angle = 45))  +
  labs(x = "Species", y = "Dh (µm)") +
  annotate("text", x = c(1.5, 3.5, 5.5, 7.5),
           y = max(Hydraulic_data$HydraulicDiameter) * 1.2,
           label = c("*,#","#","*,#", "*,#"), size =10) +
  annotate("text", x = 1:8,
           y = max(Hydraulic_data$HydraulicDiameter)*1.15,
           label = paste("n=", hdcount$HydraulicDiameter), size =8 )+
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  guides(fill = "none", color = "none")

 vd_plot <- Hydraulic_data %>%
  ggplot(aes(x = ssp, y = VesselDensity, fill = parasitism)) +
  geom_boxplot(outlier.colour = "red", outlier.shape = 8,
               notch = TRUE, varwidth = TRUE, size = 0.8, staplewidth = 0.5, alpha = 0.5) +
  geom_violin(color = "black", alpha = 0.01, adjust = 2) +
  stat_summary(fun.data = mean_cl_boot, geom = "point", color = "blue", size =3,show.legend = FALSE) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", color = "blue", lwd = 0.8,show.legend = FALSE) +
  geom_vline(xintercept = c(2.5, 4.5, 6.5)) +
  scale_fill_manual(values = c("Parasite" = "firebrick", "Host" = "grey"), name = "Parasitism") +
  scale_x_discrete(labels = short_names, guide = guide_axis(angle = 45)) +
  labs(x = "Species", y = "VD (vessels·mm⁻¹)") +
  annotate("text", x = c(1.5, 3.5, 5.5, 7.5),
           y = max(Hydraulic_data$VesselDensity) * 1.25,
           label = c("ns","ns", "#","ns"), size =10) +
  annotate("text", x = 1:8,
           y = max(Hydraulic_data$VesselDensity) * 1.05,
           label = paste("n=", hdcount$VesselDensity), size =8)+
   theme(legend.position = "none",
         axis.title.x = element_blank(),
        axis.text.x = element_blank(),
         axis.ticks.x = element_blank()) 


# Vessel Fraction Plot (fv_plot)
    fv_plot <- Hydraulic_data %>%
   ggplot(aes(x = ssp, y = VesselFraction, fill = parasitism)) +
   geom_boxplot(outlier.colour = "red", outlier.shape = 8,
                notch = TRUE, varwidth = TRUE, size = 0.8, staplewidth = 0.5, alpha = 0.5) +
   geom_violin(color = "black", alpha = 0.01, adjust = 2) +
   stat_summary(fun.data = mean_cl_boot, geom = "point", color = "blue", size =3,show.legend = F) +
   stat_summary(fun.data = mean_cl_boot, geom = "errorbar", color = "blue", lwd = 0.8,show.legend = F) +
   geom_vline(xintercept = c(2.5, 4.5, 6.5)) +
   scale_fill_manual(values = c("Parasite" = "firebrick", "Host" = "grey"),
                     name = "Parasitism") +
   scale_x_discrete(labels = short_names, guide = guide_axis(angle = 45)) +
   labs(x = "Species", y = "Fv (%)") +
   annotate("text", x = c(1.5, 3.5, 5.5, 7.5),
            y = max(Hydraulic_data$VesselFraction) * 1.25,
            label = c("*,#", "*,#", "#", "*,#"), size =10) +
   annotate("text", x = 1:8,
            y = max(Hydraulic_data$VesselFraction) * 1.1,
            label = paste("n=", hdcount$VesselFraction), size =8)+
      theme(legend.position = "none",
            axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank()) 
    

# Kmax Plot (kmax_plot)
kmax_plot <- Hydraulic_data %>%
  ggplot(aes(x = ssp, y = Kmax, fill = parasitism)) +
  geom_boxplot(outlier.colour = "red", outlier.shape = 8,
               notch = TRUE, varwidth = TRUE, size = 0.8, staplewidth = 0.5, alpha = 0.5) +
  geom_violin(color = "black", alpha = 0.01, adjust = 2) +
  stat_summary(fun.data = mean_cl_boot, geom = "point", color = "blue", size =3) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", color = "blue", lwd = 0.8) +
  geom_vline(xintercept = c(2.5, 4.5, 6.5)) +
  scale_fill_manual(values = c("Parasite" = "firebrick", "Host" = "grey"),
                    name = "Parasitism") +
  scale_x_discrete(labels = short_names, guide = guide_axis(angle = 45)) +
  scale_y_log10(expand = c(0, 0.2)) +
  labs(x = "Species", y = "Kmax (kg·s⁻¹·MPa⁻¹·m⁻¹)",size=8) +
  annotate("text", x = c(1.5, 3.5, 5.5, 7.5),
           y = max(Hydraulic_data$Kmax) * 3,
           label = c("*,#", "*,#", "#", "*,#"), size =10) +
  annotate("text", x = 1:8,
           y = max(Hydraulic_data$Kmax) * 1.5,
           label = paste("n=", hdcount$Kmax), size =8) +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  guides(fill = "none", color = "none")


# Pit Diameter Plot (pd_plot)
pd_plot <- PitDiOp_data %>%
  ggplot(aes(x = ssp, y = PitDiameter, fill = parasitism)) +
  geom_boxplot(outlier.colour = "red", outlier.shape = 8,
               notch = TRUE, varwidth = TRUE, size = 0.8, 
               staplewidth = 0.5, alpha = 0.5) +
  geom_violin(color = "black", alpha = 0.01, adjust = 2) +
  stat_summary(fun.data = mean_cl_boot, geom = "point", color = "blue", size =3) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", color = "blue", lwd = 0.8) +
  geom_vline(xintercept = c(2.5, 4.5, 6.5)) +
  scale_fill_manual(values = c("Parasite" = "firebrick", "Host" = "grey"),
                    name = "Parasitism") +
  scale_x_discrete(labels = short_names, guide = guide_axis(angle = 45)) +
  labs(x = "Species", y = "Dpit (µm)") +
  annotate("text", x = c(1.5, 3.5, 5.5, 7.5),
           y = max(PitDiOp_data$PitDiameter, na.rm = TRUE) * 1.25,
           label = rep("*,#", 4), size =10) +
  annotate("text", x = 1:8,
           y = max(PitDiOp_data$PitDiameter, na.rm = TRUE) * 1.1,
           label = paste("n=", pitcount$PitDiameter), size = 8) +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  guides(fill = "none", color = "none") 

# Pit Opening Plot (po_plot)
po_plot <- PitDiOp_data %>%
  ggplot(aes(x = ssp, y = PitOpening, fill = parasitism)) +
  geom_boxplot(outlier.colour = "red", outlier.shape = 8,
               notch = TRUE, varwidth = TRUE, size = 0.8, 
               staplewidth = 0.5, alpha = 0.5) +
  geom_violin(color = "black", alpha = 0.01, adjust = 2) +
  stat_summary(fun.data = mean_cl_boot, geom = "point", color = "blue", size =3,show.legend = FALSE) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", color = "blue", lwd = 0.8,show.legend = FALSE) +
  geom_vline(xintercept = c(2.5, 4.5, 6.5)) +
  scale_fill_manual(values = c("Parasite" = "firebrick", "Host" = "grey"),
                    name = "Parasitism") +
  scale_x_discrete(labels = short_names, guide = guide_axis(angle = 45)) +
  labs(x = "Species", y = "Dpa (µm)") +
  annotate("text", x = c(1.5, 3.5, 5.5, 7.5),
           y = max(PitDiOp_data$PitOpening, na.rm = TRUE) * 1.25,
           label = c("#", "ns", "*,#", "*,#"), size =10) +
  annotate("text", x = 1:8,
           y = max(PitDiOp_data$PitOpening, na.rm = TRUE) * 1.1,
           label =paste("n=", pitcount$PitOpening), size = 8) 
  

# Pit Fraction Plot (fp_plot)
fp_plot <- PitFraction_data %>%
  ggplot(aes(x = ssp, y = PitFraction, fill = parasitism)) +
  geom_boxplot(outlier.colour = "red", outlier.shape = 8,
               notch = TRUE, varwidth = TRUE, size = 0.8, 
               staplewidth = 0.5, alpha = 0.5) +
  geom_violin(color = "black", alpha = 0.01, adjust = 2) +
  stat_summary(fun.data = mean_cl_boot, geom = "point", color = "blue", size =3,show.legend = FALSE) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", color = "blue", lwd = 0.8,show.legend = FALSE) +
  geom_vline(xintercept = c(2.5, 4.5, 6.5)) +
  scale_fill_manual(values = c("Parasite" = "firebrick", "Host" = "grey"),
                    name = "Parasitism") +
  scale_x_discrete(labels = short_names, guide = guide_axis(angle = 45)) +
  labs(x = "Species", y = "Fp (%)") +
  annotate("text", x = c(1.5, 3.5, 5.5, 7.5),
           y = max(PitFraction_data$PitFraction, na.rm = TRUE) * 1.25,
           label = c("#", "ns", "*,#", "*,#"), size =10) +
  annotate("text", x = 1:8,
           y = max(PitFraction_data$PitFraction, na.rm = TRUE) * 1.1,
           label = paste("n=", fpitcount$PitFraction), size =8)
 



pcd_plot <- PitMembrane_data %>%
  ggplot(aes(x = ssp, y = pcd, fill = parasitism)) +
  geom_boxplot(outlier.colour = "red", outlier.shape = 8,
               notch = TRUE, varwidth = TRUE, size = 0.8, 
               staplewidth = 0.5, alpha = 0.5) +
  geom_violin(color = "black", alpha = 0.01, adjust = 2) +
  stat_summary(fun.data = mean_cl_boot, geom = "point", color = "blue", size =3,show.legend = FALSE) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", color = "blue", lwd = 0.8,show.legend = FALSE) +
  geom_vline(xintercept = c(2.5, 4.5, 6.5)) +
  scale_fill_manual(values = c("Parasite" = "firebrick", "Host" = "grey"),
                    name = "Parasitism") +
  scale_x_discrete(labels = short_names, guide = guide_axis(angle = 45)) +
  labs(x = "Species", y = "Hpit (nm)") +
  annotate("text", x = c(1.5, 3.5, 5.5, 7.5),
           y = max(PitMembrane_data$pcd, na.rm = TRUE) * 1.25,
           label = c("*,#"), size =10) +
  annotate("text", x = 1:8,
           y = max(PitMembrane_data$pcd, na.rm = TRUE) * 1.1,
           label = paste("n=", pitmcount$pcd), size =8)
  
  
  
  


Tpm_plot <- PitMembrane_data %>%
  ggplot(aes(x = ssp, y = Tpm, fill = parasitism)) +
  geom_boxplot(outlier.colour = "red", outlier.shape = 8,
               notch = TRUE, varwidth = TRUE, size = 0.8, 
               staplewidth = 0.5, alpha = 0.5) +
  geom_violin(color = "black", alpha = 0.01, adjust = 2) +
  stat_summary(fun.data = mean_cl_boot, geom = "point", color = "blue", size =3,show.legend = FALSE) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", color = "blue", lwd = 0.8,show.legend = FALSE) +
  geom_vline(xintercept = c(2.5, 4.5, 6.5)) +
  scale_fill_manual(values = c("Parasite" = "firebrick", "Host" = "grey"),
                    name = "Parasitism") +
  scale_x_discrete(labels = short_names, guide = guide_axis(angle = 45)) +
  labs(x = "Species", y = "Tpm (nm)",size=8) +
  annotate("text", x = c(1.5, 3.5, 5.5, 7.5),
           y = max(PitMembrane_data$Tpm, na.rm = TRUE) * 1.25,
           label = c("*,#"), size =10) +
  annotate("text", x = 1:8,
           y = max(PitMembrane_data$Tpm, na.rm = TRUE) * 1.1,
           label = paste("n=", pitmcount$Tpm), size =8) 


# Define axis_scale theme to apply to all plots
axis_scale <- theme(
  axis.title.x = element_text(size = 25), 
  axis.title.y = element_text(size = 35), 
  axis.text.x = element_text(size = 20), 
  axis.text.y = element_text(size = 20),
  legend.text = element_text(size = 30), 
  legend.title = element_text(size = 30),
  legend.position = "bottom",
  plot.margin = margin(8, 8, 8, 8)
)

# Apply axis_scale to all the plots
d_plot <- d_plot + axis_scale + theme(axis.text.x = element_blank(),
                                      axis.title.x=element_blank())
dtop_plot <- dtop_plot + axis_scale + theme(axis.text.x = element_blank(),
                                            axis.title.x=element_blank())
hd_plot <- hd_plot + axis_scale + theme(axis.text.x = element_blank(),
                                        axis.title.x=element_blank())
vd_plot <- vd_plot + axis_scale + theme(axis.text.x = element_blank(),
                                        axis.title.x=element_blank())
fv_plot <- fv_plot + axis_scale + theme(axis.text.x = element_blank(),
                                        axis.title.x=element_blank())
kmax_plot <- kmax_plot + axis_scale + theme(axis.text.x = element_blank(),
                                            axis.title.x=element_blank())
w_plot <- w_plot + axis_scale+ theme(axis.text.x = element_blank(),                                        
                                     axis.title.x=element_blank())
pd_plot <- pd_plot + axis_scale+ theme(axis.text.x = element_blank(),                                        
                                       axis.title.x=element_blank())
po_plot <- po_plot + axis_scale
fp_plot <- fp_plot + axis_scale
pcd_plot <- pcd_plot + axis_scale
Tpm_plot <- Tpm_plot + axis_scale

combined_plot <- 
  ( hd_plot + dtop_plot) /
  (pd_plot+ kmax_plot) /
  ( pcd_plot + Tpm_plot) +
  plot_layout(guides = "collect") &  plot_annotation(tag_levels = 'a')&
  theme(plot.tag = element_text(face = 'bold', size = 35))&
  theme(
    legend.position = "bottom"
  )

# Display the combined plot
print(combined_plot)
combined_plot2 <- (d_plot + vd_plot)/
   (fv_plot+w_plot) /
  (po_plot+fp_plot) + plot_layout(guides = "collect") &  plot_annotation(tag_levels = 'a')&
  theme(plot.tag = element_text(face = 'bold', size = 35))&
  theme(
    legend.position = "bottom"
  )
print(combined_plot2)



ggsave(file=here("outputs","figs","Fig4.png"),
       combined_plot,units = "in",dpi = 600,height = 20,width = 20)

ggsave(file=here("outputs","figs","trait_violin_plot2.png"),
       combined_plot2,units = "in",dpi = 600,height = 20,width = 20)


