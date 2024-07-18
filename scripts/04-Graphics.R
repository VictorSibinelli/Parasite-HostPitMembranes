######################################################################
#
#Victor Sibinelli (victor.sibinelli@usp.br / sibinelli95@gmail.com)
#13/07/2024
# Scrpt 03 -Graphics
#################################################################
source("scripts/03-DataAnalysis.R")






###Pit graphs

#grouped parasite and host pit membrane
png(here("outputs","figs","p-h_pitmembrane.png"), width = 9, height = 7, units = "in", res = 600)
boxplot(data=pitdata_clean,pitavg*1000~parasitism,
        ylab="Pit Membrane Thickness (nm)", xlab="", tcl=T, xaxt="n", at=c(1,3),
        col=c("grey","firebrick"), ylim=c(0,1000), cex.lab=1.5, cex.axis=1.5,whisklwd=3,staplelwd=3)
legend("topleft",cex=2, bty="n", fill=c("grey", "firebrick"), legend = c("Host", "Parasite"))
text(x=c(2), y=c(0.3120534,0.5938832)*1000, "--", cex=2, col="black")
segments(x0=2,y0=0.3120534*1000,x1=2,y1=0.5938832*1000,lwd=3)
text(x=2.2,y=(0.2818298/2+0.3120534)*1000,"0.28", cex=2)
text(x=c(1,3),y=c(2,1.2)-0.1, c("A","B"), cex=1.5)
dev.off()

#grouped parasite host pcd

png(here("outputs","figs","p-h_pitchamber.png"), width = 9, height = 7, units = "in", res = 600)
par(mar = c(5, 6, 4, 2) + 0.1, mgp = c(3, 1, 0))
boxplot(data=pitdata_clean, pcd*1000~parasitism, na.rm=T, las=2,ylim=c(0,2000),cex.lab=1.5,at=c(1,3),
        ylab="", xlab="", tcl=T, xaxt="n",cex.axis=1.5,
        col=c("grey","firebrick"),whisklwd=3,staplelwd=3, boxwex=0.5)
legend(cex=2,"topright", bty="n", fill=c("grey", "firebrick"), legend = c("Host", "Parasite"))
segments(x0=2,y0=0.5342017*1000,x1=2,y1=0.8891008*1000, lwd=3)
text(cex=2,x=2,
     y=1000*c(mean(pitdata_clean$pcd[pitdata_clean$parasitism=="h"], na.rm = T),mean(pitdata_clean$pcd[pitdata_clean$parasitism=="p"], na.rm = T)), "--")
text(x=2.2, y=(0.3608233/2+0.5342017)*1000,0.36, cex=2)
text(x=c(1,3),y=c(2,1.2), c("A","B"), cex=1.5)
title(ylab="Pit Chamber depth (nm)", line=4, cex.lab=1.5)
dev.off()


#pairwise pit membrane thickness and pcd
png(here("outputs","figs","pairwise_pit.png"), width =24, height = 9, units = "in", res = 600)
par(mar = c(7, 6, 1, 2) + 0.1, mgp = c(3, 1, 0),mfrow=c(1,2))

boxplot(data=pitdata_clean, pitavg*1000~ssp, na.rm=T,las=2,cex=1,
        ylab="", xlab=NA, tcl=T,xaxt="n",
        col=c("firebrick","grey"), at=c(1,3,5,7,9,11,13,15),
        ylim=c(0,1100), cex.lab=1.5,cex.axis=1.5,whisklwd=2,staplelwd=2)
title(ylab="Pit Membrane Thickness (nm)", line=4, cex.lab=1.5)
text(x=c(1,3,5,7,9,11,13,15), y=-150, c("P. robustus", "V. thyrsoidea","P. perrottetti", "T. guianensis", "S. rhynchophyllus", "T. tipu", "V. album","P. nigra"),
     xpd=NA  , cex=1.5,srt = 35, col="black", adj=0.50, font = 3)
axis(1,at=c(1,3,5,7,9,11,13,15), labels = NA)
text(x=c(1,3,5,7,9,11,13,15),y= 100, labels = c("A","B"), cex=2)
abline(v= c(4,8,12),lty=2)
text(x=c(2,2,6,6,10,10,14,14), y=pit_means*1000, "-", cex=2, col="black")
segments(x0=c(2,6,10,14),y0=pit_mean_diff$SP1mean*1000,y1=pit_mean_diff$SP2mean*1000,lwd=2)
text(x=c(2,6,10,14)+0.3, y=(pit_mean_diff$SP1mean + pit_mean_diff$SP2mean) * 500,round(pit_mean_diff$MeanDifference,3)*1000, cex=1.2)


boxplot(data=pitdata_clean, pcd*1000~ssp, na.rm=T,las=2,cex=1,
        ylab="", xlab=NA, tcl=T,xaxt="n",
        col=c("firebrick","grey"), at=c(1,3,5,7,9,11,13,15),
        ylim=c(0,2000), cex.lab=1.5,cex.axis=1.5,whisklwd=2,staplelwd=2)
title(ylab="Pit Chamber depth (nm)", line=4, cex.lab=1.5)
text(x=c(1,3,5,7,9,11,13,15), y=-300, c("P. robustus", "V. thyrsoidea","P. perrottetti", "T. guianensis", "S. rhynchophyllus", "T. tipu", "V. album","P. nigra"),
     xpd=NA  , cex=1.5,srt = 35, col="black", adj=0.50, font = 3)
axis(1,at=c(1,3,5,7,9,11,13,15), labels = NA)
text(x=c(1,3,5,7,9,11,13,15),y= 100, labels = c("A","B"), cex=1.5)
abline(v= c(4,8,12),lty=2)
text(x=c(2,2,6,6,10,10,14,14), y=pcd_measn*1000, "-", cex=2, col="black")
segments(x0=c(2,6,10,14),y0=pcd_measn[seq(1,8,by=2)]*1000,y1=pcd_measn[seq(2,8,by=2)]*1000,lwd=2)
text(x=c(2,6,10,14)+0.3, y=(pcd_measn[seq(1,8,by=2)]+(pcd_results$MeanDifference/2))*1000,round(pcd_results$MeanDifference,3)*1000, cex=1.2)
dev.off()

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
