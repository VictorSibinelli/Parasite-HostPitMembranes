vdata
boxplot(vdata$VD~vdata$ssp)
boxplot(vdensity~ssp, data=vadata)


HydraulicD <-  (tapply(((vdata$VD*vdata$Circ.)^4),INDEX = vdata$label,FUN=sum)/
  (tapply(vdata$VD, INDEX = vdata$label, FUN=length)))^(1/4)
HydraulicD

vdens <- tapply(vadata$vdensity,INDEX = vadata$label,FUN = mean)
vdens
data.frame(Label = names(HydraulicD),
           HydraulicD = HydraulicD,
           vdensity = vdens[match(names(HydraulicD), names(vdens))])


