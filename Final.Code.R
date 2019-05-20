Code for manuscript titled: Alternative designs and tropical tree seedling growth performance landscapes

Authors: Samantha J. Worthy, Daniel C. Laughlin, María Natalia Umaña, Caicai Zhang, Luxiang Lin, Min Cao, and Nathan G. Swenson

# Fixing the trait data
Read in Trait Data

trait.data=read.csv("growth.trait.data.csv", header=T, row.names=1)

# First remove the tag number from the end of plot name.

trait.data$Tag=gsub("_.*", "", trait.data$Tag)

# Calculate LMA by taking 1/SLA

trait.data[,11]=1/trait.data[,2]
colnames(trait.data)[11]="lma"

# Log and scale traits

traits.scaled=apply(log(traits), MARGIN=2, scale)

# Read in Soil Data

soil.data=read.csv("soil.data.csv", header=T)

# Log and scale the soil data

soil.data.scaled=apply(log(soil.data), MARGIN=2, scale)

# Read in Light Data

light.data=read.csv("LightforSam.csv", header=T)

# Remove the underscore after the plot number

light.data$plot=gsub("_", "", light.data$plot)

# Log and scale light data

light.data[,3]=scale(log(light.data$perc.canopy.open))

# Merge trait and light data by plot name

trait.light=merge(light.data, trait.data, by.x="plot", by.y="Tag")

# Log and scale LMA

trait.light[,14]=scale(log(trait.light$lma))

# Log and scale LMF

trait.light[,15]=scale(log(trait.light$lmf))

# log and scale RGR.h, first add 1 to RGR.h

trait.light[,16]=trait.light$rgr.h+1
trait.light[,17]=scale(log(trait.light$new.rgr.h))

# Run pca on scaled soil data
## Put scores for axes 1-3 into a new object
## Make a biplot of the pca output; Figure S11; Table S1

pc.soil=princomp(soil.data.scaled)
summary(pc.soil)
pc.soil.scores=pc.soil$scores[,1:3]
loadings.pc=with(pc.soil, unclass(loadings))
biplot(pc.soil,xlabs=c(rep("o",nrow(soil.data.scaled))),col=c(1,4),cex=0.7)

# Merge pc.soil.scores with trait.light.scaled

all.data=merge(pc.soil.scores.df, trait.light, by.x="V4", by.y="plot")

# Slopes for two-way interaction models

lma.light.output=lma.light[2,2]+lma.light[5,2]*final.data$log.light
rmf.light.output=rmf.light[2,2]+rmf.light[5,2]*final.data$log.light

# Plotting slopes against environment

plot(final.data$log.light, lma.light.output)
plot(final.data$log.light, rmf.light.output)

# Slopes for three-way interaction models

lma.rmf.light.output=coef(summary(lma.rmf.light))[2,1]+coef(summary(lma.rmf.light))[5,1]*all.data$log.rmf+coef(summary(lma.rmf.light))[6,1]*all.data$log.light+coef(summary(lma.rmf.light))[8,1]*(all.data$log.rmf*all.data$log.light)
lma.rmf.pc1.output=coef(summary(lma.rmf.pc1))[2,1]+coef(summary(lma.rmf.pc1))[5,1]*all.data$log.rmf+coef(summary(lma.rmf.pc1))[6,1]*all.data$Comp.1+coef(summary(lma.rmf.pc1))[8,1]*(all.data$log.rmf*all.data$Comp.1)
lma.rmf.pc2.output=coef(summary(lma.rmf.pc2))[2,1]+coef(summary(lma.rmf.pc2))[5,1]*all.data$log.rmf+coef(summary(lma.rmf.pc2))[6,1]*all.data$Comp.2+coef(summary(lma.rmf.pc2))[8,1]*(all.data$log.rmf*all.data$Comp.2)
lma.mean.thick.pc1.output=coef(summary(lma.mean.thick.pc1))[2,1]+coef(summary(lma.mean.thick.pc1))[5,1]*all.data$log.mean.thick+coef(summary(lma.mean.thick.pc1))[6,1]*all.data$Comp.1+coef(summary(lma.mean.thick.pc1))[8,1]*(all.data$log.mean.thick*all.data$Comp.1)
thick.ssl.pc2.output=coef(summary(thick.ssl.pc2))[2,1]+coef(summary(thick.ssl.pc2))[5,1]*all.data$log.ssl+coef(summary(thick.ssl.pc2))[6,1]*all.data$Comp.2+coef(summary(thick.ssl.pc2))[8,1]*(all.data$log.ssl*all.data$Comp.2)
thick.lar1.pc2.output=coef(summary(thick.lar1.pc2))[2,1]+coef(summary(thick.lar1.pc2))[5,1]*all.data$log.lar1+coef(summary(thick.lar1.pc2))[6,1]*all.data$Comp.2+coef(summary(thick.lar1.pc2))[8,1]*(all.data$log.lar1*all.data$Comp.2)
thick.smf.pc1.output=coef(summary(thick.smf.pc1))[2,1]+coef(summary(thick.smf.pc1))[5,1]*all.data$log.smf+coef(summary(thick.smf.pc1))[6,1]*all.data$Comp.1+coef(summary(thick.smf.pc1))[8,1]*(all.data$log.smf*all.data$Comp.1)
thick.rmf.pc1.output=coef(summary(thick.rmf.pc1))[2,1]+coef(summary(thick.rmf.pc1))[5,1]*all.data$log.rmf+coef(summary(thick.rmf.pc1))[6,1]*all.data$Comp.1+coef(summary(thick.rmf.pc1))[8,1]*(all.data$log.rmf*all.data$Comp.1)
lar.lmf.pc2.output=coef(summary(lar.lmf.pc2))[2,1]+coef(summary(lar.lmf.pc2))[5,1]*all.data$log.lmf+coef(summary(lar.lmf.pc2))[6,1]*all.data$Comp.2+coef(summary(lar.lmf.pc2))[8,1]*(all.data$log.lmf*all.data$Comp.2)
ssl.smf.light.output=coef(summary(ssl.smf.light))[2,1]+coef(summary(ssl.smf.light))[5,1]*all.data$log.smf+coef(summary(ssl.smf.light))[6,1]*all.data$log.light+coef(summary(ssl.smf.light))[8,1]*(all.data$log.smf*all.data$log.light)

# Plotting slopes against environment

plot((all.data$log.rmf*all.data$log.light),lma.rmf.light.output)
plot((all.data$log.rmf*all.data$Comp.1),lma.rmf.pc1.output)
plot((all.data$log.rmf*all.data$Comp.2),lma.rmf.pc2.output)
plot((all.data$log.mean.thick*all.data$Comp.1),lma.mean.thick.pc1.output)
plot((all.data$log.ssl*all.data$Comp.2),thick.ssl.pc2.output)
plot((all.data$log.lar1*all.data$Comp.2),thick.lar1.pc2.output)
plot((all.data$log.smf*all.data$Comp.1),thick.smf.pc1.output)
plot((all.data$log.rmf*all.data$Comp.1),thick.rmf.pc1.output)
plot((all.data$log.lmf*all.data$Comp.2), lar.lmf.pc2.output)
plot((all.data$log.smf*all.data$log.light),ssl.smf.light.output)

# Plots of two-way interactions; Figure 2

library(visreg)
par(mfrow=c(2,2))
visreg2d(lma.light.model, "log.lma", "log.light", plot.type="image", xlab="Leaf Mass per Area", ylab="Light", main="")
visreg2d(lma.light.model, "log.lma", "log.light", plot.type="persp", ylab="Light", zlab="\nRelative Growth Rate", xlab="Leaf Mass per Area", main="Relative Growth Rate = f(Leaf Mass per Area*Light)",cex.main=1,nn=99, cex.lab=1,lwd=0.5, border="grey40", col=adjustcolor("blue",alpha.f=.5))
visreg2d(rmf.light.model, "log.rmf", "log.light", plot.type="image", xlab="Root Mass Fraction", ylab="Light", main="")
visreg2d(rmf.light.model, "log.rmf", "log.light", plot.type="persp", ylab="Light", zlab="\nRelative Growth Rate", xlab="Root Mass Fraction", main="Relative Growth Rate = f(Root Mass Fraction*Light)",cex.main=1,nn=99, cex.lab=1,lwd=0.5, border="grey40", col=adjustcolor("blue",alpha.f=.5))

# Conceptual model, Figure 1

library(fields)

grid.l <- list(seq(-3,3,length=100), seq(-3,3,length=100))

te <- make.surface.grid(grid.l)

par(mfrow=c(1,2),mar=c(2,2,2,2))

te.linear.pred1 <- as.surface(te, 10*te[,1]*te[,2]*-2 + te[,1]*te[,2] + te[,2]*-2 + te[,1]*-2 + te[,1] + te[,2] + -2)

persp(te.linear.pred1$x,te.linear.pred1$y,te.linear.pred1$z, theta = 30, phi = 20,

      xlab="\nTrait X",ylab="\nTrait Y",zlab="\nPerformance",

      main="(a) Low end of environmental gradient",cex.lab=1,lwd=0.5, border="grey40", col=adjustcolor("blue",alpha.f=.5))


te.linear.pred2 <- as.surface(te, 10*te[,1]*te[,2]*2 + te[,1]*te[,2] + te[,2]*2 + te[,1]*2 + te[,1] + te[,2] + 2)

persp(te.linear.pred2$x,te.linear.pred2$y,te.linear.pred2$z, theta = 30, phi = 20,

      xlab="\nTrait X",ylab="\nTrait Y",zlab="\nPerformance",

      main="(b) High end of environmental gradient",cex.lab=1,lwd=0.5, border="grey40", col=adjustcolor("blue",alpha.f=.5))

results.lma.rmf.light <- readRDS("LMA.RMF.Light.jags_results.RData")
lma.rmf.light=read.csv("lma.rmf.light.Parameters.csv", header=T, row.names=1)
gelman.diag(results.lma.rmf.light)
lma.rmf.light.output=output[2,2]+output[3,2]*final.data$log.rmf+output[1,2]*final.data$log.light+output[7,2]*(final.data$log.rmf*final.data$log.light)
plot((final.data$log.rmf*final.data$log.light),lma.rmf.light.output)

results.lar.lmf.pc2 <- readRDS("lar.lmf.pc2.jags_results.RData")
lar.lmf.pc2=read.csv("lar.lmf.pc2.rgr.Parameters.csv", header=T, row.names=1)
gelman.diag(results.lar.lmf.pc2)
lar.lmf.pc2.output=lar.lmf.pc2[2,2]+lar.lmf.pc2[3,2]*final.data$log.lmf+lar.lmf.pc2[1,2]*final.data$Comp.2+lar.lmf.pc2[7,2]*(final.data$log.lmf*final.data$Comp.2)
plot((final.data$log.lmf*final.data$Comp.2),lar.lmf.pc2.output)

results.lma.rmf.pc1 <- readRDS("LMA.RMF.PC1.jags_results.RData")
lma.rmf.pc1=read.csv("lma.rmf.pc1.Parameters.csv", header=T, row.names=1)
gelman.diag(results.lma.rmf.pc1)
lma.rmf.pc1.output=lma.rmf.pc1[2,2]+lma.rmf.pc1[3,2]*final.data$log.rmf+lma.rmf.pc1[1,2]*final.data$Comp.1+lma.rmf.pc1[7,2]*(final.data$log.rmf*final.data$Comp.1)
plot((final.data$log.rmf*final.data$Comp.1),lma.rmf.pc1.output)

results.lma.rmf.pc2 <- readRDS("LMA.RMF.PC2.jags_results.RData")
lma.rmf.pc2=read.csv("lma.rmf.pc2.rgr.Parameters.csv", header=T, row.names=1)
gelman.diag(results.lma.rmf.pc2)
lma.rmf.pc2.output=lma.rmf.pc2[2,2]+lma.rmf.pc2[3,2]*final.data$log.rmf+lma.rmf.pc2[5,2]*final.data$Comp.2+lma.rmf.pc2[7,2]*(final.data$log.rmf*final.data$Comp.2)
plot((final.data$log.rmf*final.data$Comp.2),lma.rmf.pc2.output)

results.lma.smf.light <- readRDS("LMA.SMF.light.jags_results.RData")
lma.smf.light=read.csv("lma.smf.light.rgr.Parameters.csv", header=T, row.names=1)
gelman.diag(results.lma.smf.light)
lma.smf.light.output=lma.smf.light[2,2]+lma.smf.light[3,2]*final.data$log.smf+lma.smf.light[1,2]*final.data$log.light+lma.smf.light[7,2]*(final.data$log.smf*final.data$log.light)
plot((final.data$log.smf*final.data$log.light),lma.smf.light.output)

results.mean.thick.rmf.pc1 <- readRDS("mean.thick.rmf.pc1.jags_results.RData")
mean.thick.rmf.pc1=read.csv("mean.thick.rmf.pc1.rgr.Parameters.csv", header=T, row.names=1)
gelman.diag(results.mean.thick.rmf.pc1 )
mean.thick.rmf.pc1.output=mean.thick.rmf.pc1[2,2]+mean.thick.rmf.pc1[3,2]*final.data$log.rmf+mean.thick.rmf.pc1[1,2]*final.data$Comp.1+mean.thick.rmf.pc1[7,2]*(final.data$log.rmf*final.data$Comp.1)
plot((final.data$log.rmf*final.data$Comp.1),mean.thick.rmf.pc1.output)

results.mean.thick.smf.pc1 <- readRDS("THICK.SMF.pc1.jags_results.RData")
mean.thick.smf.pc1=read.csv("thick.smf.pc1.rgr.Parameters.csv", header=T, row.names=1)
gelman.diag(results.mean.thick.smf.pc1)
mean.thick.smf.pc1.output=mean.thick.smf.pc1[2,2]+mean.thick.smf.pc1[3,2]*final.data$log.smf+mean.thick.smf.pc1[1,2]*final.data$Comp.1+mean.thick.smf.pc1[7,2]*(final.data$log.smf*final.data$Comp.1)
plot((final.data$log.smf*final.data$Comp.1),mean.thick.smf.pc1.output)

results.mean.thick.ssl.pc2 <- readRDS("thick.ssl.pc2.jags_results.RData")
mean.thick.ssl.pc2=read.csv("thick.ssl.pc2.rgr.Parameters.csv", header=T, row.names=1)
gelman.diag(results.mean.thick.ssl.pc2)
mean.thick.ssl.pc2.output=mean.thick.ssl.pc2[2,2]+mean.thick.ssl.pc2[3,2]*final.data$log.ssl+mean.thick.ssl.pc2[1,2]*final.data$Comp.2+mean.thick.ssl.pc2[7,2]*(final.data$log.ssl*final.data$Comp.2)
plot((final.data$log.ssl*final.data$Comp.2), mean.thick.ssl.pc2.output)