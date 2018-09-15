Code for manuscript titled: Alternative designs in tropical tree seedling performance landscapes

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

# Packages to run mixed-effects models, get significance of effects, get marginal and conditional R-square values
library(lme4)
library(lmerTest)
library(piecewiseSEM)

# Single-Trait Models; Table 1

lma.rgr=lmer(log.rgr~log.lma+(1|Plot)+(1|sp),all.data)
rmf.rgr=lmer(log.rgr~log.rmf+(1|Plot)+(1|sp),all.data)
lma.light=lmer(log.lma~log.light+(1|Plot)+(1|sp),all.data)
rmf.light=lmer(log.rmf~log.light+(1|Plot)+(1|sp),all.data)

sem.model.fits(list(lma.rgr, rmf.rgr, lma.light, rmf.light))

library(lattice)

dotplot(ranef(lma.rgr, condVar=TRUE), scales=list(cex=0.2))
dotplot(ranef(rmf.rgr, condVar=TRUE), scales=list(cex=0.2))
dotplot(ranef(lma.light, condVar=TRUE), scales=list(cex=0.2))
dotplot(ranef(rmf.light, condVar=TRUE), scales=list(cex=0.2))


# Two-way interaction models; Table 2

lma.light.model=lmer(log.rgr~log.lma+log.light+log.lma*log.light+(1|Plot)+(1|sp), all.data)
rmf.light.model=lmer(log.rgr~log.rmf+log.light+log.rmf*log.light+(1|Plot)+(1|sp), all.data)

sem.model.fits(list(rmf.light.model,lma.light.model))

dotplot(ranef(lma.light.model, condVar=TRUE),scales=list(cex=0.2))
dotplot(ranef(rmf.light.model, condVar=TRUE),scales=list(cex=0.2))


# Three-way interaction models; Table 3

lma.rmf.light=lmer(log.rgr~log.lma+log.rmf+log.light+log.lma*log.rmf*log.light+(1|Plot)+(1|sp), all.data)
lma.rmf.pc1=lmer(log.rgr~log.lma+log.rmf+Comp.1+log.lma*log.rmf*Comp.1+(1|Plot)+(1|sp), all.data)
lma.rmf.pc2=lmer(log.rgr~log.lma+log.rmf+Comp.2+log.lma*log.rmf*Comp.2+(1|Plot)+(1|sp), all.data)
lma.mean.thick.pc1=lmer(log.rgr~log.lma+log.mean.thick+Comp.1+log.lma*log.mean.thick*Comp.1+(1|Plot)+(1|sp), all.data)
thick.ssl.pc2=lmer(log.rgr~log.mean.thick+log.ssl+Comp.2+log.mean.thick*log.ssl*Comp.2+(1|Plot)+(1|sp), all.data)
thick.lar1.pc2=lmer(log.rgr~log.mean.thick+log.lar1+Comp.2+log.mean.thick*log.lar1*Comp.2+(1|Plot)+(1|sp), all.data)
thick.smf.pc1=lmer(log.rgr~log.mean.thick+log.smf+Comp.1+log.mean.thick*log.smf*Comp.1+(1|Plot)+(1|sp), all.data)
thick.rmf.pc1=lmer(log.rgr~log.mean.thick+log.rmf+Comp.1+log.mean.thick*log.rmf*Comp.1+(1|Plot)+(1|sp), all.data)
lar.lmf.pc2=lmer(log.rgr~log.lar1+log.lmf+Comp.2+log.lar1*log.lmf*Comp.2+(1|Plot)+(1|sp), all.data)
ssl.smf.light=lmer(log.rgr~log.ssl+log.smf+log.light+log.ssl*log.smf*log.light+(1|Plot)+(1|sp), all.data)

sem.model.fits(list(lma.rmf.light,lma.rmf.pc1,lma.rmf.pc2,lma.mean.thick.pc1,thick.ssl.pc2,thick.lar1.pc2,thick.smf.pc1,thick.rmf.pc1,lar.lmf.pc2,ssl.smf.light))

dotplot(ranef(lma.rmf.light, condVar=TRUE),scales=list(cex=0.2))
dotplot(ranef(lma.rmf.pc1, condVar=TRUE),scales=list(cex=0.2))
dotplot(ranef(lma.rmf.pc2, condVar=TRUE),scales=list(cex=0.2))
dotplot(ranef(lma.mean.thick.pc1, condVar=TRUE),scales=list(cex=0.2))
dotplot(ranef(thick.ssl.pc2, condVar=TRUE),scales=list(cex=0.2))
dotplot(ranef(thick.lar1.pc2, condVar=TRUE),scales=list(cex=0.2))
dotplot(ranef(thick.smf.pc1, condVar=TRUE),scales=list(cex=0.2))
dotplot(ranef(thick.rmf.pc1, condVar=TRUE),scales=list(cex=0.2))
dotplot(ranef(lar.lmf.pc2, condVar=TRUE),scales=list(cex=0.2))
dotplot(ranef(ssl.smf.light, condVar=TRUE),scales=list(cex=0.2))


# Four-way interaction models; Table S2

lma.rmf.ltpc1=lmer(log.rgr~log.lma*log.rmf*log.light*Comp.1+(1|Plot)+(1|sp), all.data)
lma.rmf.ltpc2=lmer(log.rgr~log.lma*log.rmf*log.light*Comp.2+(1|Plot)+(1|sp), all.data)
lma.rmf.ltpc3=lmer(log.rgr~log.lma*log.rmf*log.light*Comp.3+(1|Plot)+(1|sp), all.data)
lma.mean.thick.ltpc1=lmer(log.rgr~log.lma*log.mean.thick*log.light*Comp.1+(1|Plot)+(1|sp), all.data)
thick.ssl.ltpc2=lmer(log.rgr~log.mean.thick*log.ssl*log.light*Comp.2+(1|Plot)+(1|sp),all.data)
thick.lar1.ltpc2=lmer(log.rgr~log.mean.thick*log.lar1*log.light*Comp.2+(1|Plot)+(1|sp),all.data)
thick.smf.ltpc1=lmer(log.rgr~log.mean.thick*log.smf*log.light*Comp.1+(1|Plot)+(1|sp),all.data)
thick.rmf.ltpc1=lmer(log.rgr~log.mean.thick*log.rmf*log.light*Comp.1+(1|Plot)+(1|sp),all.data)
lar1.lmf.ltpc2=lmer(log.rgr~log.lar1*log.lmf*log.light*Comp.2+(1|Plot)+(1|sp),all.data)
ssl.smf.ltpc1=lmer(log.rgr~log.ssl*log.smf*log.light*Comp.1+(1|Plot)+(1|sp),all.data)
ssl.smf.ltpc2=lmer(log.rgr~log.ssl*log.smf*log.light*Comp.2+(1|Plot)+(1|sp),all.data)
ssl.smf.ltpc3=lmer(log.rgr~log.ssl*log.smf*log.light*Comp.3+(1|Plot)+(1|sp),all.data)


sem.model.fits(list(lma.rmf.ltpc1,lma.rmf.ltpc2,lma.rmf.ltpc3, lma.mean.thick.ltpc1, 
thick.ssl.ltpc2, thick.lar1.ltpc2, thick.smf.ltpc1, thick.rmf.ltpc1, lar1.lmf.ltpc2, ssl.smf.ltpc1, ssl.smf.ltpc2, ssl.smf.ltpc3))


dotplot(ranef(lma.rmf.ltpc1, condVar=TRUE),scales=list(cex=0.2))
dotplot(ranef(lma.rmf.ltpc2, condVar=TRUE),scales=list(cex=0.2))
dotplot(ranef(lma.rmf.ltpc3, condVar=TRUE),scales=list(cex=0.2))
dotplot(ranef(lma.mean.thick.ltpc1, condVar=TRUE),scales=list(cex=0.2))
dotplot(ranef(thick.ssl.ltpc2, condVar=TRUE),scales=list(cex=0.2))
dotplot(ranef(thick.lar1.ltpc2, condVar=TRUE),scales=list(cex=0.2))
dotplot(ranef(thick.smf.ltpc1, condVar=TRUE),scales=list(cex=0.2))
dotplot(ranef(thick.rmf.ltpc1, condVar=TRUE),scales=list(cex=0.2))
dotplot(ranef(lar1.lmf.ltpc2, condVar=TRUE),scales=list(cex=0.2))
dotplot(ranef(ssl.smf.ltpc1, condVar=TRUE),scales=list(cex=0.2))
dotplot(ranef(ssl.smf.ltpc2, condVar=TRUE),scales=list(cex=0.2))
dotplot(ranef(ssl.smf.ltpc3, condVar=TRUE),scales=list(cex=0.2))


# Slopes for two-way interaction models

lma.light.output=coef(summary(lma.light.model))[2,1]+coef(summary(lma.light.model))[4,1]*all.data$log.light
rmf.light.output=coef(summary(rmf.light.model))[2,1]+coef(summary(rmf.light.model))[4,1]*all.data$log.light

# Plotting slopes against environment

plot(all.data$log.light, lma.light.output)
plot(all.data$log.light, rmf.light.output)

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
visreg2d(rmf.light.model, "log.rmf", "log.light", plot.type="image", xlab="Root Mass Fraction)", ylab="Light", main="")
visreg2d(rmf.light.model, "log.rmf", "log.light", plot.type="persp", ylab="Light", zlab="\nRelative Growth Rate", xlab="Root Mass Fraction)", main="Relative Growth Rate = f(Root Mass Fraction*Light)",cex.main=1,nn=99, cex.lab=1,lwd=0.5, border="grey40", col=adjustcolor("blue",alpha.f=.5))

# Conceptual model, Figure 1

library(fields)

grid.l <- list(seq(-3,3,length=100), seq(-3,3,length=100))

te <- make.surface.grid(grid.l)

par(mfrow=c(1,2),mar=c(2,2,2,2))

te.linear.pred1 <- as.surface(te, 10*te[,1]*te[,2]*-2 + te[,1]*te[,2] + te[,2]*-2 + te[,1]*-2 + te[,1] + te[,2] + -2)

persp(te.linear.pred1$x,te.linear.pred1$y,te.linear.pred1$z, theta = 30, phi = 20,

      xlab="\nTrait X",ylab="\nTrait Y",zlab="\nPerformance",

      main="(A) Low end of environmental gradient",cex.lab=1,lwd=0.5, border="grey40", col=adjustcolor("blue",alpha.f=.5))


te.linear.pred2 <- as.surface(te, 10*te[,1]*te[,2]*2 + te[,1]*te[,2] + te[,2]*2 + te[,1]*2 + te[,1] + te[,2] + 2)

persp(te.linear.pred2$x,te.linear.pred2$y,te.linear.pred2$z, theta = 30, phi = 20,

      xlab="\nTrait X",ylab="\nTrait Y",zlab="\nPerformance",

      main="(B) High end of environmental gradient",cex.lab=1,lwd=0.5, border="grey40", col=adjustcolor("blue",alpha.f=.5))

# Level plots; Figure S1-S7

# LMA x RMF x Light; Figure 3 and Figure S1

# Vector of the range of light by 0.1

vec=seq(-2.6,3.0, by=0.1)

# 10% quantiles for the range of light

quantile(vec, probs=seq(0,1,0.1))

# Expand the range of the variables in the model while holding light constant at each 10% of its range

data.low=with(all.data,expand.grid(log.lma=seq(-3.5,3.2,length.out=10),log.rmf=seq(-5.5,3.4,length.out=10),log.light=-2.6))
data.10=with(all.data,expand.grid(log.lma=seq(-3.5,3.2,length.out=10),log.rmf=seq(-5.5,3.4,length.out=10),log.light=-2.04))
data.20=with(all.data,expand.grid(log.lma=seq(-3.5,3.2,length.out=10),log.rmf=seq(-5.5,3.4,length.out=10),log.light=-1.48))
data.30=with(all.data,expand.grid(log.lma=seq(-3.5,3.2,length.out=10),log.rmf=seq(-5.5,3.4,length.out=10),log.light=-0.92))
data.40=with(all.data,expand.grid(log.lma=seq(-3.5,3.2,length.out=10),log.rmf=seq(-5.5,3.4,length.out=10),log.light=-0.36))
data.med=with(all.data,expand.grid(log.lma=seq(-3.5,3.2,length.out=10),log.rmf=seq(-5.5,3.4,length.out=10),log.light=0.20))
data.60=with(all.data,expand.grid(log.lma=seq(-3.5,3.2,length.out=10),log.rmf=seq(-5.5,3.4,length.out=10),log.light=0.76))
data.70=with(all.data,expand.grid(log.lma=seq(-3.5,3.2,length.out=10),log.rmf=seq(-5.5,3.4,length.out=10),log.light=1.32))
data.80=with(all.data,expand.grid(log.lma=seq(-3.5,3.2,length.out=10),log.rmf=seq(-5.5,3.4,length.out=10),log.light=1.88))
data.90=with(all.data,expand.grid(log.lma=seq(-3.5,3.2,length.out=10),log.rmf=seq(-5.5,3.4,length.out=10),log.light=2.44))
data.high=with(all.data,expand.grid(log.lma=seq(-3.5,3.2,length.out=10),log.rmf=seq(-5.5,3.4,length.out=10),log.light=3.0))

# Predict new values from the fitted model with the expanded data

p.low=predict(lma.rmf.light,data.low, re.form=NA)
p.10=predict(lma.rmf.light,data.10, re.form=NA)
p.20=predict(lma.rmf.light,data.20, re.form=NA)
p.30=predict(lma.rmf.light,data.30, re.form=NA)
p.40=predict(lma.rmf.light,data.40, re.form=NA)
p.med=predict(lma.rmf.light,data.med, re.form=NA)
p.60=predict(lma.rmf.light,data.60, re.form=NA)
p.70=predict(lma.rmf.light,data.70, re.form=NA)
p.80=predict(lma.rmf.light,data.80, re.form=NA)
p.90=predict(lma.rmf.light,data.90, re.form=NA)
p.high=predict(lma.rmf.light,data.high, re.form=NA)

# Put the two trait values and the predicted values into a data.frame

light.low.df=data.low[,1:2]
light.low.df[,3]=p.low
colnames(light.low.df)[3]="rgr"
light.10.df=data.10[,1:2]
light.10.df[,3]=p.10
light.20.df=data.20[,1:2]
light.20.df[,3]=p.20
light.30.df=data.30[,1:2]
light.30.df[,3]=p.30
light.40.df=data.40[,1:2]
light.40.df[,3]=p.40
light.med.df=data.med[,1:2]
light.med.df[,3]=p.med
colnames(light.med.df)[3]="rgr"
light.60.df=data.60[,1:2]
light.60.df[,3]=p.60
light.70.df=data.70[,1:2]
light.70.df[,3]=p.70
light.80.df=data.80[,1:2]
light.80.df[,3]=p.80
light.90.df=data.90[,1:2]
light.90.df[,3]=p.90
light.high.df=data.high[,1:2]
light.high.df[,3]=p.high
colnames(light.high.df)[3]="rgr"

# Fit linear models of predicted data for graphical formatting

new.model.low=lm(rgr~log.lma*log.rmf, data=light.low.df)
new.model.med=lm(rgr~log.lma*log.rmf, data=light.med.df)
new.model.high=lm(rgr~log.lma*log.rmf, data=light.high.df)

# Figure 3

par(mfrow=c(1,3))
visreg2d(new.model.low, "log.lma", "log.rmf", plot.type="persp", ylab="\nRoot Mass Fraction", zlab="\nRelative Growth Rate", xlab="\nLeaf Mass per Area", cex.main=1, theta=30,cex.lab=1,lwd=0.5, border="grey40", col=adjustcolor("blue",alpha.f=.5))
visreg2d(new.model.med, "log.lma", "log.rmf", plot.type="persp", ylab="\nRoot Mass Fraction", zlab="\nRelative Growth Rate", xlab="\nLeaf Mass per Area", cex.main=1, theta=30, cex.lab=1,lwd=0.5, border="grey40", col=adjustcolor("blue",alpha.f=.5))
visreg2d(new.model.high, "log.lma", "log.rmf", plot.type="persp", ylab="\nRoot Mass Fraction", zlab="\nRelative Growth Rate", xlab="\nLeaf Mass per Area",cex.main=1, theta=30, cex.lab=1,lwd=0.5, border="grey40", col=adjustcolor("blue",alpha.f=.5))

# Figure S1

library(lattice)
library(gridExtra)
library(colorspace)

light.low=levelplot(rgr~log.lma*log.rmf, data=light.low.df, col.regions=diverge_hcl(50),at=seq(-2.5,6.5, 1), main=list(label="0%", cex=.75), ylab=list(label="Root Mass Fraction", cex=.75),xlab=list(label="Leaf Mass per Area", cex=.75))
light.10=levelplot(V3~log.lma*log.rmf, data=light.10.df, col.regions=diverge_hcl(50),at=seq(-2.5,6.5, 1), main=list(label="10%", cex=.75), ylab="",xlab=list(label="Leaf Mass per Area", cex=.75))
light.20=levelplot(V3~log.lma*log.rmf, data=light.20.df,col.regions=diverge_hcl(50),at=seq(-2.5,6.5, 1), main=list(label="20%", cex=.75), ylab="",xlab=list(label="Leaf Mass per Area", cex=.75))
light.30=levelplot(V3~log.lma*log.rmf, data=light.30.df, col.regions=diverge_hcl(50),at=seq(-2.5,6.5, 1), main=list(label="30%", cex=.75),ylab="",xlab=list(label="Leaf Mass per Area", cex=.75))
light.40=levelplot(V3~log.lma*log.rmf, data=light.40.df, col.regions=diverge_hcl(50),at=seq(-2.5,6.5, 1), main=list(label="40%", cex=.75),ylab=list(label="Root Mass Fraction", cex=.75),xlab=list(label="Leaf Mass per Area", cex=.75))
light.med=levelplot(rgr~log.lma*log.rmf, data=light.med.df, col.regions=diverge_hcl(50),at=seq(-2.5,6.5, 1), main=list(label="50%", cex=.75),ylab="",xlab=list(label="Leaf Mass per Area", cex=.75))
light.60=levelplot(V3~log.lma*log.rmf, data=light.60.df, col.regions=diverge_hcl(50),at=seq(-2.5,6.5, 1), main=list(label="60%", cex=.75),ylab="",xlab=list(label="Leaf Mass per Area", cex=.75))
light.70=levelplot(V3~log.lma*log.rmf, data=light.70.df, col.regions=diverge_hcl(50),at=seq(-2.5,6.5, 1), main=list(label="70%", cex=.75),ylab="",xlab=list(label="Leaf Mass per Area", cex=.75))
light.80=levelplot(V3~log.lma*log.rmf, data=light.80.df, col.regions=diverge_hcl(50),at=seq(-2.5,6.5, 1), main=list(label="80%", cex=.75),ylab=list(label="Root Mass Fraction",cex=.75), xlab=list(label="Leaf Mass per Area", cex=.75))
light.90=levelplot(V3~log.lma*log.rmf, data=light.90.df, col.regions=diverge_hcl(50),at=seq(-2.5,6.5, 1), main=list(label="90%", cex=.75), xlab=list(label="Leaf Mass per Area", cex=.75),ylab="")
light.high=levelplot(rgr~log.lma*log.rmf, data=light.high.df, col.regions=diverge_hcl(50),at=seq(-2.5,6.5, 1), main=list(label="100%", cex=.75), xlab=list(label="Leaf Mass per Area", cex=.75),ylab="")
grid.arrange(light.low,light.10,light.20,light.30,light.40,light.med,light.60,light.70,light.80,light.90,light.high, nrow=3,ncol=4)

# Figure S2

#  mean leaf thickness x SMF x soil PC1

# Vector of the range of soil PC1 by 0.3

vec=seq(-6.0,6.1, by=0.3)

# 10% quantiles for the range of soil PC1

quantile(vec, probs=seq(0,1,0.1))

# Expand the range of the variables in the model while holding soil PC1 constant at each 10% of its range

data.low=with(all.data,expand.grid(log.mean.thick=seq(-2.1,3.2,length.out=10),log.smf=seq(-6.5,3.0,length.out=10), Comp.1=-6.0))
data.10=with(all.data,expand.grid(log.mean.thick=seq(-2.1,3.2,length.out=10),log.smf=seq(-6.5,3.0,length.out=10), Comp.1=-4.8))
data.20=with(all.data,expand.grid(log.mean.thick=seq(-2.1,3.2,length.out=10),log.smf=seq(-6.5,3.0,length.out=10), Comp.1=-3.6))
data.30=with(all.data,expand.grid(log.mean.thick=seq(-2.1,3.2,length.out=10),log.smf=seq(-6.5,3.0,length.out=10), Comp.1=-2.4))
data.40=with(all.data,expand.grid(log.mean.thick=seq(-2.1,3.2,length.out=10),log.smf=seq(-6.5,3.0,length.out=10), Comp.1=-1.2))
data.med=with(all.data,expand.grid(log.mean.thick=seq(-2.1,3.2,length.out=10),log.smf=seq(-6.5,3.0,length.out=10), Comp.1=0.0))
data.60=with(all.data,expand.grid(log.mean.thick=seq(-2.1,3.2,length.out=10),log.smf=seq(-6.5,3.0,length.out=10), Comp.1=1.2))
data.70=with(all.data,expand.grid(log.mean.thick=seq(-2.1,3.2,length.out=10),log.smf=seq(-6.5,3.0,length.out=10), Comp.1=2.4))
data.80=with(all.data,expand.grid(log.mean.thick=seq(-2.1,3.2,length.out=10),log.smf=seq(-6.5,3.0,length.out=10), Comp.1=3.6))
data.90=with(all.data,expand.grid(log.mean.thick=seq(-2.1,3.2,length.out=10),log.smf=seq(-6.5,3.0,length.out=10), Comp.1=4.8))
data.high=with(all.data,expand.grid(log.mean.thick=seq(-2.1,3.2,length.out=10),log.smf=seq(-6.5,3.0,length.out=10), Comp.1=6.1))

# Predict new values from the fitted model with the expanded data

p.low=predict(thick.smf.pc1,data.low, re.form=NA)
p.10=predict(thick.smf.pc1,data.10, re.form=NA)
p.20=predict(thick.smf.pc1,data.20, re.form=NA)
p.30=predict(thick.smf.pc1,data.30, re.form=NA)
p.40=predict(thick.smf.pc1,data.40, re.form=NA)
p.med=predict(thick.smf.pc1,data.med, re.form=NA)
p.60=predict(thick.smf.pc1,data.60, re.form=NA)
p.70=predict(thick.smf.pc1,data.70, re.form=NA)
p.80=predict(thick.smf.pc1,data.80, re.form=NA)
p.90=predict(thick.smf.pc1,data.90, re.form=NA)
p.high=predict(thick.smf.pc1,data.high, re.form=NA)

# Put the two trait values and the predicted values into a data.frame

comp1.low.df=data.low[,1:2]
comp1.low.df[,3]=p.low
comp1.10.df=data.10[,1:2]
comp1.10.df[,3]=p.10
comp1.20.df=data.20[,1:2]
comp1.20.df[,3]=p.20
comp1.30.df=data.30[,1:2]
comp1.30.df[,3]=p.30
comp1.40.df=data.40[,1:2]
comp1.40.df[,3]=p.40
comp1.med.df=data.med[,1:2]
comp1.med.df[,3]=p.med
comp1.60.df=data.60[,1:2]
comp1.60.df[,3]=p.60
comp1.70.df=data.70[,1:2]
comp1.70.df[,3]=p.70
comp1.80.df=data.80[,1:2]
comp1.80.df[,3]=p.80
comp1.90.df=data.90[,1:2]
comp1.90.df[,3]=p.90
comp1.high.df=data.high[,1:2]
comp1.high.df[,3]=p.high

# Figure S2

comp1.low=levelplot(V3~log.mean.thick*log.smf, data=comp1.low.df, col.regions=diverge_hcl(50),at=seq(-2.5,5.5, 1), main=list(label="0%", cex=.75), ylab=list(label="Stem Mass Fraction", cex=.75), xlab=list(label="Mean Leaf Thickness", cex=.75))
comp1.10=levelplot(V3~log.mean.thick*log.smf, data=comp1.10.df, col.regions=diverge_hcl(50),at=seq(-2.5,5.5, 1),main=list(label="10%", cex=.75),xlab=list(label="Mean Leaf Thickness", cex=.75), ylab="")
comp1.20=levelplot(V3~log.mean.thick*log.smf, data=comp1.20.df,col.regions=diverge_hcl(50),at=seq(-2.5,5.5, 1), main=list(label="20%", cex=.75),xlab=list(label="Mean Leaf Thickness", cex=.75),ylab="")
comp1.30=levelplot(V3~log.mean.thick*log.smf, data=comp1.30.df, col.regions=diverge_hcl(50),at=seq(-2.5,5.5, 1), main=list(label="30%", cex=.75),xlab=list(label="Mean Leaf Thickness", cex=.75),ylab="")
comp1.40=levelplot(V3~log.mean.thick*log.smf, data=comp1.40.df, col.regions=diverge_hcl(50),at=seq(-2.5,5.5, 1), main=list(label="40%", cex=.75),ylab=list(label="Stem Mass Fraction", cex=.75), xlab=list(label="Mean Leaf Thickness", cex=.75))
comp1.med=levelplot(V3~log.mean.thick*log.smf, data=comp1.med.df, col.regions=diverge_hcl(50),at=seq(-2.5,5.5, 1),main=list(label="50%", cex=.75),xlab=list(label="Mean Leaf Thickness", cex=.75),ylab="")
comp1.60=levelplot(V3~log.mean.thick*log.smf, data=comp1.60.df, col.regions=diverge_hcl(50),at=seq(-2.5,5.5, 1), main=list(label="60%", cex=.75),xlab=list(label="Mean Leaf Thickness", cex=.75),ylab="")
comp1.70=levelplot(V3~log.mean.thick*log.smf, data=comp1.70.df, col.regions=diverge_hcl(50),at=seq(-2.5,5.5, 1), main=list(label="70%", cex=.75),xlab=list(label="Mean Leaf Thickness", cex=.75),ylab="")
comp1.80=levelplot(V3~log.mean.thick*log.smf, data=comp1.80.df, col.regions=diverge_hcl(50),at=seq(-2.5,5.5, 1), main=list(label="80%", cex=.75),ylab=list(label="Stem Mass Fraction", cex=.75), xlab=list(label="Mean Leaf Thickness", cex=.75))
comp1.90=levelplot(V3~log.mean.thick*log.smf, data=comp1.90.df, col.regions=diverge_hcl(50),at=seq(-2.5,5.5, 1), main=list(label="90%", cex=.75),xlab=list(label="Mean Leaf Thickness", cex=.75), ylab="")
comp1.high=levelplot(V3~log.mean.thick*log.smf, data=comp1.high.df, col.regions=diverge_hcl(50),at=seq(-2.5,5.5, 1), main=list(label="100%", cex=.75),xlab=list(label="Mean Leaf Thickness", cex=.75), ylab="")
grid.arrange(comp1.low,comp1.10,comp1.20,comp1.30,comp1.40,comp1.med,comp1.60,comp1.70,comp1.80,comp1.90,comp1.high, nrow=3,ncol=4)


# Figure S3

#  mean leaf thickness x RMF x soil PC1

# Vector of the range of soil PC1 by 0.3

vec=seq(-6.0,6.1, by=0.3)

# 10% quantiles for the range of soil PC1

quantile(vec, probs=seq(0,1,0.1))

# Expand the range of the variables in the model while holding soil PC1 constant at each 10% of its range

data.low=with(all.data,expand.grid(log.mean.thick=seq(-2.1,3.2,length.out=10),log.rmf=seq(-5.5,3.4,length.out=10), Comp.1=-6.0))
data.10=with(all.data,expand.grid(log.mean.thick=seq(-2.1,3.2,length.out=10),log.rmf=seq(-5.5,3.4,length.out=10), Comp.1=-4.8))
data.20=with(all.data,expand.grid(log.mean.thick=seq(-2.1,3.2,length.out=10),log.rmf=seq(-5.5,3.4,length.out=10), Comp.1=-3.6))
data.30=with(all.data,expand.grid(log.mean.thick=seq(-2.1,3.2,length.out=10),log.rmf=seq(-5.5,3.4,length.out=10), Comp.1=-2.4))
data.40=with(all.data,expand.grid(log.mean.thick=seq(-2.1,3.2,length.out=10),log.rmf=seq(-5.5,3.4,length.out=10), Comp.1=-1.2))
data.med=with(all.data,expand.grid(log.mean.thick=seq(-2.1,3.2,length.out=10),log.rmf=seq(-5.5,3.4,length.out=10), Comp.1=0.0))
data.60=with(all.data,expand.grid(log.mean.thick=seq(-2.1,3.2,length.out=10),log.rmf=seq(-5.5,3.4,length.out=10), Comp.1=1.2))
data.70=with(all.data,expand.grid(log.mean.thick=seq(-2.1,3.2,length.out=10),log.rmf=seq(-5.5,3.4,length.out=10), Comp.1=2.4))
data.80=with(all.data,expand.grid(log.mean.thick=seq(-2.1,3.2,length.out=10),log.rmf=seq(-5.5,3.4,length.out=10), Comp.1=3.6))
data.90=with(all.data,expand.grid(log.mean.thick=seq(-2.1,3.2,length.out=10),log.rmf=seq(-5.5,3.4,length.out=10), Comp.1=4.8))
data.high=with(all.data,expand.grid(log.mean.thick=seq(-2.1,3.2,length.out=10),log.rmf=seq(-5.5,3.4,length.out=10), Comp.1=6.1))

# Predict new values from the fitted model with the expanded data

p.low=predict(thick.rmf.pc1,data.low, re.form=NA)
p.10=predict(thick.rmf.pc1,data.10, re.form=NA)
p.20=predict(thick.rmf.pc1,data.20, re.form=NA)
p.30=predict(thick.rmf.pc1,data.30, re.form=NA)
p.40=predict(thick.rmf.pc1,data.40, re.form=NA)
p.med=predict(thick.rmf.pc1,data.med, re.form=NA)
p.60=predict(thick.rmf.pc1,data.60, re.form=NA)
p.70=predict(thick.rmf.pc1,data.70, re.form=NA)
p.80=predict(thick.rmf.pc1,data.80, re.form=NA)
p.90=predict(thick.rmf.pc1,data.90, re.form=NA)
p.high=predict(thick.rmf.pc1,data.high, re.form=NA)

# Put the two trait values and the predicted values into a data.frame

comp1.low.df=data.low[,1:2]
comp1.low.df[,3]=p.low
comp1.10.df=data.10[,1:2]
comp1.10.df[,3]=p.10
comp1.20.df=data.20[,1:2]
comp1.20.df[,3]=p.20
comp1.30.df=data.30[,1:2]
comp1.30.df[,3]=p.30
comp1.40.df=data.40[,1:2]
comp1.40.df[,3]=p.40
comp1.med.df=data.med[,1:2]
comp1.med.df[,3]=p.med
comp1.60.df=data.60[,1:2]
comp1.60.df[,3]=p.60
comp1.70.df=data.70[,1:2]
comp1.70.df[,3]=p.70
comp1.80.df=data.80[,1:2]
comp1.80.df[,3]=p.80
comp1.90.df=data.90[,1:2]
comp1.90.df[,3]=p.90
comp1.high.df=data.high[,1:2]
comp1.high.df[,3]=p.high

# Figure S3

comp1.low=levelplot(V3~log.mean.thick*log.rmf, data=comp1.low.df, col.regions=diverge_hcl(50),at=seq(-4.5,4.5, 1), main=list(label="0%", cex=.75), ylab=list(label="Root Mass Fraction", cex=.75), xlab=list(label="Mean Leaf Thickness", cex=.75))
comp1.10=levelplot(V3~log.mean.thick*log.rmf, data=comp1.10.df, col.regions=diverge_hcl(50),at=seq(-4.5,4.5, 1),main=list(label="10%", cex=.75),xlab=list(label="Mean Leaf Thickness", cex=.75), ylab="")
comp1.20=levelplot(V3~log.mean.thick*log.rmf, data=comp1.20.df,col.regions=diverge_hcl(50),at=seq(-4.5,4.5, 1), main=list(label="20%", cex=.75),xlab=list(label="Mean Leaf Thickness", cex=.75),ylab="")
comp1.30=levelplot(V3~log.mean.thick*log.rmf, data=comp1.30.df, col.regions=diverge_hcl(50),at=seq(-4.5,4.5, 1), main=list(label="30%", cex=.75),xlab=list(label="Mean Leaf Thickness", cex=.75),ylab="")
comp1.40=levelplot(V3~log.mean.thick*log.rmf, data=comp1.40.df, col.regions=diverge_hcl(50),at=seq(-4.5,4.5, 1), main=list(label="40%", cex=.75),ylab=list(label="Root Mass Fraction", cex=.75), xlab=list(label="Mean Leaf Thickness", cex=.75))
comp1.med=levelplot(V3~log.mean.thick*log.rmf, data=comp1.med.df, col.regions=diverge_hcl(50),at=seq(-4.5,4.5, 1),main=list(label="50%", cex=.75),xlab=list(label="Mean Leaf Thickness", cex=.75),ylab="")
comp1.60=levelplot(V3~log.mean.thick*log.rmf, data=comp1.60.df, col.regions=diverge_hcl(50),at=seq(-4.5,4.5, 1), main=list(label="60%", cex=.75),xlab=list(label="Mean Leaf Thickness", cex=.75),ylab="")
comp1.70=levelplot(V3~log.mean.thick*log.rmf, data=comp1.70.df, col.regions=diverge_hcl(50),at=seq(-4.5,4.5, 1), main=list(label="70%", cex=.75),xlab=list(label="Mean Leaf Thickness", cex=.75),ylab="")
comp1.80=levelplot(V3~log.mean.thick*log.rmf, data=comp1.80.df, col.regions=diverge_hcl(50),at=seq(-4.5,4.5, 1), main=list(label="80%", cex=.75),ylab=list(label="Root Mass Fraction", cex=.75), xlab=list(label="Mean Leaf Thickness", cex=.75))
comp1.90=levelplot(V3~log.mean.thick*log.rmf, data=comp1.90.df, col.regions=diverge_hcl(50),at=seq(-4.5,4.5, 1), main=list(label="90%", cex=.75),xlab=list(label="Mean Leaf Thickness", cex=.75), ylab="")
comp1.high=levelplot(V3~log.mean.thick*log.rmf, data=comp1.high.df, col.regions=diverge_hcl(50),at=seq(-4.5,4.5, 1), main=list(label="100%", cex=.75),xlab=list(label="Mean Leaf Thickness", cex=.75), ylab="")
grid.arrange(comp1.low,comp1.10,comp1.20,comp1.30,comp1.40,comp1.med,comp1.60,comp1.70,comp1.80,comp1.90,comp1.high, nrow=3,ncol=4)

# Figure S4

# LMA x RMF x Soil PC2

# Vector of the range of soil PC2 by 0.1

vec=seq(-5.4,4.0, by=0.1)

# 10% quantiles for the range of soil PC2

quantile(vec, probs=seq(0,1,0.1))

# Expand the range of the variables in the model while holding soil PC2 constant at each 10% of its range

comp.2.data.low=with(all.data,expand.grid(log.lma=seq(-3.5,3.2,length.out=10),log.rmf=seq(-5.5,3.4,length.out=10), Comp.2=-5.40))
comp.2.data.10=with(all.data,expand.grid(log.lma=seq(-3.5,3.2,length.out=10),log.rmf=seq(-5.5,3.4,length.out=10), Comp.2=-4.46))
comp.2.data.20=with(all.data,expand.grid(log.lma=seq(-3.5,3.2,length.out=10),log.rmf=seq(-5.5,3.4,length.out=10), Comp.2=-3.52))
comp.2.data.30=with(all.data,expand.grid(log.lma=seq(-3.5,3.2,length.out=10),log.rmf=seq(-5.5,3.4,length.out=10), Comp.2=-2.58))
comp.2.data.40=with(all.data,expand.grid(log.lma=seq(-3.5,3.2,length.out=10),log.rmf=seq(-5.5,3.4,length.out=10), Comp.2=-1.64))
comp.2.data.med=with(all.data,expand.grid(log.lma=seq(-3.5,3.2,length.out=10),log.rmf=seq(-5.5,3.4,length.out=10), Comp.2=-0.70))
comp.2.data.60=with(all.data,expand.grid(log.lma=seq(-3.5,3.2,length.out=10),log.rmf=seq(-5.5,3.4,length.out=10), Comp.2=0.24))
comp.2.data.70=with(all.data,expand.grid(log.lma=seq(-3.5,3.2,length.out=10),log.rmf=seq(-5.5,3.4,length.out=10), Comp.2=1.18))
comp.2.data.80=with(all.data,expand.grid(log.lma=seq(-3.5,3.2,length.out=10),log.rmf=seq(-5.5,3.4,length.out=10), Comp.2=2.12))
comp.2.data.90=with(all.data,expand.grid(log.lma=seq(-3.5,3.2,length.out=10),log.rmf=seq(-5.5,3.4,length.out=10), Comp.2=3.06))
comp.2.data.high=with(all.data,expand.grid(log.lma=seq(-3.5,3.2,length.out=10),log.rmf=seq(-5.5,3.4,length.out=10), Comp.2=4.0))

# Predict new values from the fitted model with the expanded data

comp.2.p.low=predict(lma.rmf.pc2,comp.2.data.low, re.form=NA)
comp.2.p.10=predict(lma.rmf.pc2,comp.2.data.10, re.form=NA)
comp.2.p.20=predict(lma.rmf.pc2,comp.2.data.20, re.form=NA)
comp.2.p.30=predict(lma.rmf.pc2,comp.2.data.30, re.form=NA)
comp.2.p.40=predict(lma.rmf.pc2,comp.2.data.40, re.form=NA)
comp.2.p.med=predict(lma.rmf.pc2,comp.2.data.med, re.form=NA)
comp.2.p.60=predict(lma.rmf.pc2,comp.2.data.60, re.form=NA)
comp.2.p.70=predict(lma.rmf.pc2,comp.2.data.70, re.form=NA)
comp.2.p.80=predict(lma.rmf.pc2,comp.2.data.80, re.form=NA)
comp.2.p.90=predict(lma.rmf.pc2,comp.2.data.90, re.form=NA)
comp.2.p.high=predict(lma.rmf.pc2,comp.2.data.high, re.form=NA)

# Put the two trait values and the predicted values into a data.frame

comp2.low.df=comp.2.data.low[,1:2]
comp2.low.df[,3]=comp.2.p.low
comp2.10.df=comp.2.data.10[,1:2]
comp2.10.df[,3]=comp.2.p.10
comp2.20.df=comp.2.data.20[,1:2]
comp2.20.df[,3]=comp.2.p.20
comp2.30.df=comp.2.data.30[,1:2]
comp2.30.df[,3]=comp.2.p.30
comp2.40.df=comp.2.data.40[,1:2]
comp2.40.df[,3]=comp.2.p.40
comp2.med.df=comp.2.data.med[,1:2]
comp2.med.df[,3]=comp.2.p.med
comp2.60.df=comp.2.data.60[,1:2]
comp2.60.df[,3]=comp.2.p.60
comp2.70.df=comp.2.data.70[,1:2]
comp2.70.df[,3]=comp.2.p.70
comp2.80.df=comp.2.data.80[,1:2]
comp2.80.df[,3]=comp.2.p.80
comp2.90.df=comp.2.data.90[,1:2]
comp2.90.df[,3]=comp.2.p.90
comp2.high.df=comp.2.data.high[,1:2]
comp2.high.df[,3]=comp.2.p.high

# Figure S4

comp2.low=levelplot(V3~log.lma*log.rmf, data=comp2.low.df, col.regions=diverge_hcl(50),at=seq(-3.0,5.0, 1), main=list(label="0%", cex=.75), ylab=list(label="Root Mass Fraction", cex=.75), xlab=list(label="Leaf Mass per Area", cex=.75))
comp2.10=levelplot(V3~log.lma*log.rmf, data=comp2.10.df, col.regions=diverge_hcl(50),at=seq(-3.0,5.0, 1), main=list(label="10%", cex=.75),xlab=list(label="Leaf Mass per Area", cex=.75),ylab="")
comp2.20=levelplot(V3~log.lma*log.rmf, data=comp2.20.df,col.regions=diverge_hcl(50),at=seq(-3.0,5.0, 1), main=list(label="20%", cex=.75),xlab=list(label="Leaf Mass per Area", cex=.75),ylab="")
comp2.30=levelplot(V3~log.lma*log.rmf, data=comp2.30.df, col.regions=diverge_hcl(50),at=seq(-3.0,5.0, 1), main=list(label="30%", cex=.75),xlab=list(label="Leaf Mass per Area", cex=.75),ylab="")
comp2.40=levelplot(V3~log.lma*log.rmf, data=comp2.40.df, col.regions=diverge_hcl(50),at=seq(-3.0,5.0, 1), main=list(label="40%", cex=.75),ylab=list(label="Root Mass Fraction", cex=.75), xlab=list(label="Leaf Mass per Area", cex=.75))
comp2.med=levelplot(V3~log.lma*log.rmf, data=comp2.med.df, col.regions=diverge_hcl(50),at=seq(-3.0,5.0, 1), main=list(label="50%", cex=.75),xlab=list(label="Leaf Mass per Area", cex=.75),ylab="")
comp2.60=levelplot(V3~log.lma*log.rmf, data=comp2.60.df, col.regions=diverge_hcl(50),at=seq(-3.0,5.0, 1), main=list(label="60%", cex=.75),xlab=list(label="Leaf Mass per Area", cex=.75),ylab="")
comp2.70=levelplot(V3~log.lma*log.rmf, data=comp2.70.df, col.regions=diverge_hcl(50),at=seq(-3.0,5.0, 1), main=list(label="70%", cex=.75),xlab=list(label="Leaf Mass per Area", cex=.75),ylab="")
comp2.80=levelplot(V3~log.lma*log.rmf, data=comp2.80.df, col.regions=diverge_hcl(50),at=seq(-3.0,5.0, 1),  main=list(label="80%", cex=.75),ylab=list(label="Root Mass Fraction", cex=.75), xlab=list(label="Leaf Mass per Area", cex=.75))
comp2.90=levelplot(V3~log.lma*log.rmf, data=comp2.90.df, col.regions=diverge_hcl(50),at=seq(-3.0,5.0, 1), main=list(label="90%", cex=.75),xlab=list(label="Leaf Mass per Area", cex=.75), ylab="")
comp2.high=levelplot(V3~log.lma*log.rmf, data=comp2.high.df, col.regions=diverge_hcl(50),at=seq(-3.0,5.0, 1), main=list(label="100%", cex=.75),xlab=list(label="Leaf Mass per Area", cex=.75), ylab="")
grid.arrange(comp2.low,comp2.10,comp2.20,comp2.30,comp2.40,comp2.med,comp2.60,comp2.70,comp2.80,comp2.90,comp2.high, nrow=3,ncol=4)

# Figure S5

# mean leaf thickness x SSL x soil PC2

# Vector of the range of soil PC2 by 0.1

vec=seq(-5.4,4.0, by=0.1)

# 10% quantiles for the range of soil PC2

quantile(vec, probs=seq(0,1,0.1))

#Expand the range of the variables in the model while holding soil PC2 constant at each 10% of its range

data.low=with(all.data,expand.grid(log.mean.thick=seq(-2.1,3.2,length.out=10),log.ssl=seq(-3.3,2.9,length.out=10), Comp.2=-5.40))
data.10=with(all.data,expand.grid(log.mean.thick=seq(-2.1,3.2,length.out=10),log.ssl=seq(-3.3,2.9,length.out=10), Comp.2=-4.46))
data.20=with(all.data,expand.grid(log.mean.thick=seq(-2.1,3.2,length.out=10),log.ssl=seq(-3.3,2.9,length.out=10), Comp.2=-3.52))
data.30=with(all.data,expand.grid(log.mean.thick=seq(-2.1,3.2,length.out=10),log.ssl=seq(-3.3,2.9,length.out=10), Comp.2=-2.58))
data.40=with(all.data,expand.grid(log.mean.thick=seq(-2.1,3.2,length.out=10),log.ssl=seq(-3.3,2.9,length.out=10), Comp.2=-1.64))
data.med=with(all.data,expand.grid(log.mean.thick=seq(-2.1,3.2,length.out=10),log.ssl=seq(-3.3,2.9,length.out=10), Comp.2=-0.70))
data.60=with(all.data,expand.grid(log.mean.thick=seq(-2.1,3.2,length.out=10),log.ssl=seq(-3.3,2.9,length.out=10), Comp.2=0.24))
data.70=with(all.data,expand.grid(log.mean.thick=seq(-2.1,3.2,length.out=10),log.ssl=seq(-3.3,2.9,length.out=10), Comp.2=1.18))
data.80=with(all.data,expand.grid(log.mean.thick=seq(-2.1,3.2,length.out=10),log.ssl=seq(-3.3,2.9,length.out=10), Comp.2=2.12))
data.90=with(all.data,expand.grid(log.mean.thick=seq(-2.1,3.2,length.out=10),log.ssl=seq(-3.3,2.9,length.out=10), Comp.2=3.06))
data.high=with(all.data,expand.grid(log.mean.thick=seq(-2.1,3.2,length.out=10),log.ssl=seq(-3.3,2.9,length.out=10), Comp.2=4.0))

# Predict new values from the fitted model with the expanded data

p.low=predict(thick.ssl.pc2,data.low, re.form=NA)
p.10=predict(thick.ssl.pc2,data.10, re.form=NA)
p.20=predict(thick.ssl.pc2,data.20, re.form=NA)
p.30=predict(thick.ssl.pc2,data.30, re.form=NA)
p.40=predict(thick.ssl.pc2,data.40, re.form=NA)
p.med=predict(thick.ssl.pc2,data.med, re.form=NA)
p.60=predict(thick.ssl.pc2,data.60, re.form=NA)
p.70=predict(thick.ssl.pc2,data.70, re.form=NA)
p.80=predict(thick.ssl.pc2,data.80, re.form=NA)
p.90=predict(thick.ssl.pc2,data.90, re.form=NA)
p.high=predict(thick.ssl.pc2,data.high, re.form=NA)

# Put the two trait values and the predicted values into a data.frame

comp2.low.df=data.low[,1:2]
comp2.low.df[,3]=p.low
comp2.10.df=data.10[,1:2]
comp2.10.df[,3]=p.10
comp2.20.df=data.20[,1:2]
comp2.20.df[,3]=p.20
comp2.30.df=data.30[,1:2]
comp2.30.df[,3]=p.30
comp2.40.df=data.40[,1:2]
comp2.40.df[,3]=p.40
comp2.med.df=data.med[,1:2]
comp2.med.df[,3]=p.med
comp2.60.df=data.60[,1:2]
comp2.60.df[,3]=p.60
comp2.70.df=data.70[,1:2]
comp2.70.df[,3]=p.70
comp2.80.df=data.80[,1:2]
comp2.80.df[,3]=p.80
comp2.90.df=data.90[,1:2]
comp2.90.df[,3]=p.90
comp2.high.df=data.high[,1:2]
comp2.high.df[,3]=p.high

# Figure S5

comp2.low=levelplot(V3~log.mean.thick*log.ssl, data=comp2.low.df, col.regions=diverge_hcl(50),at=seq(-3.5,3.5, 1), main=list(label="0%", cex=.75), ylab=list(label="Stem Specific Length", cex=.75), xlab=list(label="Mean Leaf thickness", cex=.75))
comp2.10=levelplot(V3~log.mean.thick*log.ssl, data=comp2.10.df, col.regions=diverge_hcl(50),at=seq(-3.5,3.5, 1), main=list(label="10%", cex=.75),xlab=list(label="Mean Leaf Thickness", cex=.75),ylab="")
comp2.20=levelplot(V3~log.mean.thick*log.ssl, data=comp2.20.df,col.regions=diverge_hcl(50),at=seq(-3.5,3.5, 1), main=list(label="20%", cex=.75),xlab=list(label="Mean Leaf Thickness", cex=.75),ylab="")
comp2.30=levelplot(V3~log.mean.thick*log.ssl, data=comp2.30.df, col.regions=diverge_hcl(50),at=seq(-3.5,3.5, 1), main=list(label="30%", cex=.75),xlab=list(label="Mean Leaf Thickness", cex=.75),ylab="")
comp2.40=levelplot(V3~log.mean.thick*log.ssl, data=comp2.40.df, col.regions=diverge_hcl(50),at=seq(-3.5,3.5, 1), main=list(label="40%", cex=.75),ylab=list(label="Stem Specific Length", cex=.75), xlab=list(label="Mean Leaf Thickness", cex=.75))
comp2.med=levelplot(V3~log.mean.thick*log.ssl, data=comp2.med.df, col.regions=diverge_hcl(50),at=seq(-3.5,3.5, 1), main=list(label="50%", cex=.75),xlab=list(label="Mean Leaf Thickness", cex=.75),ylab="")
comp2.60=levelplot(V3~log.mean.thick*log.ssl, data=comp2.60.df, col.regions=diverge_hcl(50),at=seq(-3.5,3.5, 1), main=list(label="60%", cex=.75),xlab=list(label="Mean Leaf Thickness", cex=.75),ylab="")
comp2.70=levelplot(V3~log.mean.thick*log.ssl, data=comp2.70.df, col.regions=diverge_hcl(50),at=seq(-3.5,3.5, 1), main=list(label="70%", cex=.75),xlab=list(label="Mean Leaf Thickness", cex=.75),ylab="")
comp2.80=levelplot(V3~log.mean.thick*log.ssl, data=comp2.80.df, col.regions=diverge_hcl(50),at=seq(-3.5,3.5, 1), main=list(label="80%", cex=.75),ylab=list(label="Stem Specific Length", cex=.75), xlab=list(label="Mean Leaf Thickness", cex=.75))
comp2.90=levelplot(V3~log.mean.thick*log.ssl, data=comp2.90.df, col.regions=diverge_hcl(50),at=seq(-3.5,3.5, 1),main=list(label="90%", cex=.75),xlab=list(label="Mean Leaf Thickness", cex=.75), ylab="")
comp2.high=levelplot(V3~log.mean.thick*log.ssl, data=comp2.high.df, col.regions=diverge_hcl(50),at=seq(-3.5,3.5, 1), main=list(label="100%", cex=.75),xlab=list(label="Mean Leaf Thickness", cex=.75), ylab="")
grid.arrange(comp2.low,comp2.10,comp2.20,comp2.30,comp2.40,comp2.med,comp2.60,comp2.70,comp2.80,comp2.90,comp2.high, nrow=3,ncol=4)

# Figure S6

# mean leaf thickness x LAR x soil PC2

# Vector of the range of soil PC2 by 0.1

vec=seq(-5.4,4.0, by=0.1)

# 10% quantiles for the range of soil PC2

quantile(vec, probs=seq(0,1,0.1))

#Expand the range of the variables in the model while holding soil PC2 constant at each 10% of its range

data.low=with(all.data,expand.grid(log.mean.thick=seq(-2.1,3.2,length.out=10),log.lar1=seq(-6.3,2.5,length.out=10), Comp.2=-5.40))
data.10=with(all.data,expand.grid(log.mean.thick=seq(-2.1,3.2,length.out=10),log.lar1=seq(-6.3,2.5,length.out=10), Comp.2=-4.46))
data.20=with(all.data,expand.grid(log.mean.thick=seq(-2.1,3.2,length.out=10),log.lar1=seq(-6.3,2.5,length.out=10), Comp.2=-3.52))
data.30=with(all.data,expand.grid(log.mean.thick=seq(-2.1,3.2,length.out=10),log.lar1=seq(-6.3,2.5,length.out=10), Comp.2=-2.58))
data.40=with(all.data,expand.grid(log.mean.thick=seq(-2.1,3.2,length.out=10),log.lar1=seq(-6.3,2.5,length.out=10), Comp.2=-1.64))
data.med=with(all.data,expand.grid(log.mean.thick=seq(-2.1,3.2,length.out=10),log.lar1=seq(-6.3,2.5,length.out=10), Comp.2=-0.70))
data.60=with(all.data,expand.grid(log.mean.thick=seq(-2.1,3.2,length.out=10),log.lar1=seq(-6.3,2.5,length.out=10), Comp.2=0.24))
data.70=with(all.data,expand.grid(log.mean.thick=seq(-2.1,3.2,length.out=10),log.lar1=seq(-6.3,2.5,length.out=10), Comp.2=1.18))
data.80=with(all.data,expand.grid(log.mean.thick=seq(-2.1,3.2,length.out=10),log.lar1=seq(-6.3,2.5,length.out=10), Comp.2=2.12))
data.90=with(all.data,expand.grid(log.mean.thick=seq(-2.1,3.2,length.out=10),log.lar1=seq(-6.3,2.5,length.out=10), Comp.2=3.06))
data.high=with(all.data,expand.grid(log.mean.thick=seq(-2.1,3.2,length.out=10),log.lar1=seq(-6.3,2.5,length.out=10), Comp.2=4.0))

# Predict new values from the fitted model with the expanded data

p.low=predict(thick.lar1.pc2,data.low, re.form=NA)
p.10=predict(thick.lar1.pc2,data.10, re.form=NA)
p.20=predict(thick.lar1.pc2,data.20, re.form=NA)
p.30=predict(thick.lar1.pc2,data.30, re.form=NA)
p.40=predict(thick.lar1.pc2,data.40, re.form=NA)
p.med=predict(thick.lar1.pc2,data.med, re.form=NA)
p.60=predict(thick.lar1.pc2,data.60, re.form=NA)
p.70=predict(thick.lar1.pc2,data.70, re.form=NA)
p.80=predict(thick.lar1.pc2,data.80, re.form=NA)
p.90=predict(thick.lar1.pc2,data.90, re.form=NA)
p.high=predict(thick.lar1.pc2,data.high, re.form=NA)

# Put the two trait values and the predicted values into a data.frame

comp2.low.df=data.low[,1:2]
comp2.low.df[,3]=p.low
comp2.10.df=data.10[,1:2]
comp2.10.df[,3]=p.10
comp2.20.df=data.20[,1:2]
comp2.20.df[,3]=p.20
comp2.30.df=data.30[,1:2]
comp2.30.df[,3]=p.30
comp2.40.df=data.40[,1:2]
comp2.40.df[,3]=p.40
comp2.med.df=data.med[,1:2]
comp2.med.df[,3]=p.med
comp2.60.df=data.60[,1:2]
comp2.60.df[,3]=p.60
comp2.70.df=data.70[,1:2]
comp2.70.df[,3]=p.70
comp2.80.df=data.80[,1:2]
comp2.80.df[,3]=p.80
comp2.90.df=data.90[,1:2]
comp2.90.df[,3]=p.90
comp2.high.df=data.high[,1:2]
comp2.high.df[,3]=p.high

# Figure S6

comp2.low=levelplot(V3~log.mean.thick*log.lar1, data=comp2.low.df, col.regions=diverge_hcl(50),at=seq(-5.5,3.0, 1), main=list(label="0%", cex=.75), ylab=list(label="Leaf Area Ratio", cex=.75), xlab=list(label="Mean Leaf Thickness", cex=.75))
comp2.10=levelplot(V3~log.mean.thick*log.lar1, data=comp2.10.df, col.regions=diverge_hcl(50),at=seq(-5.5,3.0, 1), main=list(label="10%", cex=.75),xlab=list(label="Mean Leaf Thickness", cex=.75),ylab="")
comp2.20=levelplot(V3~log.mean.thick*log.lar1, data=comp2.20.df,col.regions=diverge_hcl(50),at=seq(-5.5,3.0, 1),  main=list(label="20%", cex=.75),xlab=list(label="Mean Leaf Thickness", cex=.75),ylab="")
comp2.30=levelplot(V3~log.mean.thick*log.lar1, data=comp2.30.df, col.regions=diverge_hcl(50),at=seq(-5.5,3.0, 1), main=list(label="30%", cex=.75),xlab=list(label="Mean Leaf Thickness", cex=.75),ylab="")
comp2.40=levelplot(V3~log.mean.thick*log.lar1, data=comp2.40.df, col.regions=diverge_hcl(50),at=seq(-5.5,3.0, 1),  main=list(label="40%", cex=.75),ylab=list(label="Leaf Area Ratio", cex=.75), xlab=list(label="Mean Leaf Thickness", cex=.75))
comp2.med=levelplot(V3~log.mean.thick*log.lar1, data=comp2.med.df, col.regions=diverge_hcl(50),at=seq(-5.5,3.0, 1), main=list(label="50%", cex=.75),xlab=list(label="Mean Leaf Thickness", cex=.75),ylab="")
comp2.60=levelplot(V3~log.mean.thick*log.lar1, data=comp2.60.df, col.regions=diverge_hcl(50),at=seq(-5.5,3.0, 1), main=list(label="60%", cex=.75),xlab=list(label="Mean Leaf Thickness", cex=.75),ylab="")
comp2.70=levelplot(V3~log.mean.thick*log.lar1, data=comp2.70.df, col.regions=diverge_hcl(50),at=seq(-5.5,3.0, 1),  main=list(label="70%", cex=.75),xlab=list(label="Mean Leaf Thickness", cex=.75),ylab="")
comp2.80=levelplot(V3~log.mean.thick*log.lar1, data=comp2.80.df, col.regions=diverge_hcl(50),at=seq(-5.5,3.0, 1),   main=list(label="80%", cex=.75),ylab=list(label="Leaf Area Ratio", cex=.75), xlab=list(label="Mean Leaf Thickness", cex=.75))
comp2.90=levelplot(V3~log.mean.thick*log.lar1, data=comp2.90.df, col.regions=diverge_hcl(50),at=seq(-5.5,3.0, 1), main=list(label="90%", cex=.75),xlab=list(label="Mean Leaf Thickness", cex=.75), ylab="")
comp2.high=levelplot(V3~log.mean.thick*log.lar1, data=comp2.high.df, col.regions=diverge_hcl(50),at=seq(-5.5,3.0, 1),main=list(label="100%", cex=.75),xlab=list(label="Mean Leaf Thickness", cex=.75), ylab="")
grid.arrange(comp2.low,comp2.10,comp2.20,comp2.30,comp2.40,comp2.med,comp2.60,comp2.70,comp2.80,comp2.90,comp2.high, nrow=3,ncol=4)

# Figure S7

# LAR x LMF x soil PC2

# Vector of the range of soil PC2 by 0.1

vec=seq(-5.4,4.0, by=0.1)

# 10% quantiles for the range of soil PC2

quantile(vec, probs=seq(0,1,0.1))

#Expand the range of the variables in the model while holding soil PC2 constant at each 10% of its range

data.low=with(all.data,expand.grid(log.lar1=seq(-6.3,2.5,length.out=10),log.lmf=seq(-6.0,1.7,length.out=10), Comp.2=-5.40))
data.10=with(all.data,expand.grid(log.lar1=seq(-6.3,2.5,length.out=10),log.lmf=seq(-6.0,1.7,length.out=10), Comp.2=-4.46))
data.20=with(all.data,expand.grid(log.lar1=seq(-6.3,2.5,length.out=10),log.lmf=seq(-6.0,1.7,length.out=10), Comp.2=-3.52))
data.30=with(all.data,expand.grid(log.lar1=seq(-6.3,2.5,length.out=10),log.lmf=seq(-6.0,1.7,length.out=10), Comp.2=-2.58))
data.40=with(all.data,expand.grid(log.lar1=seq(-6.3,2.5,length.out=10),log.lmf=seq(-6.0,1.7,length.out=10), Comp.2=-1.64))
data.med=with(all.data,expand.grid(log.lar1=seq(-6.3,2.5,length.out=10),log.lmf=seq(-6.0,1.7,length.out=10), Comp.2=-0.70))
data.60=with(all.data,expand.grid(log.lar1=seq(-6.3,2.5,length.out=10),log.lmf=seq(-6.0,1.7,length.out=10), Comp.2=0.24))
data.70=with(all.data,expand.grid(log.lar1=seq(-6.3,2.5,length.out=10),log.lmf=seq(-6.0,1.7,length.out=10), Comp.2=1.18))
data.80=with(all.data,expand.grid(log.lar1=seq(-6.3,2.5,length.out=10),log.lmf=seq(-6.0,1.7,length.out=10),Comp.2=2.12))
data.90=with(all.data,expand.grid(log.lar1=seq(-6.3,2.5,length.out=10),log.lmf=seq(-6.0,1.7,length.out=10), Comp.2=3.06))
data.high=with(all.data,expand.grid(log.lar1=seq(-6.3,2.5,length.out=10),log.lmf=seq(-6.0,1.7,length.out=10), Comp.2=4.0))

# Predict new values from the fitted model with the expanded data

p.low=predict(lar.lmf.pc2,data.low, re.form=NA)
p.10=predict(lar.lmf.pc2,data.10, re.form=NA)
p.20=predict(lar.lmf.pc2,data.20, re.form=NA)
p.30=predict(lar.lmf.pc2,data.30, re.form=NA)
p.40=predict(lar.lmf.pc2,data.40, re.form=NA)
p.med=predict(lar.lmf.pc2,data.med, re.form=NA)
p.60=predict(lar.lmf.pc2,data.60, re.form=NA)
p.70=predict(lar.lmf.pc2,data.70, re.form=NA)
p.80=predict(lar.lmf.pc2,data.80, re.form=NA)
p.90=predict(lar.lmf.pc2,data.90, re.form=NA)
p.high=predict(lar.lmf.pc2,data.high, re.form=NA)
#Put the two trait values and the predicted values into a data.frame

comp2.low.df=data.low[,1:2]
comp2.low.df[,3]=p.low
comp2.10.df=data.10[,1:2]
comp2.10.df[,3]=p.10
comp2.20.df=data.20[,1:2]
comp2.20.df[,3]=p.20
comp2.30.df=data.30[,1:2]
comp2.30.df[,3]=p.30
comp2.40.df=data.40[,1:2]
comp2.40.df[,3]=p.40
comp2.med.df=data.med[,1:2]
comp2.med.df[,3]=p.med
comp2.60.df=data.60[,1:2]
comp2.60.df[,3]=p.60
comp2.70.df=data.70[,1:2]
comp2.70.df[,3]=p.70
comp2.80.df=data.80[,1:2]
comp2.80.df[,3]=p.80
comp2.90.df=data.90[,1:2]
comp2.90.df[,3]=p.90
comp2.high.df=data.high[,1:2]
comp2.high.df[,3]=p.high

# Figure S7

comp2.low=levelplot(V3~log.lar1*log.lmf, data=comp2.low.df, col.regions=diverge_hcl(50),at=seq(-5,3, 1),  main=list(label="0%", cex=.75), ylab=list(label="Leaf Mass Fraction", cex=.75), xlab=list(label="Leaf Area Ratio", cex=.75))
comp2.10=levelplot(V3~log.lar1*log.lmf, data=comp2.10.df, col.regions=diverge_hcl(50),at=seq(-5,3, 1), main=list(label="10%", cex=.75),xlab=list(label="Leaf Area Ratio", cex=.75),ylab="")
comp2.20=levelplot(V3~log.lar1*log.lmf, data=comp2.20.df,col.regions=diverge_hcl(50),at=seq(-5,3, 1),  main=list(label="20%", cex=.75),xlab=list(label="Leaf Area Ratio", cex=.75),ylab="")
comp2.30=levelplot(V3~log.lar1*log.lmf, data=comp2.30.df, col.regions=diverge_hcl(50),at=seq(-5,3, 1),  main=list(label="30%", cex=.75),xlab=list(label="Leaf Area Ratio", cex=.75),ylab="")
comp2.40=levelplot(V3~log.lar1*log.lmf, data=comp2.40.df, col.regions=diverge_hcl(50),at=seq(-5,3, 1),  main=list(label="40%", cex=.75),ylab=list(label="Leaf Mass Fraction", cex=.75), xlab=list(label="Leaf Area Ratio", cex=.75))
comp2.med=levelplot(V3~log.lar1*log.lmf, data=comp2.med.df, col.regions=diverge_hcl(50),at=seq(-5,3, 1),  main=list(label="50%", cex=.75),xlab=list(label="Leaf Area Ratio", cex=.75),ylab="")
comp2.60=levelplot(V3~log.lar1*log.lmf, data=comp2.60.df, col.regions=diverge_hcl(50),at=seq(-5,3, 1),main=list(label="60%", cex=.75),xlab=list(label="Leaf Area Ratio", cex=.75),ylab="")
comp2.70=levelplot(V3~log.lar1*log.lmf, data=comp2.70.df, col.regions=diverge_hcl(50),at=seq(-5,3, 1),  main=list(label="70%", cex=.75),xlab=list(label="Leaf Area Ratio", cex=.75),ylab="")
comp2.80=levelplot(V3~log.lar1*log.lmf, data=comp2.80.df, col.regions=diverge_hcl(50),at=seq(-5,3, 1),   main=list(label="80%", cex=.75),ylab=list(label="Leaf Mass Fraction", cex=.75), xlab=list(label="Leaf Area Ratio", cex=.75))
comp2.90=levelplot(V3~log.lar1*log.lmf, data=comp2.90.df, col.regions=diverge_hcl(50),at=seq(-5,3, 1),main=list(label="90%", cex=.75),xlab=list(label="Leaf Area Ratios", cex=.75), ylab="")
comp2.high=levelplot(V3~log.lar1*log.lmf, data=comp2.high.df, col.regions=diverge_hcl(50),at=seq(-5,3, 1), main=list(label="100%", cex=.75),xlab=list(label="Leaf Area Ratio", cex=.75), ylab="")
grid.arrange(comp2.low,comp2.10,comp2.20,comp2.30,comp2.40,comp2.med,comp2.60,comp2.70,comp2.80,comp2.90,comp2.high, nrow=3,ncol=4)


# Figure S8

#  LMA x RMF x soil PC1

# Vector of the range of soil PC1 by 0.3

vec=seq(-6.0,6.1, by=0.3)

# 10% quantiles for the range of soil PC1

quantile(vec, probs=seq(0,1,0.1))

# Expand the range of the variables in the model while holding soil PC1 constant at each 10% of its range

data.low=with(all.data,expand.grid(log.lma=seq(-3.5,3.2,length.out=10),log.rmf=seq(-5.5,3.4,length.out=10), Comp.1=-6.0))
data.10=with(all.data,expand.grid(log.lma=seq(-3.5,3.2,length.out=10),log.rmf=seq(-5.5,3.4,length.out=10), Comp.1=-4.8))
data.20=with(all.data,expand.grid(log.lma=seq(-3.5,3.2,length.out=10),log.rmf=seq(-5.5,3.4,length.out=10), Comp.1=-3.6))
data.30=with(all.data,expand.grid(log.lma=seq(-3.5,3.2,length.out=10),log.rmf=seq(-5.5,3.4,length.out=10), Comp.1=-2.4))
data.40=with(all.data,expand.grid(log.lma=seq(-3.5,3.2,length.out=10),log.rmf=seq(-5.5,3.4,length.out=10), Comp.1=-1.2))
data.med=with(all.data,expand.grid(log.lma=seq(-3.5,3.2,length.out=10),log.rmf=seq(-5.5,3.4,length.out=10), Comp.1=0.0))
data.60=with(all.data,expand.grid(log.lma=seq(-3.5,3.2,length.out=10),log.rmf=seq(-5.5,3.4,length.out=10), Comp.1=1.2))
data.70=with(all.data,expand.grid(log.lma=seq(-3.5,3.2,length.out=10),log.rmf=seq(-5.5,3.4,length.out=10), Comp.1=2.4))
data.80=with(all.data,expand.grid(log.lma=seq(-3.5,3.2,length.out=10),log.rmf=seq(-5.5,3.4,length.out=10), Comp.1=3.6))
data.90=with(all.data,expand.grid(log.lma=seq(-3.5,3.2,length.out=10),log.rmf=seq(-5.5,3.4,length.out=10), Comp.1=4.8))
data.high=with(all.data,expand.grid(log.lma=seq(-3.5,3.2,length.out=10),log.rmf=seq(-5.5,3.4,length.out=10), Comp.1=6.1))

# Predict new values from the fitted model with the expanded data

p.low=predict(lma.rmf.pc1,data.low, re.form=NA)
p.10=predict(lma.rmf.pc1,data.10, re.form=NA)
p.20=predict(lma.rmf.pc1,data.20, re.form=NA)
p.30=predict(lma.rmf.pc1,data.30, re.form=NA)
p.40=predict(lma.rmf.pc1,data.40, re.form=NA)
p.med=predict(lma.rmf.pc1,data.med, re.form=NA)
p.60=predict(lma.rmf.pc1,data.60, re.form=NA)
p.70=predict(lma.rmf.pc1,data.70, re.form=NA)
p.80=predict(lma.rmf.pc1,data.80, re.form=NA)
p.90=predict(lma.rmf.pc1,data.90, re.form=NA)
p.high=predict(lma.rmf.pc1,data.high, re.form=NA)

# Put the two trait values and the predicted values into a data.frame

comp1.low.df=data.low[,1:2]
comp1.low.df[,3]=p.low
comp1.10.df=data.10[,1:2]
comp1.10.df[,3]=p.10
comp1.20.df=data.20[,1:2]
comp1.20.df[,3]=p.20
comp1.30.df=data.30[,1:2]
comp1.30.df[,3]=p.30
comp1.40.df=data.40[,1:2]
comp1.40.df[,3]=p.40
comp1.med.df=data.med[,1:2]
comp1.med.df[,3]=p.med
comp1.60.df=data.60[,1:2]
comp1.60.df[,3]=p.60
comp1.70.df=data.70[,1:2]
comp1.70.df[,3]=p.70
comp1.80.df=data.80[,1:2]
comp1.80.df[,3]=p.80
comp1.90.df=data.90[,1:2]
comp1.90.df[,3]=p.90
comp1.high.df=data.high[,1:2]
comp1.high.df[,3]=p.high


comp1.low=levelplot(V3~log.lma*log.rmf, data=comp1.low.df, col.regions=diverge_hcl(50),at=seq(-3.5,5.5, 1), main=list(label="0%", cex=.75), ylab=list(label="Root Mass Fraction", cex=.75), xlab=list(label="Leaf Mass per Area", cex=.75))
comp1.10=levelplot(V3~log.lma*log.rmf, data=comp1.10.df, col.regions=diverge_hcl(50),at=seq(-3.5,5.5, 1),main=list(label="10%", cex=.75),xlab=list(label="Leaf Mass per Area", cex=.75), ylab="")
comp1.20=levelplot(V3~log.lma*log.rmf, data=comp1.20.df,col.regions=diverge_hcl(50),at=seq(-3.5,5.5, 1), main=list(label="20%", cex=.75),xlab=list(label="Leaf Mass per Area", cex=.75),ylab="")
comp1.30=levelplot(V3~log.lma*log.rmf, data=comp1.30.df, col.regions=diverge_hcl(50),at=seq(-3.5,5.5, 1), main=list(label="30%", cex=.75),xlab=list(label="Leaf Mass per Area", cex=.75),ylab="")
comp1.40=levelplot(V3~log.lma*log.rmf, data=comp1.40.df, col.regions=diverge_hcl(50),at=seq(-3.5,5.5, 1), main=list(label="40%", cex=.75),ylab=list(label="Root Mass Fraction", cex=.75), xlab=list(label="Leaf Mass per Area", cex=.75))
comp1.med=levelplot(V3~log.lma*log.rmf, data=comp1.med.df, col.regions=diverge_hcl(50),at=seq(-3.5,5.5, 1),main=list(label="50%", cex=.75),xlab=list(label="Leaf Mass per Area", cex=.75),ylab="")
comp1.60=levelplot(V3~log.lma*log.rmf, data=comp1.60.df, col.regions=diverge_hcl(50),at=seq(-3.5,5.5, 1), main=list(label="60%", cex=.75),xlab=list(label="Leaf Mass per Area", cex=.75),ylab="")
comp1.70=levelplot(V3~log.lma*log.rmf, data=comp1.70.df, col.regions=diverge_hcl(50),at=seq(-3.5,5.5, 1), main=list(label="70%", cex=.75),xlab=list(label="Leaf Mass per Area", cex=.75),ylab="")
comp1.80=levelplot(V3~log.lma*log.rmf, data=comp1.80.df, col.regions=diverge_hcl(50),at=seq(-3.5,5.5, 1), main=list(label="80%", cex=.75),ylab=list(label="Root Mass Fraction", cex=.75), xlab=list(label="Leaf Mass per Area", cex=.75))
comp1.90=levelplot(V3~log.lma*log.rmf, data=comp1.90.df, col.regions=diverge_hcl(50),at=seq(-3.5,5.5, 1), main=list(label="90%", cex=.75),xlab=list(label="Leaf Mass per Area", cex=.75), ylab="")
comp1.high=levelplot(V3~log.lma*log.rmf, data=comp1.high.df, col.regions=diverge_hcl(50),at=seq(-3.5,5.5, 1), main=list(label="100%", cex=.75),xlab=list(label="Leaf Mass per Area", cex=.75), ylab="")
grid.arrange(comp1.low,comp1.10,comp1.20,comp1.30,comp1.40,comp1.med,comp1.60,comp1.70,comp1.80,comp1.90,comp1.high, nrow=3,ncol=4)


# Figure S9

#  LMA x mean leaf thickness x soil PC1

# Vector of the range of soil PC1 by 0.3

vec=seq(-6.0,6.1, by=0.3)

# 10% quantiles for the range of soil PC1

quantile(vec, probs=seq(0,1,0.1))

# Expand the range of the variables in the model while holding soil PC1 constant at each 10% of its range

data.low=with(all.data,expand.grid(log.lma=seq(-3.5,3.2,length.out=10),log.mean.thick=seq(-2.1,3.2,length.out=10), Comp.1=-6.0))
data.10=with(all.data,expand.grid(log.lma=seq(-3.5,3.2,length.out=10),log.mean.thick=seq(-2.1,3.2,length.out=10), Comp.1=-4.8))
data.20=with(all.data,expand.grid(log.lma=seq(-3.5,3.2,length.out=10),log.mean.thick=seq(-2.1,3.2,length.out=10), Comp.1=-3.6))
data.30=with(all.data,expand.grid(log.lma=seq(-3.5,3.2,length.out=10),log.mean.thick=seq(-2.1,3.2,length.out=10), Comp.1=-2.4))
data.40=with(all.data,expand.grid(log.lma=seq(-3.5,3.2,length.out=10),log.mean.thick=seq(-2.1,3.2,length.out=10), Comp.1=-1.2))
data.med=with(all.data,expand.grid(log.lma=seq(-3.5,3.2,length.out=10),log.mean.thick=seq(-2.1,3.2,length.out=10), Comp.1=0.0))
data.60=with(all.data,expand.grid(log.lma=seq(-3.5,3.2,length.out=10),log.mean.thick=seq(-2.1,3.2,length.out=10), Comp.1=1.2))
data.70=with(all.data,expand.grid(log.lma=seq(-3.5,3.2,length.out=10),log.mean.thick=seq(-2.1,3.2,length.out=10), Comp.1=2.4))
data.80=with(all.data,expand.grid(log.lma=seq(-3.5,3.2,length.out=10),log.mean.thick=seq(-2.1,3.2,length.out=10), Comp.1=3.6))
data.90=with(all.data,expand.grid(log.lma=seq(-3.5,3.2,length.out=10),log.mean.thick=seq(-2.1,3.2,length.out=10), Comp.1=4.8))
data.high=with(all.data,expand.grid(log.lma=seq(-3.5,3.2,length.out=10),log.mean.thick=seq(-2.1,3.2,length.out=10), Comp.1=6.1))

# Predict new values from the fitted model with the expanded data

p.low=predict(lma.mean.thick.pc1,data.low, re.form=NA)
p.10=predict(lma.mean.thick.pc1,data.10, re.form=NA)
p.20=predict(lma.mean.thick.pc1,data.20, re.form=NA)
p.30=predict(lma.mean.thick.pc1,data.30, re.form=NA)
p.40=predict(lma.mean.thick.pc1,data.40, re.form=NA)
p.med=predict(lma.mean.thick.pc1,data.med, re.form=NA)
p.60=predict(lma.mean.thick.pc1,data.60, re.form=NA)
p.70=predict(lma.mean.thick.pc1,data.70, re.form=NA)
p.80=predict(lma.mean.thick.pc1,data.80, re.form=NA)
p.90=predict(lma.mean.thick.pc1,data.90, re.form=NA)
p.high=predict(lma.mean.thick.pc1,data.high, re.form=NA)

# Put the two trait values and the predicted values into a data.frame

comp1.low.df=data.low[,1:2]
comp1.low.df[,3]=p.low
comp1.10.df=data.10[,1:2]
comp1.10.df[,3]=p.10
comp1.20.df=data.20[,1:2]
comp1.20.df[,3]=p.20
comp1.30.df=data.30[,1:2]
comp1.30.df[,3]=p.30
comp1.40.df=data.40[,1:2]
comp1.40.df[,3]=p.40
comp1.med.df=data.med[,1:2]
comp1.med.df[,3]=p.med
comp1.60.df=data.60[,1:2]
comp1.60.df[,3]=p.60
comp1.70.df=data.70[,1:2]
comp1.70.df[,3]=p.70
comp1.80.df=data.80[,1:2]
comp1.80.df[,3]=p.80
comp1.90.df=data.90[,1:2]
comp1.90.df[,3]=p.90
comp1.high.df=data.high[,1:2]
comp1.high.df[,3]=p.high


comp1.low=levelplot(V3~log.mean.thick*log.rmf, data=comp1.low.df, col.regions=diverge_hcl(50),at=seq(-4.5,4.5, 1), main=list(label="0%", cex=.75), ylab=list(label="Root Mass Fraction", cex=.75), xlab=list(label="Mean Leaf Thickness", cex=.75))
comp1.10=levelplot(V3~log.mean.thick*log.rmf, data=comp1.10.df, col.regions=diverge_hcl(50),at=seq(-4.5,4.5, 1),main=list(label="10%", cex=.75),xlab=list(label="Mean Leaf Thickness", cex=.75), ylab="")
comp1.20=levelplot(V3~log.mean.thick*log.rmf, data=comp1.20.df,col.regions=diverge_hcl(50),at=seq(-4.5,4.5, 1), main=list(label="20%", cex=.75),xlab=list(label="Mean Leaf Thickness", cex=.75),ylab="")
comp1.30=levelplot(V3~log.mean.thick*log.rmf, data=comp1.30.df, col.regions=diverge_hcl(50),at=seq(-4.5,4.5, 1), main=list(label="30%", cex=.75),xlab=list(label="Mean Leaf Thickness", cex=.75),ylab="")
comp1.40=levelplot(V3~log.mean.thick*log.rmf, data=comp1.40.df, col.regions=diverge_hcl(50),at=seq(-4.5,4.5, 1), main=list(label="40%", cex=.75),ylab=list(label="Root Mass Fraction", cex=.75), xlab=list(label="Mean Leaf Thickness", cex=.75))
comp1.med=levelplot(V3~log.mean.thick*log.rmf, data=comp1.med.df, col.regions=diverge_hcl(50),at=seq(-4.5,4.5, 1),main=list(label="50%", cex=.75),xlab=list(label="Mean Leaf Thickness", cex=.75),ylab="")
comp1.60=levelplot(V3~log.mean.thick*log.rmf, data=comp1.60.df, col.regions=diverge_hcl(50),at=seq(-4.5,4.5, 1), main=list(label="60%", cex=.75),xlab=list(label="Mean Leaf Thickness", cex=.75),ylab="")
comp1.70=levelplot(V3~log.mean.thick*log.rmf, data=comp1.70.df, col.regions=diverge_hcl(50),at=seq(-4.5,4.5, 1), main=list(label="70%", cex=.75),xlab=list(label="Mean Leaf Thickness", cex=.75),ylab="")
comp1.80=levelplot(V3~log.mean.thick*log.rmf, data=comp1.80.df, col.regions=diverge_hcl(50),at=seq(-4.5,4.5, 1), main=list(label="80%", cex=.75),ylab=list(label="Root Mass Fraction", cex=.75), xlab=list(label="Mean Leaf Thickness", cex=.75))
comp1.90=levelplot(V3~log.mean.thick*log.rmf, data=comp1.90.df, col.regions=diverge_hcl(50),at=seq(-4.5,4.5, 1), main=list(label="90%", cex=.75),xlab=list(label="Mean Leaf Thickness", cex=.75), ylab="")
comp1.high=levelplot(V3~log.mean.thick*log.rmf, data=comp1.high.df, col.regions=diverge_hcl(50),at=seq(-4.5,4.5, 1), main=list(label="100%", cex=.75),xlab=list(label="Mean Leaf Thickness", cex=.75), ylab="")
grid.arrange(comp1.low,comp1.10,comp1.20,comp1.30,comp1.40,comp1.med,comp1.60,comp1.70,comp1.80,comp1.90,comp1.high, nrow=3,ncol=4)


# Figure S10

# SSL x SMF x Light

# Vector of the range of light by 0.1

vec=seq(-2.6,3.0, by=0.1)

# 10% quantiles for the range of light

quantile(vec, probs=seq(0,1,0.1))

# Expand the range of the variables in the model while holding light constant at each 10% of its range

data.low=with(all.data,expand.grid(log.ssl=seq(-3.3,2.9,length.out=10),log.smf=seq(-6.5,3.0,length.out=10),log.light=-2.6))
data.10=with(all.data,expand.grid(log.ssl=seq(-3.3,2.9,length.out=10),log.smf=seq(-6.5,3.0,length.out=10),log.light=-2.04))
data.20=with(all.data,expand.grid(log.ssl=seq(-3.3,2.9,length.out=10),log.smf=seq(-6.5,3.0,length.out=10),log.light=-1.48))
data.30=with(all.data,expand.grid(log.ssl=seq(-3.3,2.9,length.out=10),log.smf=seq(-6.5,3.0,length.out=10),log.light=-0.92))
data.40=with(all.data,expand.grid(log.ssl=seq(-3.3,2.9,length.out=10),log.smf=seq(-6.5,3.0,length.out=10),log.light=-0.36))
data.med=with(all.data,expand.grid(log.ssl=seq(-3.3,2.9,length.out=10),log.smf=seq(-6.5,3.0,length.out=10),log.light=0.20))
data.60=with(all.data,expand.grid(log.ssl=seq(-3.3,2.9,length.out=10),log.smf=seq(-6.5,3.0,length.out=10),log.light=0.76))
data.70=with(all.data,expand.grid(log.ssl=seq(-3.3,2.9,length.out=10),log.smf=seq(-6.5,3.0,length.out=10),log.light=1.32))
data.80=with(all.data,expand.grid(log.ssl=seq(-3.3,2.9,length.out=10),log.smf=seq(-6.5,3.0,length.out=10),log.light=1.88))
data.90=with(all.data,expand.grid(log.ssl=seq(-3.3,2.9,length.out=10),log.smf=seq(-6.5,3.0,length.out=10),log.light=2.44))
data.high=with(all.data,expand.grid(log.ssl=seq(-3.3,2.9,length.out=10),log.smf=seq(-6.5,3.0,length.out=10),log.light=3.0))

# Predict new values from the fitted model with the expanded data

p.low=predict(ssl.smf.light,data.low, re.form=NA)
p.10=predict(ssl.smf.light,data.10, re.form=NA)
p.20=predict(ssl.smf.light,data.20, re.form=NA)
p.30=predict(ssl.smf.light,data.30, re.form=NA)
p.40=predict(ssl.smf.light,data.40, re.form=NA)
p.med=predict(ssl.smf.light,data.med, re.form=NA)
p.60=predict(ssl.smf.light,data.60, re.form=NA)
p.70=predict(ssl.smf.light,data.70, re.form=NA)
p.80=predict(ssl.smf.light,data.80, re.form=NA)
p.90=predict(ssl.smf.light,data.90, re.form=NA)
p.high=predict(ssl.smf.light,data.high, re.form=NA)

# Put the two trait values and the predicted values into a data.frame

light.low.df=data.low[,1:2]
light.low.df[,3]=p.low
light.10.df=data.10[,1:2]
light.10.df[,3]=p.10
light.20.df=data.20[,1:2]
light.20.df[,3]=p.20
light.30.df=data.30[,1:2]
light.30.df[,3]=p.30
light.40.df=data.40[,1:2]
light.40.df[,3]=p.40
light.med.df=data.med[,1:2]
light.med.df[,3]=p.med
light.60.df=data.60[,1:2]
light.60.df[,3]=p.60
light.70.df=data.70[,1:2]
light.70.df[,3]=p.70
light.80.df=data.80[,1:2]
light.80.df[,3]=p.80
light.90.df=data.90[,1:2]
light.90.df[,3]=p.90
light.high.df=data.high[,1:2]
light.high.df[,3]=p.high

light.low=levelplot(V3~log.ssl*log.smf, data=light.low.df, col.regions=diverge_hcl(50),at=seq(-3,5, 1), main=list(label="0%", cex=.75), ylab=list(label="Stem Mass Fraction", cex=.75),xlab=list(label="Stem Specific Length", cex=.75))
light.10=levelplot(V3~log.ssl*log.smf, data=light.10.df, col.regions=diverge_hcl(50),at=seq(-3,5, 1), main=list(label="10%", cex=.75), ylab="",xlab=list(label="Stem Specific Length", cex=.75))
light.20=levelplot(V3~log.ssl*log.smf, data=light.20.df,col.regions=diverge_hcl(50),at=seq(-3,5, 1), main=list(label="20%", cex=.75), ylab="",xlab=list(label="Stem Specific Length", cex=.75))
light.30=levelplot(V3~log.ssl*log.smf, data=light.30.df, col.regions=diverge_hcl(50),at=seq(-3,5, 1), main=list(label="30%", cex=.75),ylab="",xlab=list(label="Stem Specific Length", cex=.75))
light.40=levelplot(V3~log.ssl*log.smf, data=light.40.df, col.regions=diverge_hcl(50),at=seq(-3,5, 1), main=list(label="40%", cex=.75),ylab=list(label="Stem Mass Fraction", cex=.75),xlab=list(label="Stem Specific Length", cex=.75))
light.med=levelplot(V3~log.ssl*log.smf, data=light.med.df, col.regions=diverge_hcl(50),at=seq(-3,5, 1), main=list(label="50%", cex=.75),ylab="",xlab=list(label="Stem Specific Length", cex=.75))
light.60=levelplot(V3~log.ssl*log.smf, data=light.60.df, col.regions=diverge_hcl(50),at=seq(-3,5, 1), main=list(label="60%", cex=.75),ylab="",xlab=list(label="Stem Specific Length", cex=.75))
light.70=levelplot(V3~log.ssl*log.smf, data=light.70.df, col.regions=diverge_hcl(50),at=seq(-3,5, 1), main=list(label="70%", cex=.75),ylab="",xlab=list(label="Stem Specific Length", cex=.75))
light.80=levelplot(V3~log.ssl*log.smf, data=light.80.df, col.regions=diverge_hcl(50),at=seq(-3,5, 1), main=list(label="80%", cex=.75),ylab=list(label="Stem Mass Fraction",cex=.75), xlab=list(label="Stem Specific Length", cex=.75))
light.90=levelplot(V3~log.ssl*log.smf, data=light.90.df, col.regions=diverge_hcl(50),at=seq(-3,5, 1), main=list(label="90%", cex=.75), xlab=list(label="Stem Specific Length", cex=.75),ylab="")
light.high=levelplot(V3~log.ssl*log.smf, data=light.high.df, col.regions=diverge_hcl(50),at=seq(-3,5, 1), main=list(label="100%", cex=.75), xlab=list(label="Stem Specific Length", cex=.75),ylab="")
grid.arrange(light.low,light.10,light.20,light.30,light.40,light.med,light.60,light.70,light.80,light.90,light.high, nrow=3,ncol=4)


# 3-way models without species as a random effect

# Three-way interaction models; Table 3

lma.rmf.light=lmer(log.rgr~log.lma+log.rmf+log.light+log.lma*log.rmf*log.light+(1|Plot), all.data)
lma.rmf.pc1=lmer(log.rgr~log.lma+log.rmf+Comp.1+log.lma*log.rmf*Comp.1+(1|Plot), all.data)
lma.rmf.pc2=lmer(log.rgr~log.lma+log.rmf+Comp.2+log.lma*log.rmf*Comp.2+(1|Plot), all.data)
lma.mean.thick.pc1=lmer(log.rgr~log.lma+log.mean.thick+Comp.1+log.lma*log.mean.thick*Comp.1+(1|Plot), all.data)
thick.ssl.pc2=lmer(log.rgr~log.mean.thick+log.ssl+Comp.2+log.mean.thick*log.ssl*Comp.2+(1|Plot), all.data)
thick.lar1.pc2=lmer(log.rgr~log.mean.thick+log.lar1+Comp.2+log.mean.thick*log.lar1*Comp.2+(1|Plot), all.data)
thick.smf.pc1=lmer(log.rgr~log.mean.thick+log.smf+Comp.1+log.mean.thick*log.smf*Comp.1+(1|Plot), all.data)
thick.rmf.pc1=lmer(log.rgr~log.mean.thick+log.rmf+Comp.1+log.mean.thick*log.rmf*Comp.1+(1|Plot), all.data)
lar.lmf.pc2=lmer(log.rgr~log.lar1+log.lmf+Comp.2+log.lar1*log.lmf*Comp.2+(1|Plot), all.data)
ssl.smf.light=lmer(log.rgr~log.ssl+log.smf+log.light+log.ssl*log.smf*log.light+(1|Plot), all.data)

sem.model.fits(list(lma.rmf.light,lma.rmf.pc1,lma.rmf.pc2,lma.mean.thick.pc1,thick.ssl.pc2,thick.lar1.pc2,thick.smf.pc1,thick.rmf.pc1,lar.lmf.pc2,ssl.smf.light))

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
