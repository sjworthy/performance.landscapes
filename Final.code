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
## Make a biplot of the pca output; Figure S8; Table S1

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

lma.rgr=lmer(log.rgr.h~log.lma+(1|Plot)+(1|sp),all.data)
rmf.rgr=lmer(log.rgr.h~log.rmf+(1|Plot)+(1|sp),all.data)
lma.light=lmer(log.lma~log.light+(1|Plot)+(1|sp),all.data)
rmf.light=lmer(log.rmf~log.light+(1|Plot)+(1|sp),all.data)

sem.model.fits(list(lma.rgr, rmf.rgr, lma.light, rmf.light))

# Figure S9

dotplot(ranef(lma.rgr, condVar=TRUE), scales=list(cex=0.2))
dotplot(ranef(rmf.rgr, condVar=TRUE), scales=list(cex=0.2))
dotplot(ranef(lma.light, condVar=TRUE), scales=list(cex=0.2))
dotplot(ranef(rmf.light, condVar=TRUE), scales=list(cex=0.2))

# Table S2

extractAIC(lma.rgr)
extractAIC(rmf.rgr)
extractAIC(lma.light)
extractAIC(rmf.light)

# Two-way interaction models; Table 2

lma.light.model=lmer(log.rgr.h~log.lma*log.light+log.lma*Comp.1+log.lma*Comp.2+log.lma*Comp.3+(1|Plot)+(1|sp), all.data)
rmf.light.model=lmer(log.rgr.h~log.rmf*log.light+log.rmf*Comp.1+log.rmf*Comp.2+log.rmf*Comp.3+(1|Plot)+(1|sp), all.data)

sem.model.fits(list(rmf.light.model,lma.light.model))

# Figure S9

dotplot(ranef(lma.light.model, condVar=TRUE),scales=list(cex=0.2))
dotplot(ranef(rmf.light.model, condVar=TRUE),scales=list(cex=0.2))

# Table S2

extractAIC(lma.light.model)
extractAIC(rmf.light.model)

# Three-way interaction models; Table 3

lma.rmf.light=lmer(log.rgr.h~log.lma*log.rmf*log.light+log.lma*log.rmf*Comp.1+log.lma*log.rmf*Comp.2+log.lma*log.rmf*Comp.3+(1|Plot)+(1|sp), all.data)
lma.rmf.comp2=lmer(log.rgr.h~log.lma*log.rmf*log.light+log.lma*log.rmf*Comp.1+log.lma*log.rmf*Comp.2+log.lma*log.rmf*Comp.3+(1|Plot)+(1|sp), all.data)
thick.ssl.comp2=lmer(log.rgr.h~log.mean.thick*log.ssl*log.light+log.mean.thick*log.ssl*Comp.1+log.mean.thick*log.ssl*Comp.2+log.mean.thick*log.ssl*Comp.3+(1|Plot)+(1|sp),all.data)
thick.lar1.comp2=lmer(log.rgr.h~log.mean.thick*log.lar1*log.light+log.mean.thick*log.lar1*Comp.1+log.mean.thick*log.lar1*Comp.2+log.mean.thick*log.lar1*Comp.3+(1|Plot)+(1|sp),all.data)
thick.smf.comp1=lmer(log.rgr.h~log.mean.thick*log.smf*log.light+log.mean.thick*log.smf*Comp.1+log.mean.thick*log.smf*Comp.2+log.mean.thick*log.smf*Comp.3+(1|Plot)+(1|sp),all.data)
thick.rmf.comp1=lmer(log.rgr.h~log.mean.thick*log.rmf*log.light+log.mean.thick*log.rmf*Comp.1+log.mean.thick*log.rmf*Comp.2+log.mean.thick*log.rmf*Comp.3+(1|Plot)+(1|sp),all.data)
lar1.lmf.comp.2=lmer(log.rgr.h~log.lar1*log.lmf*log.light+log.lar1*log.lmf*Comp.1+log.lar1*log.lmf*Comp.2+log.lar1*log.lmf*Comp.3+(1|Plot)+(1|sp),all.data)

sem.model.fits(list(lma.rmf.light,lma.rmf.comp2,thick.ssl.comp2,thick.lar1.comp2,thick.smf.comp1,thick.rmf.comp1,lar1.lmf.comp.2))

# Figure S9

dotplot(ranef(lma.rmf.light, condVar=TRUE),scales=list(cex=0.2))
dotplot(ranef(lma.rmf.comp2, condVar=TRUE),scales=list(cex=0.2))
dotplot(ranef(thick.ssl.comp2, condVar=TRUE),scales=list(cex=0.2))
dotplot(ranef(thick.lar1.comp2, condVar=TRUE),scales=list(cex=0.2))
dotplot(ranef(thick.smf.comp1, condVar=TRUE),scales=list(cex=0.2))
dotplot(ranef(thick.rmf.comp1, condVar=TRUE),scales=list(cex=0.2))
dotplot(ranef(lar1.lmf.comp.2, condVar=TRUE),scales=list(cex=0.2))

# Table S2

extractAIC(lma.rmf.light)
extractAIC(lma.rmf.comp2)
extractAIC(thick.ssl.comp2)
extractAIC(thick.lar1.comp2)
extractAIC(thick.smf.comp1)
extractAIC(thick.rmf.comp1)
extractAIC(lar1.lmf.comp.2)

# Four-way interaction models; Table S3

lma.rmf=lmer(log.rgr.h~log.lma*log.rmf*log.light*Comp.1+log.lma*log.rmf*log.light*Comp.2+log.lma*log.rmf*log.light*Comp.3+(1|Plot)+(1|sp), all.data)
thick.ssl=lmer(log.rgr.h~log.mean.thick*log.ssl*log.light*Comp.1+log.mean.thick*log.ssl*log.light*Comp.2+log.mean.thick*log.ssl*log.light*Comp.3+(1|Plot)+(1|sp),all.data)
thick.lar1=lmer(log.rgr.h~log.mean.thick*log.lar1*log.light*Comp.1+log.mean.thick*log.lar1*log.light*Comp.2+log.mean.thick*log.lar1*log.light*Comp.3+(1|Plot)+(1|sp),all.data)
thick.smf=lmer(log.rgr.h~log.mean.thick*log.smf*log.light*Comp.1+log.mean.thick*log.smf*log.light*Comp.2+log.mean.thick*log.smf*log.light*Comp.3+(1|Plot)+(1|sp),all.data)
thick.rmf=lmer(log.rgr.h~log.mean.thick*log.rmf*log.light*Comp.1+log.mean.thick*log.rmf*log.light*Comp.2+log.mean.thick*log.rmf*log.light*Comp.3+(1|Plot)+(1|sp),all.data)
lar1.lmf=lmer(log.rgr.h~log.lar1*log.lmf*log.light*Comp.1+log.lar1*log.lmf*log.light*Comp.2+log.lar1*log.lmf*log.light*Comp.3+(1|Plot)+(1|sp),all.data)

sem.model.fits(list(lma.rmf, thick.ssl, thick.lar1, thick.smf, thick.rmf, lar1.lmf))

# Figure S9

dotplot(ranef(lma.rmf, condVar=TRUE),scales=list(cex=0.2))
dotplot(ranef(thick.ssl, condVar=TRUE),scales=list(cex=0.2))
dotplot(ranef(thick.lar1, condVar=TRUE),scales=list(cex=0.2))
dotplot(ranef(thick.smf, condVar=TRUE),scales=list(cex=0.2))
dotplot(ranef(thick.rmf, condVar=TRUE),scales=list(cex=0.2))
dotplot(ranef(lar1.lmf, condVar=TRUE),scales=list(cex=0.2))

# Table S2

extractAIC(lma.rmf)
extractAIC(thick.ssl)
extractAIC(thick.lar1)
extractAIC(thick.smf)
extractAIC(thick.rmf)
extractAIC(lar1.lmf)

# Slopes for two-way interaction models

lma.light.output=coef(summary(lma.light.model))[2,1]+coef(summary(lma.light.model))[7,1]*all.data$log.light
rmf.light.output=coef(summary(rmf.light.model))[2,1]+coef(summary(rmf.light.model))[7,1]*all.data$log.light

# Plotting slopes against environment

plot(all.data$log.light, lma.light.output)
plot(all.data$log.light, rmf.light.output)

# Slopes for three-way interaction models

lma.rmf.light.output=coef(summary(lma.rmf.light))[2,1]+coef(summary(lma.rmf.light))[8,1]*all.data$log.rmf+coef(summary(lma.rmf.light))[9,1]*all.data$log.light+coef(summary(lma.rmf.light))[17,1]*(all.data$log.rmf*all.data$log.light)
thick.ssl.comp.2.output=coef(summary(thick.ssl.comp2))[2,1]+coef(summary(thick.ssl.comp2))[8,1]*all.data$log.ssl+coef(summary(thick.ssl.comp2))[13,1]*all.data$Comp.2+coef(summary(thick.ssl.comp2))[19,1]*(all.data$log.ssl*all.data$Comp.2)
thick.lar1.comp.2.output=coef(summary(thick.lar1.comp2))[2,1]+coef(summary(thick.lar1.comp2))[8,1]*all.data$log.lar1+coef(summary(thick.lar1.comp2))[13,1]*all.data$Comp.2+coef(summary(thick.lar1.comp2))[19,1]*(all.data$log.lar1*all.data$Comp.2)
lma.rmf.comp.2.output=coef(summary(lma.rmf.light))[2,1]+coef(summary(lma.rmf.light))[8,1]*all.data$log.rmf+coef(summary(lma.rmf.light))[13,1]*all.data$Comp.2+coef(summary(lma.rmf.light))[19,1]*(all.data$log.rmf*all.data$Comp.2)
lar1.lmf.comp.2.output=coef(summary(lar1.lmf.comp.2))[2,1]+coef(summary(lar1.lmf.comp.2))[8,1]*all.data$log.lmf+coef(summary(lar1.lmf.comp.2))[13,1]*all.data$Comp.2+coef(summary(lar1.lmf.comp.2))[19,1]*(all.data$log.lmf*all.data$Comp.2)
thick.smf.comp.1.output=coef(summary(thick.smf.comp1))[2,1]+coef(summary(thick.smf.comp1))[8,1]*all.data$log.smf+coef(summary(thick.smf.comp1))[11,1]*all.data$Comp.1+coef(summary(thick.smf.comp1))[18,1]*(all.data$log.smf*all.data$Comp.2)
thick.rmf.comp.1.output=coef(summary(thick.rmf.comp1))[2,1]+coef(summary(thick.rmf.comp1))[8,1]*all.data$log.rmf+coef(summary(thick.rmf.comp1))[11,1]*all.data$Comp.1+coef(summary(thick.rmf.comp1))[18,1]*(all.data$log.rmf*all.data$Comp.2)

# Plotting slopes against environment

plot((all.data$log.rmf*all.data$log.light),lma.rmf.light.output)
plot((all.data$log.ssl*all.data$Comp.2),thick.ssl.comp.2.output)
plot((all.data$log.lar1*all.data$Comp.2),thick.lar1.comp.2.output)
plot((all.data$log.rmf*all.data$Comp.2),lma.rmf.comp.2.output)
plot((all.data$log.lmf*all.data$Comp.2), lar1.lmf.comp.2.output)
plot((all.data$log.smf*all.data$Comp.1),thick.smf.comp.1.output)
plot((all.data$log.rmf*all.data$Comp.1),thick.rmf.comp.1.output)

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

vec=seq(-2.7,3.0, by=0.1)

# 10% quantiles for the range of light

quantile(vec, probs=seq(0,1,0.1))

# Expand the range of the variables in the model while holding light constant at each 10% of its range

data.low=with(all.data,expand.grid(log.lma=seq(-3.5,6.1,length.out=10),log.rmf=seq(-5.6,3.4,length.out=10),log.light=-2.7, Comp.1=seq(-6.1,6.1,length.out=10), Comp.2=seq(-5.5,3.7,length.out=10),Comp.3=seq(-3.3,4.9,length.out=10)))
data.10=with(all.data,expand.grid(log.lma=seq(-3.5,6.1,length.out=10),log.rmf=seq(-5.6,3.4,length.out=10),log.light=-2.13, Comp.1=seq(-6.1,6.1,length.out=10), Comp.2=seq(-5.5,3.7,length.out=10),Comp.3=seq(-3.3,4.9,length.out=10)))
data.20=with(all.data,expand.grid(log.lma=seq(-3.5,6.1,length.out=10),log.rmf=seq(-5.6,3.4,length.out=10),log.light=-1.56, Comp.1=seq(-6.1,6.1,length.out=10), Comp.2=seq(-5.5,3.7,length.out=10),Comp.3=seq(-3.3,4.9,length.out=10)))
data.30=with(all.data,expand.grid(log.lma=seq(-3.5,6.1,length.out=10),log.rmf=seq(-5.6,3.4,length.out=10),log.light=-0.99, Comp.1=seq(-6.1,6.1,length.out=10), Comp.2=seq(-5.5,3.7,length.out=10),Comp.3=seq(-3.3,4.9,length.out=10)))
data.40=with(all.data,expand.grid(log.lma=seq(-3.5,6.1,length.out=10),log.rmf=seq(-5.6,3.4,length.out=10),log.light=-0.42, Comp.1=seq(-6.1,6.1,length.out=10), Comp.2=seq(-5.5,3.7,length.out=10),Comp.3=seq(-3.3,4.9,length.out=10)))
data.med=with(all.data,expand.grid(log.lma=seq(-3.5,6.1,length.out=10),log.rmf=seq(-5.6,3.4,length.out=10),log.light=0.15, Comp.1=seq(-6.1,6.1,length.out=10), Comp.2=seq(-5.5,3.7,length.out=10),Comp.3=seq(-3.3,4.9,length.out=10)))
data.60=with(all.data,expand.grid(log.lma=seq(-3.5,6.1,length.out=10),log.rmf=seq(-5.6,3.4,length.out=10),log.light=0.72, Comp.1=seq(-6.1,6.1,length.out=10), Comp.2=seq(-5.5,3.7,length.out=10),Comp.3=seq(-3.3,4.9,length.out=10)))
data.70=with(all.data,expand.grid(log.lma=seq(-3.5,6.1,length.out=10),log.rmf=seq(-5.6,3.4,length.out=10),log.light=1.29, Comp.1=seq(-6.1,6.1,length.out=10), Comp.2=seq(-5.5,3.7,length.out=10),Comp.3=seq(-3.3,4.9,length.out=10)))
data.80=with(all.data,expand.grid(log.lma=seq(-3.5,6.1,length.out=10),log.rmf=seq(-5.6,3.4,length.out=10),log.light=1.86, Comp.1=seq(-6.1,6.1,length.out=10), Comp.2=seq(-5.5,3.7,length.out=10),Comp.3=seq(-3.3,4.9,length.out=10)))
data.90=with(all.data,expand.grid(log.lma=seq(-3.5,6.1,length.out=10),log.rmf=seq(-5.6,3.4,length.out=10),log.light=2.43, Comp.1=seq(-6.1,6.1,length.out=10), Comp.2=seq(-5.5,3.7,length.out=10),Comp.3=seq(-3.3,4.9,length.out=10)))
data.high=with(all.data,expand.grid(log.lma=seq(-3.5,6.1,length.out=10),log.rmf=seq(-5.6,3.4,length.out=10),log.light=3.0, Comp.1=seq(-6.1,6.1,length.out=10), Comp.2=seq(-5.5,3.7,length.out=10),Comp.3=seq(-3.3,4.9,length.out=10)))

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

light.low=levelplot(rgr~log.lma+log.rmf, data=light.low.df, col.regions=diverge_hcl(25),at=seq(-12.5,12, 1), main=list(label="0%", cex=.75), ylab=list(label="Root Mass Fraction", cex=.75),xlab=list(label="Leaf Mass per Area", cex=.75))
light.10=levelplot(V3~log.lma+log.rmf, data=light.10.df, col.regions=diverge_hcl(25),at=seq(-12.5,12, 1), main=list(label="10%", cex=.75), ylab="",xlab=list(label="Leaf Mass per Area", cex=.75))
light.20=levelplot(V3~log.lma+log.rmf, data=light.20.df,col.regions=diverge_hcl(25),at=seq(-12.5,12, 1), main=list(label="20%", cex=.75), ylab="",xlab=list(label="Leaf Mass per Area", cex=.75))
light.30=levelplot(V3~log.lma+log.rmf, data=light.30.df, col.regions=diverge_hcl(25),at=seq(-12.5,12, 1), main=list(label="30%", cex=.75),ylab="",xlab=list(label="Leaf Mass per Area", cex=.75))
light.40=levelplot(V3~log.lma+log.rmf, data=light.40.df, col.regions=diverge_hcl(25),at=seq(-12.5,12, 1), main=list(label="40%", cex=.75),ylab=list(label="Root Mass Fraction", cex=.75),xlab=list(label="Leaf Mass per Area", cex=.75))
light.med=levelplot(rgr~log.lma+log.rmf, data=light.med.df, col.regions=diverge_hcl(25),at=seq(-12.5,12, 1), main=list(label="50%", cex=.75),ylab="",xlab=list(label="Leaf Mass per Area", cex=.75))
light.60=levelplot(V3~log.lma+log.rmf, data=light.60.df, col.regions=diverge_hcl(25),at=seq(-12.5,12, 1), main=list(label="60%", cex=.75),ylab="",xlab=list(label="Leaf Mass per Area", cex=.75))
light.70=levelplot(V3~log.lma+log.rmf, data=light.70.df, col.regions=diverge_hcl(25),at=seq(-12.5,12, 1), main=list(label="70%", cex=.75),ylab="",xlab=list(label="Leaf Mass per Area", cex=.75))
light.80=levelplot(V3~log.lma+log.rmf, data=light.80.df, col.regions=diverge_hcl(25),at=seq(-12.5,12, 1), main=list(label="80%", cex=.75),ylab=list(label="Root Mass Fraction",cex=.75), xlab=list(label="Leaf Mass per Area", cex=.75))
light.90=levelplot(V3~log.lma+log.rmf, data=light.90.df, col.regions=diverge_hcl(25),at=seq(-12.5,12, 1), main=list(label="90%", cex=.75), xlab=list(label="Leaf Mass per Area", cex=.75),ylab="")
light.high=levelplot(rgr~log.lma+log.rmf, data=light.high.df, col.regions=diverge_hcl(25),at=seq(-12.5,12, 1), main=list(label="100%", cex=.75), xlab=list(label="Leaf Mass per Area", cex=.75),ylab="")
grid.arrange(light.low,light.10,light.20,light.30,light.40,light.med,light.60,light.70,light.80,light.90,light.high, nrow=3,ncol=4)

# Figure S2

#  mean leaf thickness x SMF x soil PC1

# Vector of the range of soil PC1 by 0.3

vec=seq(-6.1,6.1, by=0.3)

# 10% quantiles for the range of soil PC1

quantile(vec, probs=seq(0,1,0.1))

# Expand the range of the variables in the model while holding soil PC1 constant at each 10% of its range

data.low=with(all.data,expand.grid(log.mean.thick=seq(-2.2,3.2,length.out=10),log.smf=seq(-6.6,3.0,length.out=10),log.light=seq(-2.7,3.0,length.out=10), Comp.1=-6.1, Comp.2=seq(-5.5,3.7,length.out=10),Comp.3=seq(-3.3,4.9,length.out=10)))
data.10=with(all.data,expand.grid(log.mean.thick=seq(-2.2,3.2,length.out=10),log.smf=seq(-6.6,3.0,length.out=10),log.light=seq(-2.7,3.0,length.out=10), Comp.1=-4.9, Comp.2=seq(-5.5,3.7,length.out=10),Comp.3=seq(-3.3,4.9,length.out=10)))
data.20=with(all.data,expand.grid(log.mean.thick=seq(-2.2,3.2,length.out=10),log.smf=seq(-6.6,3.0,length.out=10),log.light=seq(-2.7,3.0,length.out=10), Comp.1=-3.7, Comp.2=seq(-5.5,3.7,length.out=10),Comp.3=seq(-3.3,4.9,length.out=10)))
data.30=with(all.data,expand.grid(log.mean.thick=seq(-2.2,3.2,length.out=10),log.smf=seq(-6.6,3.0,length.out=10),log.light=seq(-2.7,3.0,length.out=10), Comp.1=-2.5, Comp.2=seq(-5.5,3.7,length.out=10),Comp.3=seq(-3.3,4.9,length.out=10)))
data.40=with(all.data,expand.grid(log.mean.thick=seq(-2.2,3.2,length.out=10),log.smf=seq(-6.6,3.0,length.out=10),log.light=seq(-2.7,3.0,length.out=10), Comp.1=-1.3, Comp.2=seq(-5.5,3.7,length.out=10),Comp.3=seq(-3.3,4.9,length.out=10)))
data.med=with(all.data,expand.grid(log.mean.thick=seq(-2.2,3.2,length.out=10),log.smf=seq(-6.6,3.0,length.out=10),log.light=seq(-2.7,3.0,length.out=10), Comp.1=-0.1, Comp.2=seq(-5.5,3.7,length.out=10),Comp.3=seq(-3.3,4.9,length.out=10)))
data.60=with(all.data,expand.grid(log.mean.thick=seq(-2.2,3.2,length.out=10),log.smf=seq(-6.6,3.0,length.out=10),log.light=seq(-2.7,3.0,length.out=10), Comp.1=1.1, Comp.2=seq(-5.5,3.7,length.out=10),Comp.3=seq(-3.3,4.9,length.out=10)))
data.70=with(all.data,expand.grid(log.mean.thick=seq(-2.2,3.2,length.out=10),log.smf=seq(-6.6,3.0,length.out=10),log.light=seq(-2.7,3.0,length.out=10), Comp.1=2.3, Comp.2=seq(-5.5,3.7,length.out=10),Comp.3=seq(-3.3,4.9,length.out=10)))
data.80=with(all.data,expand.grid(log.mean.thick=seq(-2.2,3.2,length.out=10),log.smf=seq(-6.6,3.0,length.out=10),log.light=seq(-2.7,3.0,length.out=10), Comp.1=3.5, Comp.2=seq(-5.5,3.7,length.out=10),Comp.3=seq(-3.3,4.9,length.out=10)))
data.90=with(all.data,expand.grid(log.mean.thick=seq(-2.2,3.2,length.out=10),log.smf=seq(-6.6,3.0,length.out=10),log.light=seq(-2.7,3.0,length.out=10), Comp.1=4.7, Comp.2=seq(-5.5,3.7,length.out=10),Comp.3=seq(-3.3,4.9,length.out=10)))
data.high=with(all.data,expand.grid(log.mean.thick=seq(-2.2,3.2,length.out=10),log.smf=seq(-6.6,3.0,length.out=10),log.light=seq(-2.7,3.0,length.out=10), Comp.1=6.1, Comp.2=seq(-5.5,3.7,length.out=10),Comp.3=seq(-3.3,4.9,length.out=10)))

# Predict new values from the fitted model with the expanded data

p.low=predict(thick.smf.comp1,data.low, re.form=NA)
p.10=predict(thick.smf.comp1,data.10, re.form=NA)
p.20=predict(thick.smf.comp1,data.20, re.form=NA)
p.30=predict(thick.smf.comp1,data.30, re.form=NA)
p.40=predict(thick.smf.comp1,data.40, re.form=NA)
p.med=predict(thick.smf.comp1,data.med, re.form=NA)
p.60=predict(thick.smf.comp1,data.60, re.form=NA)
p.70=predict(thick.smf.comp1,data.70, re.form=NA)
p.80=predict(thick.smf.comp1,data.80, re.form=NA)
p.90=predict(thick.smf.comp1,data.90, re.form=NA)
p.high=predict(thick.smf.comp1,data.high, re.form=NA)

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

comp1.low=levelplot(V3~log.mean.thick+log.smf, data=comp1.low.df, col.regions=diverge_hcl(25),at=seq(-2,8, 1), main=list(label="0%", cex=.75), ylab=list(label="Stem Mass Fraction", cex=.75), xlab=list(label="Mean Leaf Thickness", cex=.75))
comp1.10=levelplot(V3~log.mean.thick+log.smf, data=comp1.10.df, col.regions=diverge_hcl(25),at=seq(-2,8, 1),main=list(label="10%", cex=.75),xlab=list(label="Mean Leaf Thickness", cex=.75), ylab="")
comp1.20=levelplot(V3~log.mean.thick+log.smf, data=comp1.20.df,col.regions=diverge_hcl(25),at=seq(-2,8, 1), main=list(label="20%", cex=.75),xlab=list(label="Mean Leaf Thickness", cex=.75),ylab="")
comp1.30=levelplot(V3~log.mean.thick+log.smf, data=comp1.30.df, col.regions=diverge_hcl(25),at=seq(-2,8, 1), main=list(label="30%", cex=.75),xlab=list(label="Mean Leaf Thickness", cex=.75),ylab="")
comp1.40=levelplot(V3~log.mean.thick+log.smf, data=comp1.40.df, col.regions=diverge_hcl(25),at=seq(-2,8, 1), main=list(label="40%", cex=.75),ylab=list(label="Stem Mass Fraction", cex=.75), xlab=list(label="Mean Leaf Thickness", cex=.75))
comp1.med=levelplot(V3~log.mean.thick+log.smf, data=comp1.med.df, col.regions=diverge_hcl(25),at=seq(-2,8, 1),main=list(label="50%", cex=.75),xlab=list(label="Mean Leaf Thickness", cex=.75),ylab="")
comp1.60=levelplot(V3~log.mean.thick+log.smf, data=comp1.60.df, col.regions=diverge_hcl(25),at=seq(-2,8, 1), main=list(label="60%", cex=.75),xlab=list(label="Mean Leaf Thickness", cex=.75),ylab="")
comp1.70=levelplot(V3~log.mean.thick+log.smf, data=comp1.70.df, col.regions=diverge_hcl(25),at=seq(-2,8, 1), main=list(label="70%", cex=.75),xlab=list(label="Mean Leaf Thickness", cex=.75),ylab="")
comp1.80=levelplot(V3~log.mean.thick+log.smf, data=comp1.80.df, col.regions=diverge_hcl(25),at=seq(-2,8, 1), main=list(label="80%", cex=.75),ylab=list(label="Stem Mass Fraction", cex=.75), xlab=list(label="Mean Leaf Thickness", cex=.75))
comp1.90=levelplot(V3~log.mean.thick+log.smf, data=comp1.90.df, col.regions=diverge_hcl(25),at=seq(-2,8, 1), main=list(label="90%", cex=.75),xlab=list(label="Mean Leaf Thickness", cex=.75), ylab="")
comp1.high=levelplot(V3~log.mean.thick+log.smf, data=comp1.high.df, col.regions=diverge_hcl(25),at=seq(-2,8, 1), main=list(label="100%", cex=.75),xlab=list(label="Mean Leaf Thickness", cex=.75), ylab="")
grid.arrange(comp1.low,comp1.10,comp1.20,comp1.30,comp1.40,comp1.med,comp1.60,comp1.70,comp1.80,comp1.90,comp1.high, nrow=3,ncol=4)

# Figure S3

#  mean leaf thickness x RMF x soil PC1

# Vector of the range of soil PC1 by 0.3

vec=seq(-6.1,6.1, by=0.3)

# 10% quantiles for the range of soil PC1

quantile(vec, probs=seq(0,1,0.1))

# Expand the range of the variables in the model while holding soil PC1 constant at each 10% of its range

data.low=with(all.data,expand.grid(log.mean.thick=seq(-2.2,3.2,length.out=10),log.rmf=seq(-5.6,3.4,length.out=10),log.light=seq(-2.7,3.0,length.out=10), Comp.1=-6.1, Comp.2=seq(-5.5,3.7,length.out=10),Comp.3=seq(-3.3,4.9,length.out=10)))
data.10=with(all.data,expand.grid(log.mean.thick=seq(-2.2,3.2,length.out=10),log.rmf=seq(-5.6,3.4,length.out=10),log.light=seq(-2.7,3.0,length.out=10), Comp.1=-4.9, Comp.2=seq(-5.5,3.7,length.out=10),Comp.3=seq(-3.3,4.9,length.out=10)))
data.20=with(all.data,expand.grid(log.mean.thick=seq(-2.2,3.2,length.out=10),log.rmf=seq(-5.6,3.4,length.out=10),log.light=seq(-2.7,3.0,length.out=10), Comp.1=-3.7, Comp.2=seq(-5.5,3.7,length.out=10),Comp.3=seq(-3.3,4.9,length.out=10)))
data.30=with(all.data,expand.grid(log.mean.thick=seq(-2.2,3.2,length.out=10),log.rmf=seq(-5.6,3.4,length.out=10),log.light=seq(-2.7,3.0,length.out=10), Comp.1=-2.5, Comp.2=seq(-5.5,3.7,length.out=10),Comp.3=seq(-3.3,4.9,length.out=10)))
data.40=with(all.data,expand.grid(log.mean.thick=seq(-2.2,3.2,length.out=10),log.rmf=seq(-5.6,3.4,length.out=10),log.light=seq(-2.7,3.0,length.out=10), Comp.1=-1.3, Comp.2=seq(-5.5,3.7,length.out=10),Comp.3=seq(-3.3,4.9,length.out=10)))
data.med=with(all.data,expand.grid(log.mean.thick=seq(-2.2,3.2,length.out=10),log.rmf=seq(-5.6,3.4,length.out=10),log.light=seq(-2.7,3.0,length.out=10), Comp.1=-0.1, Comp.2=seq(-5.5,3.7,length.out=10),Comp.3=seq(-3.3,4.9,length.out=10)))
data.60=with(all.data,expand.grid(log.mean.thick=seq(-2.2,3.2,length.out=10),log.rmf=seq(-5.6,3.4,length.out=10),log.light=seq(-2.7,3.0,length.out=10), Comp.1=1.1, Comp.2=seq(-5.5,3.7,length.out=10),Comp.3=seq(-3.3,4.9,length.out=10)))
data.70=with(all.data,expand.grid(log.mean.thick=seq(-2.2,3.2,length.out=10),log.rmf=seq(-5.6,3.4,length.out=10),log.light=seq(-2.7,3.0,length.out=10), Comp.1=2.3, Comp.2=seq(-5.5,3.7,length.out=10),Comp.3=seq(-3.3,4.9,length.out=10)))
data.80=with(all.data,expand.grid(log.mean.thick=seq(-2.2,3.2,length.out=10),log.rmf=seq(-5.6,3.4,length.out=10),log.light=seq(-2.7,3.0,length.out=10), Comp.1=3.5, Comp.2=seq(-5.5,3.7,length.out=10),Comp.3=seq(-3.3,4.9,length.out=10)))
data.90=with(all.data,expand.grid(log.mean.thick=seq(-2.2,3.2,length.out=10),log.rmf=seq(-5.6,3.4,length.out=10),log.light=seq(-2.7,3.0,length.out=10), Comp.1=4.7, Comp.2=seq(-5.5,3.7,length.out=10),Comp.3=seq(-3.3,4.9,length.out=10)))
data.high=with(all.data,expand.grid(log.mean.thick=seq(-2.2,3.2,length.out=10),log.rmf=seq(-5.6,3.4,length.out=10),log.light=seq(-2.7,3.0,length.out=10), Comp.1=6.1, Comp.2=seq(-5.5,3.7,length.out=10),Comp.3=seq(-3.3,4.9,length.out=10)))

# Predict new values from the fitted model with the expanded data

p.low=predict(thick.rmf.comp1,data.low, re.form=NA)
p.10=predict(thick.rmf.comp1,data.10, re.form=NA)
p.20=predict(thick.rmf.comp1,data.20, re.form=NA)
p.30=predict(thick.rmf.comp1,data.30, re.form=NA)
p.40=predict(thick.rmf.comp1,data.40, re.form=NA)
p.med=predict(thick.rmf.comp1,data.med, re.form=NA)
p.60=predict(thick.rmf.comp1,data.60, re.form=NA)
p.70=predict(thick.rmf.comp1,data.70, re.form=NA)
p.80=predict(thick.rmf.comp1,data.80, re.form=NA)
p.90=predict(thick.rmf.comp1,data.90, re.form=NA)
p.high=predict(thick.rmf.comp1,data.high, re.form=NA)

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

comp1.low=levelplot(V3~log.mean.thick+log.rmf, data=comp1.low.df, col.regions=diverge_hcl(25),at=seq(-2,7, 1), main=list(label="0%", cex=.75), ylab=list(label="Root Mass Fraction", cex=.75), xlab=list(label="Mean Leaf Thickness", cex=.75))
comp1.10=levelplot(V3~log.mean.thick+log.rmf, data=comp1.10.df, col.regions=diverge_hcl(25),at=seq(-2,7, 1),main=list(label="10%", cex=.75),xlab=list(label="Mean Leaf Thickness", cex=.75), ylab="")
comp1.20=levelplot(V3~log.mean.thick+log.rmf, data=comp1.20.df,col.regions=diverge_hcl(25),at=seq(-2,7, 1), main=list(label="20%", cex=.75),xlab=list(label="Mean Leaf Thickness", cex=.75),ylab="")
comp1.30=levelplot(V3~log.mean.thick+log.rmf, data=comp1.30.df, col.regions=diverge_hcl(25),at=seq(-2,7, 1), main=list(label="30%", cex=.75),xlab=list(label="Mean Leaf Thickness", cex=.75),ylab="")
comp1.40=levelplot(V3~log.mean.thick+log.rmf, data=comp1.40.df, col.regions=diverge_hcl(25),at=seq(-2,7, 1), main=list(label="40%", cex=.75),ylab=list(label="Root Mass Fraction", cex=.75), xlab=list(label="Mean Leaf Thickness", cex=.75))
comp1.med=levelplot(V3~log.mean.thick+log.rmf, data=comp1.med.df, col.regions=diverge_hcl(25),at=seq(-2,7, 1),main=list(label="50%", cex=.75),xlab=list(label="Mean Leaf Thickness", cex=.75),ylab="")
comp1.60=levelplot(V3~log.mean.thick+log.rmf, data=comp1.60.df, col.regions=diverge_hcl(25),at=seq(-2,7, 1), main=list(label="60%", cex=.75),xlab=list(label="Mean Leaf Thickness", cex=.75),ylab="")
comp1.70=levelplot(V3~log.mean.thick+log.rmf, data=comp1.70.df, col.regions=diverge_hcl(25),at=seq(-2,7, 1), main=list(label="70%", cex=.75),xlab=list(label="Mean Leaf Thickness", cex=.75),ylab="")
comp1.80=levelplot(V3~log.mean.thick+log.rmf, data=comp1.80.df, col.regions=diverge_hcl(25),at=seq(-2,7, 1), main=list(label="80%", cex=.75),ylab=list(label="Root Mass Fraction", cex=.75), xlab=list(label="Mean Leaf Thickness", cex=.75))
comp1.90=levelplot(V3~log.mean.thick+log.rmf, data=comp1.90.df, col.regions=diverge_hcl(25),at=seq(-2,7, 1), main=list(label="90%", cex=.75),xlab=list(label="Mean Leaf Thickness", cex=.75), ylab="")
comp1.high=levelplot(V3~log.mean.thick+log.rmf, data=comp1.high.df, col.regions=diverge_hcl(25),at=seq(-2,7, 1), main=list(label="100%", cex=.75),xlab=list(label="Mean Leaf Thickness", cex=.75), ylab="")
grid.arrange(comp1.low,comp1.10,comp1.20,comp1.30,comp1.40,comp1.med,comp1.60,comp1.70,comp1.80,comp1.90,comp1.high, nrow=3,ncol=4)

# Figure S4

# LMA x RMF x Soil PC2

# Vector of the range of soil PC2 by 0.1

vec=seq(-5.5,3.7, by=0.1)

# 10% quantiles for the range of soil PC2

quantile(vec, probs=seq(0,1,0.1))

# Expand the range of the variables in the model while holding soil PC2 constant at each 10% of its range

comp.2.data.low=with(all.data,expand.grid(log.lma=seq(-3.5,6.1,length.out=10),log.rmf=seq(-5.6,3.4,length.out=10),log.light=seq(-2.7,3.0,length.out=10), Comp.1=seq(-6.1,6.1,length.out=10), Comp.2=-5.50,Comp.3=seq(-3.3,4.9,length.out=10)))
comp.2.data.10=with(all.data,expand.grid(log.lma=seq(-3.5,6.1,length.out=10),log.rmf=seq(-5.6,3.4,length.out=10),log.light=seq(-2.7,3.0,length.out=10), Comp.1=seq(-6.1,6.1,length.out=10), Comp.2=-4.58,Comp.3=seq(-3.3,4.9,length.out=10)))
comp.2.data.20=with(all.data,expand.grid(log.lma=seq(-3.5,6.1,length.out=10),log.rmf=seq(-5.6,3.4,length.out=10),log.light=seq(-2.7,3.0,length.out=10), Comp.1=seq(-6.1,6.1,length.out=10), Comp.2=-3.66,Comp.3=seq(-3.3,4.9,length.out=10)))
comp.2.data.30=with(all.data,expand.grid(log.lma=seq(-3.5,6.1,length.out=10),log.rmf=seq(-5.6,3.4,length.out=10),log.light=seq(-2.7,3.0,length.out=10), Comp.1=seq(-6.1,6.1,length.out=10), Comp.2=-2.74,Comp.3=seq(-3.3,4.9,length.out=10)))
comp.2.data.40=with(all.data,expand.grid(log.lma=seq(-3.5,6.1,length.out=10),log.rmf=seq(-5.6,3.4,length.out=10),log.light=seq(-2.7,3.0,length.out=10), Comp.1=seq(-6.1,6.1,length.out=10), Comp.2=-1.82,Comp.3=seq(-3.3,4.9,length.out=10)))
comp.2.data.med=with(all.data,expand.grid(log.lma=seq(-3.5,6.1,length.out=10),log.rmf=seq(-5.6,3.4,length.out=10),log.light=seq(-2.7,3.0,length.out=10), Comp.1=seq(-6.1,6.1,length.out=10), Comp.2=-0.90,Comp.3=seq(-3.3,4.9,length.out=10)))
comp.2.data.60=with(all.data,expand.grid(log.lma=seq(-3.5,6.1,length.out=10),log.rmf=seq(-5.6,3.4,length.out=10),log.light=seq(-2.7,3.0,length.out=10), Comp.1=seq(-6.1,6.1,length.out=10), Comp.2=0.02,Comp.3=seq(-3.3,4.9,length.out=10)))
comp.2.data.70=with(all.data,expand.grid(log.lma=seq(-3.5,6.1,length.out=10),log.rmf=seq(-5.6,3.4,length.out=10),log.light=seq(-2.7,3.0,length.out=10), Comp.1=seq(-6.1,6.1,length.out=10), Comp.2=0.94,Comp.3=seq(-3.3,4.9,length.out=10)))
comp.2.data.80=with(all.data,expand.grid(log.lma=seq(-3.5,6.1,length.out=10),log.rmf=seq(-5.6,3.4,length.out=10),log.light=seq(-2.7,3.0,length.out=10), Comp.1=seq(-6.1,6.1,length.out=10), Comp.2=1.86,Comp.3=seq(-3.3,4.9,length.out=10)))
comp.2.data.90=with(all.data,expand.grid(log.lma=seq(-3.5,6.1,length.out=10),log.rmf=seq(-5.6,3.4,length.out=10),log.light=seq(-2.7,3.0,length.out=10), Comp.1=seq(-6.1,6.1,length.out=10), Comp.2=2.78,Comp.3=seq(-3.3,4.9,length.out=10)))
comp.2.data.high=with(all.data,expand.grid(log.lma=seq(-3.5,6.1,length.out=10),log.rmf=seq(-5.6,3.4,length.out=10),log.light=seq(-2.7,3.0,length.out=10), Comp.1=seq(-6.1,6.1,length.out=10), Comp.2=3.70,Comp.3=seq(-3.3,4.9,length.out=10)))

# Predict new values from the fitted model with the expanded data

comp.2.p.low=predict(lma.rmf.comp2,comp.2.data.low, re.form=NA)
comp.2.p.10=predict(lma.rmf.comp2,comp.2.data.10, re.form=NA)
comp.2.p.20=predict(lma.rmf.comp2,comp.2.data.20, re.form=NA)
comp.2.p.30=predict(lma.rmf.comp2,comp.2.data.30, re.form=NA)
comp.2.p.40=predict(lma.rmf.comp2,comp.2.data.40, re.form=NA)
comp.2.p.med=predict(lma.rmf.comp2,comp.2.data.med, re.form=NA)
comp.2.p.60=predict(lma.rmf.comp2,comp.2.data.60, re.form=NA)
comp.2.p.70=predict(lma.rmf.comp2,comp.2.data.70, re.form=NA)
comp.2.p.80=predict(lma.rmf.comp2,comp.2.data.80, re.form=NA)
comp.2.p.90=predict(lma.rmf.comp2,comp.2.data.90, re.form=NA)
comp.2.p.high=predict(lma.rmf.comp2,comp.2.data.high, re.form=NA)

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

comp2.low=levelplot(V3~log.lma+log.rmf, data=comp2.low.df, col.regions=diverge_hcl(25),at=seq(-12.5,12, 1), main=list(label="0%", cex=.75), ylab=list(label="Root Mass Fraction", cex=.75), xlab=list(label="Leaf Mass per Area", cex=.75))
comp2.10=levelplot(V3~log.lma+log.rmf, data=comp2.10.df, col.regions=diverge_hcl(25),at=seq(-12.5,12, 1), main=list(label="10%", cex=.75),xlab=list(label="Leaf Mass per Area", cex=.75),ylab="")
comp2.20=levelplot(V3~log.lma+log.rmf, data=comp2.20.df,col.regions=diverge_hcl(25),at=seq(-12.5,12, 1), main=list(label="20%", cex=.75),xlab=list(label="Leaf Mass per Area", cex=.75),ylab="")
comp2.30=levelplot(V3~log.lma+log.rmf, data=comp2.30.df, col.regions=diverge_hcl(25),at=seq(-12.5,12, 1), main=list(label="30%", cex=.75),xlab=list(label="Leaf Mass per Area", cex=.75),ylab="")
comp2.40=levelplot(V3~log.lma+log.rmf, data=comp2.40.df, col.regions=diverge_hcl(25),at=seq(-12.5,12, 1), main=list(label="40%", cex=.75),ylab=list(label="Root Mass Fraction", cex=.75), xlab=list(label="Leaf Mass per Area", cex=.75))
comp2.med=levelplot(V3~log.lma+log.rmf, data=comp2.med.df, col.regions=diverge_hcl(25),at=seq(-12.5,12, 1), main=list(label="50%", cex=.75),xlab=list(label="Leaf Mass per Area", cex=.75),ylab="")
comp2.60=levelplot(V3~log.lma+log.rmf, data=comp2.60.df, col.regions=diverge_hcl(25),at=seq(-12.5,12, 1), main=list(label="60%", cex=.75),xlab=list(label="Leaf Mass per Area", cex=.75),ylab="")
comp2.70=levelplot(V3~log.lma+log.rmf, data=comp2.70.df, col.regions=diverge_hcl(25),at=seq(-12.5,12, 1), main=list(label="70%", cex=.75),xlab=list(label="Leaf Mass per Area", cex=.75),ylab="")
comp2.80=levelplot(V3~log.lma+log.rmf, data=comp2.80.df, col.regions=diverge_hcl(25),at=seq(-12.5,12, 1),  main=list(label="80%", cex=.75),ylab=list(label="Root Mass Fraction", cex=.75), xlab=list(label="Leaf Mass per Area", cex=.75))
comp2.90=levelplot(V3~log.lma+log.rmf, data=comp2.90.df, col.regions=diverge_hcl(25),at=seq(-12.5,12, 1), main=list(label="90%", cex=.75),xlab=list(label="Leaf Mass per Area", cex=.75), ylab="")
comp2.high=levelplot(V3~log.lma+log.rmf, data=comp2.high.df, col.regions=diverge_hcl(25),at=seq(-12.5,12, 1), main=list(label="100%", cex=.75),xlab=list(label="Leaf Mass per Area", cex=.75), ylab="")
grid.arrange(comp2.low,comp2.10,comp2.20,comp2.30,comp2.40,comp2.med,comp2.60,comp2.70,comp2.80,comp2.90,comp2.high, nrow=3,ncol=4)

# Figure S5

# mean leaf thickness x SSL x soil PC2

# Vector of the range of soil PC2 by 0.1

vec=seq(-5.5,3.7, by=0.1)

# 10% quantiles for the range of soil PC2

quantile(vec, probs=seq(0,1,0.1))

#Expand the range of the variables in the model while holding soil PC2 constant at each 10% of its range

data.low=with(all.data,expand.grid(log.mean.thick=seq(-2.2,3.2,length.out=10),log.ssl=seq(-3.4,2.9,length.out=10),log.light=seq(-2.7,3.0,length.out=10), Comp.1=seq(-6.1,6.1,length.out=10), Comp.2=-5.50,Comp.3=seq(-3.3,4.9,length.out=10)))
data.10=with(all.data,expand.grid(log.mean.thick=seq(-2.2,3.2,length.out=10),log.ssl=seq(-3.4,2.9,length.out=10),log.light=seq(-2.7,3.0,length.out=10), Comp.1=seq(-6.1,6.1,length.out=10), Comp.2=-4.58,Comp.3=seq(-3.3,4.9,length.out=10)))
data.20=with(all.data,expand.grid(log.mean.thick=seq(-2.2,3.2,length.out=10),log.ssl=seq(-3.4,2.9,length.out=10),log.light=seq(-2.7,3.0,length.out=10), Comp.1=seq(-6.1,6.1,length.out=10), Comp.2=-3.66,Comp.3=seq(-3.3,4.9,length.out=10)))
data.30=with(all.data,expand.grid(log.mean.thick=seq(-2.2,3.2,length.out=10),log.ssl=seq(-3.4,2.9,length.out=10),log.light=seq(-2.7,3.0,length.out=10), Comp.1=seq(-6.1,6.1,length.out=10), Comp.2=-2.74,Comp.3=seq(-3.3,4.9,length.out=10)))
data.40=with(all.data,expand.grid(log.mean.thick=seq(-2.2,3.2,length.out=10),log.ssl=seq(-3.4,2.9,length.out=10),log.light=seq(-2.7,3.0,length.out=10), Comp.1=seq(-6.1,6.1,length.out=10), Comp.2=-1.82,Comp.3=seq(-3.3,4.9,length.out=10)))
data.med=with(all.data,expand.grid(log.mean.thick=seq(-2.2,3.2,length.out=10),log.ssl=seq(-3.4,2.9,length.out=10),log.light=seq(-2.7,3.0,length.out=10), Comp.1=seq(-6.1,6.1,length.out=10), Comp.2=-0.90,Comp.3=seq(-3.3,4.9,length.out=10)))
data.60=with(all.data,expand.grid(log.mean.thick=seq(-2.2,3.2,length.out=10),log.ssl=seq(-3.4,2.9,length.out=10),log.light=seq(-2.7,3.0,length.out=10), Comp.1=seq(-6.1,6.1,length.out=10), Comp.2=0.02,Comp.3=seq(-3.3,4.9,length.out=10)))
data.70=with(all.data,expand.grid(log.mean.thick=seq(-2.2,3.2,length.out=10),log.ssl=seq(-3.4,2.9,length.out=10),log.light=seq(-2.7,3.0,length.out=10), Comp.1=seq(-6.1,6.1,length.out=10), Comp.2=0.94,Comp.3=seq(-3.3,4.9,length.out=10)))
data.80=with(all.data,expand.grid(log.mean.thick=seq(-2.2,3.2,length.out=10),log.ssl=seq(-3.4,2.9,length.out=10),log.light=seq(-2.7,3.0,length.out=10), Comp.1=seq(-6.1,6.1,length.out=10), Comp.2=1.86,Comp.3=seq(-3.3,4.9,length.out=10)))
data.90=with(all.data,expand.grid(log.mean.thick=seq(-2.2,3.2,length.out=10),log.ssl=seq(-3.4,2.9,length.out=10),log.light=seq(-2.7,3.0,length.out=10), Comp.1=seq(-6.1,6.1,length.out=10), Comp.2=2.78,Comp.3=seq(-3.3,4.9,length.out=10)))
data.high=with(all.data,expand.grid(log.mean.thick=seq(-2.2,3.2,length.out=10),log.ssl=seq(-3.4,2.9,length.out=10),log.light=seq(-2.7,3.0,length.out=10), Comp.1=seq(-6.1,6.1,length.out=10), Comp.2=3.70,Comp.3=seq(-3.3,4.9,length.out=10)))

# Predict new values from the fitted model with the expanded data

p.low=predict(thick.ssl.comp2,data.low, re.form=NA)
p.10=predict(thick.ssl.comp2,data.10, re.form=NA)
p.20=predict(thick.ssl.comp2,data.20, re.form=NA)
p.30=predict(thick.ssl.comp2,data.30, re.form=NA)
p.40=predict(thick.ssl.comp2,data.40, re.form=NA)
p.med=predict(thick.ssl.comp2,data.med, re.form=NA)
p.60=predict(thick.ssl.comp2,data.60, re.form=NA)
p.70=predict(thick.ssl.comp2,data.70, re.form=NA)
p.80=predict(thick.ssl.comp2,data.80, re.form=NA)
p.90=predict(thick.ssl.comp2,data.90, re.form=NA)
p.high=predict(thick.ssl.comp2,data.high, re.form=NA)

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

comp2.low=levelplot(V3~log.mean.thick+log.ssl, data=comp2.low.df, col.regions=diverge_hcl(25),at=seq(-5.5,6.3, 1), main=list(label="0%", cex=.75), ylab=list(label="Stem Specific Length", cex=.75), xlab=list(label="Mean Leaf thickness", cex=.75))
comp2.10=levelplot(V3~log.mean.thick+log.ssl, data=comp2.10.df, col.regions=diverge_hcl(25),at=seq(-5.5,6.3, 1), main=list(label="10%", cex=.75),xlab=list(label="Mean Leaf Thickness", cex=.75),ylab="")
comp2.20=levelplot(V3~log.mean.thick+log.ssl, data=comp2.20.df,col.regions=diverge_hcl(25),at=seq(-5.5,6.3, 1), main=list(label="20%", cex=.75),xlab=list(label="Mean Leaf Thickness", cex=.75),ylab="")
comp2.30=levelplot(V3~log.mean.thick+log.ssl, data=comp2.30.df, col.regions=diverge_hcl(25),at=seq(-5.5,6.3, 1), main=list(label="30%", cex=.75),xlab=list(label="Mean Leaf Thickness", cex=.75),ylab="")
comp2.40=levelplot(V3~log.mean.thick+log.ssl, data=comp2.40.df, col.regions=diverge_hcl(25),at=seq(-5.5,6.3, 1), main=list(label="40%", cex=.75),ylab=list(label="Stem Specific Length", cex=.75), xlab=list(label="Mean Leaf Thickness", cex=.75))
comp2.med=levelplot(V3~log.mean.thick+log.ssl, data=comp2.med.df, col.regions=diverge_hcl(25),at=seq(-5.5,6.3, 1), main=list(label="50%", cex=.75),xlab=list(label="Mean Leaf Thickness", cex=.75),ylab="")
comp2.60=levelplot(V3~log.mean.thick+log.ssl, data=comp2.60.df, col.regions=diverge_hcl(25),at=seq(-5.5,6.3, 1), main=list(label="60%", cex=.75),xlab=list(label="Mean Leaf Thickness", cex=.75),ylab="")
comp2.70=levelplot(V3~log.mean.thick+log.ssl, data=comp2.70.df, col.regions=diverge_hcl(25),at=seq(-5.5,6.3, 1), main=list(label="70%", cex=.75),xlab=list(label="Mean Leaf Thickness", cex=.75),ylab="")
comp2.80=levelplot(V3~log.mean.thick+log.ssl, data=comp2.80.df, col.regions=diverge_hcl(25),at=seq(-5.5,6.3, 1), main=list(label="80%", cex=.75),ylab=list(label="Stem Specific Length", cex=.75), xlab=list(label="Mean Leaf Thickness", cex=.75))
comp2.90=levelplot(V3~log.mean.thick+log.ssl, data=comp2.90.df, col.regions=diverge_hcl(25),at=seq(-5.5,6.3, 1),main=list(label="90%", cex=.75),xlab=list(label="Mean Leaf Thickness", cex=.75), ylab="")
comp2.high=levelplot(V3~log.mean.thick+log.ssl, data=comp2.high.df, col.regions=diverge_hcl(25),at=seq(-5.5,6.3, 1), main=list(label="100%", cex=.75),xlab=list(label="Mean Leaf Thickness", cex=.75), ylab="")
grid.arrange(comp2.low,comp2.10,comp2.20,comp2.30,comp2.40,comp2.med,comp2.60,comp2.70,comp2.80,comp2.90,comp2.high, nrow=3,ncol=4)

# Figure S6

# mean leaf thickness x LAR x soil PC2

# Vector of the range of soil PC2 by 0.1

vec=seq(-5.5,3.7, by=0.1)

# 10% quantiles for the range of soil PC2

quantile(vec, probs=seq(0,1,0.1))

#Expand the range of the variables in the model while holding soil PC2 constant at each 10% of its range

data.low=with(all.data,expand.grid(log.mean.thick=seq(-2.2,3.2,length.out=10),log.lar1=seq(-6.4,2.5,length.out=10),log.light=seq(-2.7,3.0,length.out=10), Comp.1=seq(-6.1,6.1,length.out=10), Comp.2=-5.50,Comp.3=seq(-3.3,4.9,length.out=10)))
data.10=with(all.data,expand.grid(log.mean.thick=seq(-2.2,3.2,length.out=10),log.lar1=seq(-6.4,2.5,length.out=10),log.light=seq(-2.7,3.0,length.out=10), Comp.1=seq(-6.1,6.1,length.out=10), Comp.2=-4.58,Comp.3=seq(-3.3,4.9,length.out=10)))
data.20=with(all.data,expand.grid(log.mean.thick=seq(-2.2,3.2,length.out=10),log.lar1=seq(-6.4,2.5,length.out=10),log.light=seq(-2.7,3.0,length.out=10), Comp.1=seq(-6.1,6.1,length.out=10), Comp.2=-3.66,Comp.3=seq(-3.3,4.9,length.out=10)))
data.30=with(all.data,expand.grid(log.mean.thick=seq(-2.2,3.2,length.out=10),log.lar1=seq(-6.4,2.5,length.out=10),log.light=seq(-2.7,3.0,length.out=10), Comp.1=seq(-6.1,6.1,length.out=10), Comp.2=-2.74,Comp.3=seq(-3.3,4.9,length.out=10)))
data.40=with(all.data,expand.grid(log.mean.thick=seq(-2.2,3.2,length.out=10),log.lar1=seq(-6.4,2.5,length.out=10),log.light=seq(-2.7,3.0,length.out=10), Comp.1=seq(-6.1,6.1,length.out=10), Comp.2=-1.82,Comp.3=seq(-3.3,4.9,length.out=10)))
data.med=with(all.data,expand.grid(log.mean.thick=seq(-2.2,3.2,length.out=10),log.lar1=seq(-6.4,2.5,length.out=10),log.light=seq(-2.7,3.0,length.out=10), Comp.1=seq(-6.1,6.1,length.out=10), Comp.2=-0.90,Comp.3=seq(-3.3,4.9,length.out=10)))
data.60=with(all.data,expand.grid(log.mean.thick=seq(-2.2,3.2,length.out=10),log.lar1=seq(-6.4,2.5,length.out=10),log.light=seq(-2.7,3.0,length.out=10), Comp.1=seq(-6.1,6.1,length.out=10), Comp.2=0.02,Comp.3=seq(-3.3,4.9,length.out=10)))
data.70=with(all.data,expand.grid(log.mean.thick=seq(-2.2,3.2,length.out=10),log.lar1=seq(-6.4,2.5,length.out=10),log.light=seq(-2.7,3.0,length.out=10), Comp.1=seq(-6.1,6.1,length.out=10), Comp.2=0.94,Comp.3=seq(-3.3,4.9,length.out=10)))
data.80=with(all.data,expand.grid(log.mean.thick=seq(-2.2,3.2,length.out=10),log.lar1=seq(-6.4,2.5,length.out=10),log.light=seq(-2.7,3.0,length.out=10), Comp.1=seq(-6.1,6.1,length.out=10), Comp.2=1.86,Comp.3=seq(-3.3,4.9,length.out=10)))
data.90=with(all.data,expand.grid(log.mean.thick=seq(-2.2,3.2,length.out=10),log.lar1=seq(-6.4,2.5,length.out=10),log.light=seq(-2.7,3.0,length.out=10), Comp.1=seq(-6.1,6.1,length.out=10), Comp.2=2.78,Comp.3=seq(-3.3,4.9,length.out=10)))
data.high=with(all.data,expand.grid(log.mean.thick=seq(-2.2,3.2,length.out=10),log.lar1=seq(-6.4,2.5,length.out=10),log.light=seq(-2.7,3.0,length.out=10), Comp.1=seq(-6.1,6.1,length.out=10), Comp.2=3.70,Comp.3=seq(-3.3,4.9,length.out=10)))

# Predict new values from the fitted model with the expanded data

p.low=predict(thick.lar1.comp2,data.low, re.form=NA)
p.10=predict(thick.lar1.comp2,data.10, re.form=NA)
p.20=predict(thick.lar1.comp2,data.20, re.form=NA)
p.30=predict(thick.lar1.comp2,data.30, re.form=NA)
p.40=predict(thick.lar1.comp2,data.40, re.form=NA)
p.med=predict(thick.lar1.comp2,data.med, re.form=NA)
p.60=predict(thick.lar1.comp2,data.60, re.form=NA)
p.70=predict(thick.lar1.comp2,data.70, re.form=NA)
p.80=predict(thick.lar1.comp2,data.80, re.form=NA)
p.90=predict(thick.lar1.comp2,data.90, re.form=NA)
p.high=predict(thick.lar1.comp2,data.high, re.form=NA)

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

comp2.low=levelplot(V3~log.mean.thick+log.lar1, data=comp2.low.df, col.regions=diverge_hcl(25),at=seq(-9,5, 1), main=list(label="0%", cex=.75), ylab=list(label="Leaf Area Ratio", cex=.75), xlab=list(label="Mean Leaf Thickness", cex=.75))
comp2.10=levelplot(V3~log.mean.thick+log.lar1, data=comp2.10.df, col.regions=diverge_hcl(25),at=seq(-9,5, 1), main=list(label="10%", cex=.75),xlab=list(label="Mean Leaf Thickness", cex=.75),ylab="")
comp2.20=levelplot(V3~log.mean.thick+log.lar1, data=comp2.20.df,col.regions=diverge_hcl(25),at=seq(-9,5, 1),  main=list(label="20%", cex=.75),xlab=list(label="Mean Leaf Thickness", cex=.75),ylab="")
comp2.30=levelplot(V3~log.mean.thick+log.lar1, data=comp2.30.df, col.regions=diverge_hcl(25),at=seq(-9,5, 1), main=list(label="30%", cex=.75),xlab=list(label="Mean Leaf Thickness", cex=.75),ylab="")
comp2.40=levelplot(V3~log.mean.thick+log.lar1, data=comp2.40.df, col.regions=diverge_hcl(25),at=seq(-9,5, 1),  main=list(label="40%", cex=.75),ylab=list(label="Leaf Area Ratio", cex=.75), xlab=list(label="Mean Leaf Thickness", cex=.75))
comp2.med=levelplot(V3~log.mean.thick+log.lar1, data=comp2.med.df, col.regions=diverge_hcl(25),at=seq(-9,5, 1), main=list(label="50%", cex=.75),xlab=list(label="Mean Leaf Thickness", cex=.75),ylab="")
comp2.60=levelplot(V3~log.mean.thick+log.lar1, data=comp2.60.df, col.regions=diverge_hcl(25),at=seq(-9,5, 1), main=list(label="60%", cex=.75),xlab=list(label="Mean Leaf Thickness", cex=.75),ylab="")
comp2.70=levelplot(V3~log.mean.thick+log.lar1, data=comp2.70.df, col.regions=diverge_hcl(25),at=seq(-9,5, 1),  main=list(label="70%", cex=.75),xlab=list(label="Mean Leaf Thickness", cex=.75),ylab="")
comp2.80=levelplot(V3~log.mean.thick+log.lar1, data=comp2.80.df, col.regions=diverge_hcl(25),at=seq(-9,5, 1),   main=list(label="80%", cex=.75),ylab=list(label="Leaf Area Ratio", cex=.75), xlab=list(label="Mean Leaf Thickness", cex=.75))
comp2.90=levelplot(V3~log.mean.thick+log.lar1, data=comp2.90.df, col.regions=diverge_hcl(25),at=seq(-9,5, 1), main=list(label="90%", cex=.75),xlab=list(label="Mean Leaf Thickness", cex=.75), ylab="")
comp2.high=levelplot(V3~log.mean.thick+log.lar1, data=comp2.high.df, col.regions=diverge_hcl(25),at=seq(-9,5, 1),main=list(label="100%", cex=.75),xlab=list(label="Mean Leaf Thickness", cex=.75), ylab="")
grid.arrange(comp2.low,comp2.10,comp2.20,comp2.30,comp2.40,comp2.med,comp2.60,comp2.70,comp2.80,comp2.90,comp2.high, nrow=3,ncol=4)

# Figure S7

# LAR x LMF x soil PC2

# Vector of the range of soil PC2 by 0.1

vec=seq(-5.5,3.7, by=0.1)

# 10% quantiles for the range of soil PC2

quantile(vec, probs=seq(0,1,0.1))

#Expand the range of the variables in the model while holding soil PC2 constant at each 10% of its range

data.low=with(all.data,expand.grid(log.lar1=seq(-6.4,2.5,length.out=10),log.lmf=seq(-6.0,1.7,length.out=10),log.light=seq(-2.7,3.0,length.out=10), Comp.1=seq(-6.1,6.1,length.out=10), Comp.2=-5.50,Comp.3=seq(-3.3,4.9,length.out=10)))
data.10=with(all.data,expand.grid(log.lar1=seq(-6.4,2.5,length.out=10),log.lmf=seq(-6.0,1.7,length.out=10),log.light=seq(-2.7,3.0,length.out=10), Comp.1=seq(-6.1,6.1,length.out=10), Comp.2=-4.58,Comp.3=seq(-3.3,4.9,length.out=10)))
data.20=with(all.data,expand.grid(log.lar1=seq(-6.4,2.5,length.out=10),log.lmf=seq(-6.0,1.7,length.out=10),log.light=seq(-2.7,3.0,length.out=10), Comp.1=seq(-6.1,6.1,length.out=10), Comp.2=-3.66,Comp.3=seq(-3.3,4.9,length.out=10)))
data.30=with(all.data,expand.grid(log.lar1=seq(-6.4,2.5,length.out=10),log.lmf=seq(-6.0,1.7,length.out=10),log.light=seq(-2.7,3.0,length.out=10), Comp.1=seq(-6.1,6.1,length.out=10), Comp.2=-2.74,Comp.3=seq(-3.3,4.9,length.out=10)))
data.40=with(all.data,expand.grid(log.lar1=seq(-6.4,2.5,length.out=10),log.lmf=seq(-6.0,1.7,length.out=10),log.light=seq(-2.7,3.0,length.out=10), Comp.1=seq(-6.1,6.1,length.out=10), Comp.2=-1.82,Comp.3=seq(-3.3,4.9,length.out=10)))
data.med=with(all.data,expand.grid(log.lar1=seq(-6.4,2.5,length.out=10),log.lmf=seq(-6.0,1.7,length.out=10),log.light=seq(-2.7,3.0,length.out=10), Comp.1=seq(-6.1,6.1,length.out=10), Comp.2=-0.90,Comp.3=seq(-3.3,4.9,length.out=10)))
data.60=with(all.data,expand.grid(log.lar1=seq(-6.4,2.5,length.out=10),log.lmf=seq(-6.0,1.7,length.out=10),log.light=seq(-2.7,3.0,length.out=10), Comp.1=seq(-6.1,6.1,length.out=10), Comp.2=0.02,Comp.3=seq(-3.3,4.9,length.out=10)))
data.70=with(all.data,expand.grid(log.lar1=seq(-6.4,2.5,length.out=10),log.lmf=seq(-6.0,1.7,length.out=10),log.light=seq(-2.7,3.0,length.out=10), Comp.1=seq(-6.1,6.1,length.out=10), Comp.2=0.94,Comp.3=seq(-3.3,4.9,length.out=10)))
data.80=with(all.data,expand.grid(log.lar1=seq(-6.4,2.5,length.out=10),log.lmf=seq(-6.0,1.7,length.out=10),log.light=seq(-2.7,3.0,length.out=10), Comp.1=seq(-6.1,6.1,length.out=10), Comp.2=1.86,Comp.3=seq(-3.3,4.9,length.out=10)))
data.90=with(all.data,expand.grid(log.lar1=seq(-6.4,2.5,length.out=10),log.lmf=seq(-6.0,1.7,length.out=10),log.light=seq(-2.7,3.0,length.out=10), Comp.1=seq(-6.1,6.1,length.out=10), Comp.2=2.78,Comp.3=seq(-3.3,4.9,length.out=10)))
data.high=with(all.data,expand.grid(log.lar1=seq(-6.4,2.5,length.out=10),log.lmf=seq(-6.0,1.7,length.out=10),log.light=seq(-2.7,3.0,length.out=10), Comp.1=seq(-6.1,6.1,length.out=10), Comp.2=3.70,Comp.3=seq(-3.3,4.9,length.out=10)))

# Predict new values from the fitted model with the expanded data

p.low=predict(lar1.lmf.comp.2,data.low, re.form=NA)
p.10=predict(lar1.lmf.comp.2,data.10, re.form=NA)
p.20=predict(lar1.lmf.comp.2,data.20, re.form=NA)
p.30=predict(lar1.lmf.comp.2,data.30, re.form=NA)
p.40=predict(lar1.lmf.comp.2,data.40, re.form=NA)
p.med=predict(lar1.lmf.comp.2,data.med, re.form=NA)
p.60=predict(lar1.lmf.comp.2,data.60, re.form=NA)
p.70=predict(lar1.lmf.comp.2,data.70, re.form=NA)
p.80=predict(lar1.lmf.comp.2,data.80, re.form=NA)
p.90=predict(lar1.lmf.comp.2,data.90, re.form=NA)
p.high=predict(lar1.lmf.comp.2,data.high, re.form=NA)

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

comp2.low=levelplot(V3~log.lar1+log.lmf, data=comp2.low.df, col.regions=diverge_hcl(25),at=seq(-6.8,5, 1),  main=list(label="0%", cex=.75), ylab=list(label="Leaf Mass Fraction", cex=.75), xlab=list(label="Leaf Area Ratio", cex=.75))
comp2.10=levelplot(V3~log.lar1+log.lmf, data=comp2.10.df, col.regions=diverge_hcl(25),at=seq(-6.8,5, 1), main=list(label="10%", cex=.75),xlab=list(label="Leaf Area Ratio", cex=.75),ylab="")
comp2.20=levelplot(V3~log.lar1+log.lmf, data=comp2.20.df,col.regions=diverge_hcl(25),at=seq(-6.8,5, 1),  main=list(label="20%", cex=.75),xlab=list(label="Leaf Area Ratio", cex=.75),ylab="")
comp2.30=levelplot(V3~log.lar1+log.lmf, data=comp2.30.df, col.regions=diverge_hcl(25),at=seq(-6.8,5, 1),  main=list(label="30%", cex=.75),xlab=list(label="Leaf Area Ratio", cex=.75),ylab="")
comp2.40=levelplot(V3~log.lar1+log.lmf, data=comp2.40.df, col.regions=diverge_hcl(25),at=seq(-6.8,5, 1),  main=list(label="40%", cex=.75),ylab=list(label="Leaf Mass Fraction", cex=.75), xlab=list(label="Leaf Area Ratio", cex=.75))
comp2.med=levelplot(V3~log.lar1+log.lmf, data=comp2.med.df, col.regions=diverge_hcl(25),at=seq(-6.8,5, 1),  main=list(label="50%", cex=.75),xlab=list(label="Leaf Area Ratio", cex=.75),ylab="")
comp2.60=levelplot(V3~log.lar1+log.lmf, data=comp2.60.df, col.regions=diverge_hcl(25),at=seq(-6.8,5, 1),main=list(label="60%", cex=.75),xlab=list(label="Leaf Area Ratio", cex=.75),ylab="")
comp2.70=levelplot(V3~log.lar1+log.lmf, data=comp2.70.df, col.regions=diverge_hcl(25),at=seq(-6.8,5, 1),  main=list(label="70%", cex=.75),xlab=list(label="Leaf Area Ratio", cex=.75),ylab="")
comp2.80=levelplot(V3~log.lar1+log.lmf, data=comp2.80.df, col.regions=diverge_hcl(25),at=seq(-6.8,5, 1),   main=list(label="80%", cex=.75),ylab=list(label="Leaf Mass Fraction", cex=.75), xlab=list(label="Leaf Area Ratio", cex=.75))
comp2.90=levelplot(V3~log.lar1+log.lmf, data=comp2.90.df, col.regions=diverge_hcl(25),at=seq(-6.8,5, 1),main=list(label="90%", cex=.75),xlab=list(label="Leaf Area Ratios", cex=.75), ylab="")
comp2.high=levelplot(V3~log.lar1+log.lmf, data=comp2.high.df, col.regions=diverge_hcl(25),at=seq(-6.8,5, 1), main=list(label="100%", cex=.75),xlab=list(label="Leaf Area Ratio", cex=.75), ylab="")
grid.arrange(comp2.low,comp2.10,comp2.20,comp2.30,comp2.40,comp2.med,comp2.60,comp2.70,comp2.80,comp2.90,comp2.high, nrow=3,ncol=4)