# Effect plots for 3-way models

## Observed ranges of traits and environmental variables

range(final.data$log.rmf)
-5.53, 3.32
range(final.data$log.lmf)
-5.98, 1.69
range(final.data$log.smf)
-6.51, 2.95
range(final.data$log.lar1)
-6.35, 2.45
range(final.data$log.ssl)
-3.25, 2.84
range(final.data$log.mean.thick)
-2.14, 3.14
range(final.data$log.light)
-2.69, 2.97
range(final.data$Comp.1)
-6.06, 6.02
range(final.data$Comp.2)
-5.41, 3.90

# Read in parameter tables
# 3-way
lma.rmf.light=read.csv("lma.rmf.light.rgr.Parameters.csv", header=T, row.names=1)
lar.lmf.pc2=read.csv("lar.lmf.pc2.rgr.Parameters.csv", header=T, row.names=1)
lma.rmf.pc1=read.csv("lma.rmf.pc1.rgr.Parameters.csv", header=T, row.names=1)
lma.rmf.pc2=read.csv("lma.rmf.pc2.rgr.Parameters.csv", header=T, row.names=1)
lma.smf.light=read.csv("lma.smf.light.rgr.Parameters.csv", header=T, row.names=1)
thick.rmf.pc1=read.csv("thick.rmf.pc1.rgr.Parameters.csv", header=T, row.names=1)
thick.smf.pc1=read.csv("thick.smf.pc1.rgr.Parameters.csv", header=T, row.names=1)
thick.ssl.pc2=read.csv("thick.ssl.pc2.rgr.Parameters.csv", header=T, row.names=1)

# 2-way

lma.light=read.csv("lma.light.ini.size.Parameters.csv", header=T, row.names=1)

#Slope
lma.light.High=lma.light[2,2]+lma.light[5,2]*2.97
-0.107117
lma.light.Low=lma.light[2,2]+lma.light[5,2]*-2.69
0.017969

#Intercept
lma.light.High=lma.light[1,2]+lma.light[3,2]*2.97
0.2514939
lma.light.Low=lma.light[1,2]+lma.light[3,2]*-2.69
-0.2628303

#Plots
plot(final.data$scale.log.lma, final.data$log.rgr, type="n", xlab="Leaf Mass per Area", ylab="Relative Growth Rate", ylim=c(-3,3))
abline(0.2514939,-0.107117, col="black", lwd=4)
abline(-0.2628303,0.017969, col="blue", lwd=4)
legend("topright", legend=c("Light","High", "Low"),
col=c("NA","black", "blue"), lty=1, lwd=4)

xtick=seq(-3,3,1)
xlabels=c("0.05", "0.14", "0.37", "1.00", "2.72", "7.39", "20.09")
ytick=seq(-3,3,1)
ylabels=c("0.05", "0.14", "0.37", "1.00", "2.72", "7.39", "20.09")

plot(final.data$scale.log.lma, final.data$log.rgr, type="n", xlab="Leaf Mass per Area", ylab="Relative Growth Rate", xaxt="n", yaxt="n", cex.lab=1.25)
axis(side=1, at=xtick, labels=xlabels, cex.axis=1.15)
axis(side=2, at=ytick, labels=ylabels, cex.axis=1.15)
abline(0.2514939,-0.107117, col="black", lwd=4)
abline(-0.2628303,0.017969, col="blue", lwd=4)
legend("topright", legend=c("Light","High", "Low"),
col=c("NA","black", "blue"), lty=1, lwd=4)

dotchart(lma.light[2:5,2], labels=lma.light[2:5,7], xlim=c(-.15,.2), cex=1.5)
abline(v=0, lty=2)
lines(x=c(lma.light[2,4], lma.light[2,6]), y=c(1,1), lwd=2)
lines(x=c(lma.light[3,4], lma.light[3,6]), y=c(2,2), lwd=2)
lines(x=c(lma.light[4,4], lma.light[4,6]), y=c(3,3), lwd=2)
lines(x=c(lma.light[5,4], lma.light[5,6]), y=c(4,4), lwd=2)
points(x=lma.light[3,2], y=2, pch=19)
points(x=lma.light[4,2], y=3, pch=19)


rmf.light=read.csv("rmf.light.ini.size.Parameters.csv", header=T, row.names=1)

#Slope
rmf.light.High=rmf.light[2,2]+rmf.light[5,2]*2.97
-0.3265209
rmf.light.Low=rmf.light[2,2]+rmf.light[5,2]*-2.69
0.0298893

#Intercept
rmf.light.High=rmf.light[1,2]+rmf.light[3,2]*2.97
-2.09879
rmf.light.Low=rmf.light[1,2]+rmf.light[3,2]*-2.69
-2.630887

#Plots

xtick=seq(-5,3,1)
xlabels=c("0.01","0.02", "0.05", "0.14", "0.37", "1.00", "2.72", "7.39", "20.09")
ytick=seq(-3,3,1)
ylabels=c("0.05", "0.14", "0.37", "1.00", "2.72", "7.39", "20.09")


plot(final.data$log.rmf, final.data$log.rgr, type="n", xlab="Root Mass Fraction", ylab="Relative Growth Rate", ylim=c(-3,3), xaxt="n", yaxt="n", cex.lab=1.25)
axis(side=1, at=xtick, labels=xlabels, cex.axis=1.15)
axis(side=2, at=ytick, labels=ylabels, cex.axis=1.15)
abline(-2.09879,-0.3265209, col="black", lwd=4)
abline(-2.630887,0.0298893, col="blue", lwd=4)
legend("topright", legend=c("Light","High", "Low"),
col=c("NA","black", "blue"), lty=1, lwd=4)


dotchart(rmf.light[2:5,2], labels=rmf.light[2:5,7], xlim=c(-.2,.2), cex=1.5)
abline(v=0, lty=2)
lines(x=c(rmf.light[2,4], rmf.light[2,6]), y=c(1,1), lwd=2)
lines(x=c(rmf.light[3,4], rmf.light[3,6]), y=c(2,2), lwd=2)
lines(x=c(rmf.light[4,4], rmf.light[4,6]), y=c(3,3), lwd=2)
lines(x=c(rmf.light[5,4], rmf.light[5,6]), y=c(4,4), lwd=2)
points(x=rmf.light[2,2], y=1, pch=19)
points(x=rmf.light[3,2], y=2, pch=19)
points(x=rmf.light[4,2], y=3, pch=19)
points(x=rmf.light[5,2], y=4, pch=19)

# Slopes and Intercept

slope = coef(LMA) + coef(LMA:RMF)*RMF + coef(LMA:Light)*Light + coef(LMA:RMF:Light)*RMF*Light
intercept = coef(alpha) + coef(RMF)*RMF + coef(Light)*Light + coef(RMF:Light)*RMF*Light

# lma.rmf.light 
## Slopes

lma.rmf.light.HH=lma.rmf.light[3,2]+lma.rmf.light[5,2]*3.32+lma.rmf.light[6,2]*2.97+lma.rmf.light[8,2]*(3.32*2.97)
0.4429909
lma.rmf.light.LL=lma.rmf.light[3,2]+lma.rmf.light[5,2]*-5.53+lma.rmf.light[6,2]*-2.69+lma.rmf.light[8,2]*(-5.53*-2.69)
0.6726187
lma.rmf.light.HL=lma.rmf.light[3,2]+lma.rmf.light[5,2]*3.32+lma.rmf.light[6,2]*-2.69+lma.rmf.light[8,2]*(3.32*-2.69)
-0.4383616
lma.rmf.light.LH=lma.rmf.light[3,2]+lma.rmf.light[5,2]*-5.53+lma.rmf.light[6,2]*2.97+lma.rmf.light[8,2]*(-5.53*2.97)
-1.103356

## Intercepts

lma.rmf.light.HH=lma.rmf.light[1,2]+lma.rmf.light[4,2]*3.32+lma.rmf.light[2,2]*2.97+lma.rmf.light[7,2]*(3.32*2.97)
-0.2655351
lma.rmf.light.LL=lma.rmf.light[1,2]+lma.rmf.light[4,2]*-5.53+lma.rmf.light[2,2]*-2.69+lma.rmf.light[7,2]*(-5.53*-2.69)
0.2403631
lma.rmf.light.HL=lma.rmf.light[1,2]+lma.rmf.light[4,2]*3.32+lma.rmf.light[2,2]*-2.69+lma.rmf.light[7,2]*(3.32*-2.69)
0.2158185
lma.rmf.light.LH=lma.rmf.light[1,2]+lma.rmf.light[4,2]*-5.53+lma.rmf.light[2,2]*2.97+lma.rmf.light[7,2]*(-5.53*2.97)
2.521027

# Plot

xtick=seq(-3,3,1)
xlabels=c("0.05", "0.14", "0.37", "1.00", "2.72", "7.39", "20.09")
ytick=seq(-1,4,1)
ylabels=c("0.37", "1.00", "2.72", "7.39", "20.09", "54.60")

plot(final.data$scale.log.lma, final.data$log.rgr, type="n", xlab="Leaf Mass per Area", ylab="Relative Growth Rate", xaxt="n", yaxt="n", cex.lab=1.25)
axis(side=1, at=xtick, labels=xlabels, cex.axis=1.15)
axis(side=2, at=ytick, labels=ylabels, cex.axis=1.15)
abline(-0.2655351,0.4429909, col="black", lwd=4)
abline(2.521027,-1.103356, col="blue", lwd=4)
abline(0.2403631,0.6726187, col="coral", lwd=4)
abline(0.2158185,-0.4383616, col="forest green", lwd=4)
legend("topright", legend=c("RMF  Light","High  High", "Low   High", "Low   Low", "High  Low"),
col=c("NA","black", "blue", "coral", "forest green"), lty=1, lwd=4)


coef.vect=lma.rmf.light[,2]
lower.vect=lma.rmf.light[,4]
upper.vect=lma.rmf.light[,6]
long.names=lma.rmf.light[,7]


dotchart(lma.rmf.light[2:9,2], labels=lma.rmf.light[2:9,7], xlim=c(-.2,.2), cex=1.5)
abline(v=0, lty=2)
lines(x=c(lma.rmf.light[2,4], lma.rmf.light[2,6]), y=c(1,1),lwd=2)
lines(x=c(lma.rmf.light[3,4], lma.rmf.light[3,6]), y=c(2,2),lwd=2)
lines(x=c(lma.rmf.light[4,4], lma.rmf.light[4,6]), y=c(3,3),lwd=2)
lines(x=c(lma.rmf.light[5,4], lma.rmf.light[5,6]), y=c(4,4),lwd=2)
lines(x=c(lma.rmf.light[6,4], lma.rmf.light[6,6]), y=c(5,5),lwd=2)
lines(x=c(lma.rmf.light[7,4], lma.rmf.light[7,6]), y=c(6,6),lwd=2)
lines(x=c(lma.rmf.light[8,4], lma.rmf.light[8,6]), y=c(7,7),lwd=2)
lines(x=c(lma.rmf.light[9,4], lma.rmf.light[9,6]), y=c(8,8),lwd=2)
points(x=lma.rmf.light[2,2], y=1, pch=19)
points(x=lma.rmf.light[3,2], y=2, pch=19)
points(x=lma.rmf.light[4,2], y=3, pch=19)
points(x=lma.rmf.light[7,2], y=6, pch=19)
points(x=lma.rmf.light[8,2], y=7, pch=19)
points(x=lma.rmf.light[9,2], y=8, pch=19)

# lar.lmf.pc2

lar.lmf.pc2.HH=lar.lmf.pc2[3,2]+lar.lmf.pc2[5,2]*1.69+lar.lmf.pc2[6,2]*3.90+lar.lmf.pc2[8,2]*(1.69*3.90)
0.5854976
lar.lmf.pc2.LL=lar.lmf.pc2[3,2]+lar.lmf.pc2[5,2]*-5.98+lar.lmf.pc2[6,2]*-5.41+lar.lmf.pc2[8,2]*(-5.98*-5.41)
0.4544581
lar.lmf.pc2.HL=lar.lmf.pc2[3,2]+lar.lmf.pc2[5,2]*1.69+lar.lmf.pc2[6,2]*-5.41+lar.lmf.pc2[8,2]*(1.69*-5.41)
0.1541132
lar.lmf.pc2.LH=lar.lmf.pc2[3,2]+lar.lmf.pc2[5,2]*-5.98+lar.lmf.pc2[6,2]*3.90+lar.lmf.pc2[8,2]*(-5.98*3.90)
-0.4166339

## Intercepts
lar.lmf.pc2.HH=lar.lmf.pc2[1,2]+lar.lmf.pc2[4,2]*1.69+lar.lmf.pc2[2,2]*3.90+lar.lmf.pc2[7,2]*(1.69*3.90)
1.102746
lar.lmf.pc2.LL=lar.lmf.pc2[1,2]+lar.lmf.pc2[4,2]*-5.98+lar.lmf.pc2[2,2]*-5.41+lar.lmf.pc2[7,2]*(-5.98*-5.41)
1.44388
lar.lmf.pc2.HL=lar.lmf.pc2[1,2]+lar.lmf.pc2[4,2]*1.69+lar.lmf.pc2[2,2]*-5.41+lar.lmf.pc2[7,2]*(1.69*-5.41)
0.8226008
lar.lmf.pc2.LH=lar.lmf.pc2[1,2]+lar.lmf.pc2[4,2]*-5.98+lar.lmf.pc2[2,2]*3.90+lar.lmf.pc2[7,2]*(-5.98*3.90)
-0.226833

# Plot

xtick=seq(-6,2,2)
xlabels=c("0.002", "0.02", "0.14", "1.00", "7.39")
ytick=seq(-1,4,1)
ylabels=c("0.37", "1.00", "2.72", "7.39", "20.09", "54.60")

plot(final.data$log.lar1, final.data$log.rgr, type="n", xlab="Leaf Area Ratio", ylab="Relative Growth Rate", xaxt="n", yaxt="n", cex.lab=1.25)
axis(side=1, at=xtick, labels=xlabels, cex.axis=1.10)
axis(side=2, at=ytick, labels=ylabels, cex.axis=1.10)
abline(1.102746,0.5854976, col="black", lwd=4)
abline(-0.226833,-0.4166339, col="blue", lwd=4)
abline(1.44388,0.4544581, col="coral", lwd=4)
abline(-0.226833,0.1541132, col="forest green", lwd=4)
legend("topright", legend=c("LMF  PC2","High  High", "Low   High", "Low   Low", "High  Low"),
col=c("NA","black", "blue", "coral", "forest green"), lty=1, lwd=4)

dotchart(lar.lmf.pc2[2:9,2], labels=lar.lmf.pc2[2:9,7], xlim=c(-.1,.45), cex=1.5)
abline(v=0, lty=2)
lines(x=c(lar.lmf.pc2[2,4], lar.lmf.pc2[2,6]), y=c(1,1),lwd=2)
lines(x=c(lar.lmf.pc2[3,4], lar.lmf.pc2[3,6]), y=c(2,2),lwd=2)
lines(x=c(lar.lmf.pc2[4,4], lar.lmf.pc2[4,6]), y=c(3,3),lwd=2)
lines(x=c(lar.lmf.pc2[5,4], lar.lmf.pc2[5,6]), y=c(4,4),lwd=2)
lines(x=c(lar.lmf.pc2[6,4], lar.lmf.pc2[6,6]), y=c(5,5),lwd=2)
lines(x=c(lar.lmf.pc2[7,4], lar.lmf.pc2[7,6]), y=c(6,6),lwd=2)
lines(x=c(lar.lmf.pc2[8,4], lar.lmf.pc2[8,6]), y=c(7,7),lwd=2)
lines(x=c(lar.lmf.pc2[9,4], lar.lmf.pc2[9,6]), y=c(8,8),lwd=2)
points(x=lar.lmf.pc2[3,2], y=2, pch=19)
points(x=lar.lmf.pc2[5,2], y=4, pch=19)
points(x=lar.lmf.pc2[8,2], y=7, pch=19)
points(x=lar.lmf.pc2[9,2], y=8, pch=19)

# lma.rmf.pc1

lma.rmf.pc1.HH=lma.rmf.pc1[3,2]+lma.rmf.pc1[5,2]*3.32+lma.rmf.pc1[6,2]*6.02+lma.rmf.pc1[8,2]*(3.32*6.02)
0.4296389
lma.rmf.pc1.LL=lma.rmf.pc1[3,2]+lma.rmf.pc1[5,2]*-5.53+lma.rmf.pc1[6,2]*-6.06+lma.rmf.pc1[8,2]*(-5.53*-6.06)
0.8465451
lma.rmf.pc1.HL=lma.rmf.pc1[3,2]+lma.rmf.pc1[5,2]*3.32+lma.rmf.pc1[6,2]*-6.06+lma.rmf.pc1[8,2]*(3.32*-6.06)
-0.5350603
lma.rmf.pc1.LH=lma.rmf.pc1[3,2]+lma.rmf.pc1[5,2]*-5.53+lma.rmf.pc1[6,2]*6.02+lma.rmf.pc1[8,2]*(-5.53*6.02)
-1.070995

# Intercepts

lma.rmf.pc1.HH=lma.rmf.pc1[1,2]+lma.rmf.pc1[4,2]*3.32+lma.rmf.pc1[2,2]*6.02+lma.rmf.pc1[7,2]*(3.32*6.02)
-1.255323
lma.rmf.pc1.LL=lma.rmf.pc1[1,2]+lma.rmf.pc1[4,2]*-5.53+lma.rmf.pc1[2,2]*-6.06+lma.rmf.pc1[7,2]*(-5.53*-6.06)
0.5149377
lma.rmf.pc1.HL=lma.rmf.pc1[1,2]+lma.rmf.pc1[4,2]*3.32+lma.rmf.pc1[2,2]*-6.06+lma.rmf.pc1[7,2]*(3.32*-6.06)
-0.1132796
lma.rmf.pc1.LH=lma.rmf.pc1[1,2]+lma.rmf.pc1[4,2]*-5.53+lma.rmf.pc1[2,2]*6.02+lma.rmf.pc1[7,2]*(-5.53*6.02)
0.4687015

# Plot

xtick=seq(-3,3,1)
xlabels=c("0.05", "0.14", "0.37", "1.00", "2.72", "7.39", "20.09")
ytick=seq(-1,4,1)
ylabels=c("0.37", "1.00", "2.72", "7.39", "20.09", "54.60")

plot(final.data$scale.log.lma, final.data$log.rgr, type="n", xlab="Leaf Mass per Area", ylab="Relative Growth Rate", xaxt="n", yaxt="n", cex.lab=1.25)
axis(side=1, at=xtick, labels=xlabels, cex.axis=1.10)
axis(side=2, at=ytick, labels=ylabels, cex.axis=1.10)
abline(-1.255323,0.4296389, col="black", lwd=4)
abline(0.4687015,-1.070995, col="blue", lwd=4)
abline(0.5149377,0.8465451, col="coral", lwd=4)
abline(-0.1132796,-0.5350603, col="forest green", lwd=4)
legend("topright", legend=c("RMF  PC1","High  High", "Low   High", "Low   Low", "High  Low"),
col=c("NA","black", "blue", "coral", "forest green"), lty=1, lwd=4)

dotchart(lma.rmf.pc1[2:9,2], labels=lma.rmf.pc1[2:9,7], xlim=c(-.2,.1), cex=1.5)
abline(v=0, lty=2)
lines(x=c(lma.rmf.pc1[2,4], lma.rmf.pc1[2,6]), y=c(1,1), lwd=2)
lines(x=c(lma.rmf.pc1[3,4], lma.rmf.pc1[3,6]), y=c(2,2), lwd=2)
lines(x=c(lma.rmf.pc1[4,4], lma.rmf.pc1[4,6]), y=c(3,3), lwd=2)
lines(x=c(lma.rmf.pc1[5,4], lma.rmf.pc1[5,6]), y=c(4,4), lwd=2)
lines(x=c(lma.rmf.pc1[6,4], lma.rmf.pc1[6,6]), y=c(5,5), lwd=2)
lines(x=c(lma.rmf.pc1[7,4], lma.rmf.pc1[7,6]), y=c(6,6), lwd=2)
lines(x=c(lma.rmf.pc1[8,4], lma.rmf.pc1[8,6]), y=c(7,7), lwd=2)
lines(x=c(lma.rmf.pc1[9,4], lma.rmf.pc1[9,6]), y=c(8,8), lwd=2)
points(x=lma.rmf.pc1[2,2], y=1, pch=19)
points(x=lma.rmf.pc1[3,2], y=2, pch=19)
points(x=lma.rmf.pc1[4,2], y=3, pch=19)
points(x=lma.rmf.pc1[8,2], y=7, pch=19)
points(x=lma.rmf.pc1[9,2], y=8, pch=19)


# lma.rmf.pc2

lma.rmf.pc2.HH=lma.rmf.pc2[3,2]+lma.rmf.pc2[5,2]*3.32+lma.rmf.pc2[5,2]*3.90+lma.rmf.pc2[8,2]*(3.32*3.90)
0.3469226
lma.rmf.pc2.LL=lma.rmf.pc2[3,2]+lma.rmf.pc2[5,2]*-5.53+lma.rmf.pc2[5,2]*-5.41+lma.rmf.pc2[8,2]*(-5.53*-5.41)
0.5992843
lma.rmf.pc2.HL=lma.rmf.pc2[3,2]+lma.rmf.pc2[5,2]*3.32+lma.rmf.pc2[5,2]*-5.41+lma.rmf.pc2[8,2]*(3.32*-5.41)
-0.5527362
lma.rmf.pc2.LH=lma.rmf.pc2[3,2]+lma.rmf.pc2[5,2]*-5.53+lma.rmf.pc2[5,2]*3.90+lma.rmf.pc2[8,2]*(-5.53*3.90)
-0.6416401

# Intercepts

lma.rmf.pc2.HH=lma.rmf.pc2[1,2]+lma.rmf.pc2[4,2]*3.32+lma.rmf.pc2[2,2]*3.90+lma.rmf.pc2[7,2]*(3.32*3.90)
0.7289176
lma.rmf.pc2.LL=lma.rmf.pc2[1,2]+lma.rmf.pc2[4,2]*-5.53+lma.rmf.pc2[2,2]*-5.41+lma.rmf.pc2[7,2]*(-5.53*-5.41)
1.698563
lma.rmf.pc2.HL=lma.rmf.pc2[1,2]+lma.rmf.pc2[4,2]*3.32+lma.rmf.pc2[2,2]*-5.41+lma.rmf.pc2[7,2]*(3.32*-5.41)
1.567559
lma.rmf.pc2.LH=lma.rmf.pc2[1,2]+lma.rmf.pc2[4,2]*-5.53+lma.rmf.pc2[2,2]*3.90+lma.rmf.pc2[7,2]*(-5.53*3.90)
2.798641

# Plot
xtick=seq(-3,3,1)
xlabels=c("0.05", "0.14", "0.37", "1.00", "2.72", "7.39", "20.09")
ytick=seq(-1,4,1)
ylabels=c("0.37", "1.00", "2.72", "7.39", "20.09", "54.60")

plot(final.data$scale.log.lma, final.data$log.rgr, type="n", xlab="Leaf Mass per Area", ylab="Relative Growth Rate", xaxt="n", yaxt="n", cex.lab=1.25)
axis(side=1, at=xtick, labels=xlabels, cex.axis=1.10)
axis(side=2, at=ytick, labels=ylabels, cex.axis=1.10)
abline(0.7289176,0.3469226, col="black", lwd=4)
abline(2.798641,-0.6416401, col="blue", lwd=4)
abline(1.698563,0.5992843, col="coral", lwd=4)
abline(1.567559,-0.5527362, col="forest green", lwd=4)
legend("bottomleft", legend=c("RMF  PC2","High  High", "Low   High", "Low   Low", "High  Low"),
col=c("NA","black", "blue", "coral", "forest green"), lty=1, lwd=4)

dotchart(lma.rmf.pc2[2:9,2], labels=lma.rmf.pc2[2:9,7], xlim=c(-.2,.1), cex=1.5)
abline(v=0, lty=2)
lines(x=c(lma.rmf.pc2[2,4], lma.rmf.pc2[2,6]), y=c(1,1), lwd=2)
lines(x=c(lma.rmf.pc2[3,4], lma.rmf.pc2[3,6]), y=c(2,2), lwd=2)
lines(x=c(lma.rmf.pc2[4,4], lma.rmf.pc2[4,6]), y=c(3,3), lwd=2)
lines(x=c(lma.rmf.pc2[5,4], lma.rmf.pc2[5,6]), y=c(4,4), lwd=2)
lines(x=c(lma.rmf.pc2[6,4], lma.rmf.pc2[6,6]), y=c(5,5), lwd=2)
lines(x=c(lma.rmf.pc2[7,4], lma.rmf.pc2[7,6]), y=c(6,6), lwd=2)
lines(x=c(lma.rmf.pc2[8,4], lma.rmf.pc2[8,6]), y=c(7,7), lwd=2)
lines(x=c(lma.rmf.pc2[9,4], lma.rmf.pc2[9,6]), y=c(8,8), lwd=2)
points(x=lma.rmf.pc2[4,2], y=3, pch=19)
points(x=lma.rmf.pc2[8,2], y=7, pch=19)
points(x=lma.rmf.pc2[9,2], y=8, pch=19)

# lma.smf.light

lma.smf.light.HH=lma.smf.light[3,2]+lma.smf.light[5,2]*2.95+lma.smf.light[6,2]*2.97+lma.smf.light[8,2]*(2.95*2.97)
-0.3306204
lma.smf.light.LL=lma.smf.light[3,2]+lma.smf.light[5,2]*-6.51+lma.smf.light[6,2]*-2.69+lma.smf.light[8,2]*(-6.51*-2.69)
-1.25495
lma.smf.light.HL=lma.smf.light[3,2]+lma.smf.light[5,2]*2.95+lma.smf.light[6,2]*-2.69+lma.smf.light[8,2]*(2.95*-2.69)
0.4548235
lma.smf.light.LH=lma.smf.light[3,2]+lma.smf.light[5,2]*-6.51+lma.smf.light[6,2]*2.97+lma.smf.light[8,2]*(-6.51*2.97)
0.3626425

# Intercepts

lma.smf.light.HH=lma.smf.light[1,2]+lma.smf.light[4,2]*2.95+lma.smf.light[2,2]*2.97+lma.smf.light[7,2]*(2.95*2.97)
1.164414
lma.smf.light.LL=lma.smf.light[1,2]+lma.smf.light[4,2]*-6.51+lma.smf.light[2,2]*-2.69+lma.smf.light[7,2]*(-6.51*-2.69)
2.146015
lma.smf.light.HL=lma.smf.light[1,2]+lma.smf.light[4,2]*2.95+lma.smf.light[2,2]*-2.69+lma.smf.light[7,2]*(2.95*-2.69)
0.8296849
lma.smf.light.LH=lma.smf.light[1,2]+lma.smf.light[4,2]*-6.51+lma.smf.light[2,2]*2.97+lma.smf.light[7,2]*(-6.51*2.97)
2.883874

# Plot

xtick=seq(-3,3,1)
xlabels=c("0.05", "0.14", "0.37", "1.00", "2.72", "7.39", "20.09")
ytick=seq(-1,4,1)
ylabels=c("0.37", "1.00", "2.72", "7.39", "20.09", "54.60")

plot(final.data$scale.log.lma, final.data$log.rgr, type="n", xlab="Leaf Mass per Area", ylab="Relative Growth Rate", xaxt="n", yaxt="n", cex.lab=1.25)
axis(side=1, at=xtick, labels=xlabels, cex.axis=1.10)
axis(side=2, at=ytick, labels=ylabels, cex.axis=1.10)
abline(1.164414,-0.3306204, col="black", lwd=4)
abline(2.883874,0.3626425, col="blue", lwd=4)
abline(2.146015,-1.25495, col="coral", lwd=4)
abline(0.8296849,0.4548235, col="forest green", lwd=4)
legend("bottomleft", legend=c("SMF  Light","High  High", "Low   High", "Low   Low", "High  Low"),
col=c("NA","black", "blue", "coral", "forest green"), lty=1, lwd=4)

dotchart(lma.smf.light[2:9,2], labels=lma.smf.light[2:9,7], xlim=c(-.25,.2), cex=1.5)
abline(v=0, lty=2)
lines(x=c(lma.smf.light[2,4], lma.smf.light[2,6]), y=c(1,1), lwd=2)
lines(x=c(lma.smf.light[3,4], lma.smf.light[3,6]), y=c(2,2), lwd=2)
lines(x=c(lma.smf.light[4,4], lma.smf.light[4,6]), y=c(3,3), lwd=2)
lines(x=c(lma.smf.light[5,4], lma.smf.light[5,6]), y=c(4,4), lwd=2)
lines(x=c(lma.smf.light[6,4], lma.smf.light[6,6]), y=c(5,5), lwd=2)
lines(x=c(lma.smf.light[7,4], lma.smf.light[7,6]), y=c(6,6), lwd=2)
lines(x=c(lma.smf.light[8,4], lma.smf.light[8,6]), y=c(7,7), lwd=2)
lines(x=c(lma.smf.light[9,4], lma.smf.light[9,6]), y=c(8,8), lwd=2)
points(x=lma.smf.light[2,2], y=1, pch=19)
points(x=lma.smf.light[3,2], y=2, pch=19)
points(x=lma.smf.light[4,2], y=3, pch=19)
points(x=lma.smf.light[5,2], y=4, pch=19)
points(x=lma.smf.light[8,2], y=7, pch=19)
points(x=lma.smf.light[9,2], y=8, pch=19)

# thick.rmf.pc1

thick.rmf.pc1.HH=thick.rmf.pc1[3,2]+thick.rmf.pc1[5,2]*3.32+thick.rmf.pc1[6,2]*6.02+thick.rmf.pc1[8,2]*(3.32*6.02)
0.7411207
thick.rmf.pc1.LL=thick.rmf.pc1[3,2]+thick.rmf.pc1[5,2]*-5.53+thick.rmf.pc1[6,2]*-6.06+thick.rmf.pc1[8,2]*(-5.53*-6.06)
1.007932
thick.rmf.pc1.HL=thick.rmf.pc1[3,2]+thick.rmf.pc1[5,2]*3.32+thick.rmf.pc1[6,2]*-6.06+thick.rmf.pc1[8,2]*(3.32*-6.06)
-0.6721572
thick.rmf.pc1.LH=thick.rmf.pc1[3,2]+thick.rmf.pc1[5,2]*-5.53+thick.rmf.pc1[6,2]*6.02+thick.rmf.pc1[8,2]*(-5.53*6.02)
-1.428547

# Intercept

thick.rmf.pc1.HH=thick.rmf.pc1[1,2]+thick.rmf.pc1[4,2]*3.32+thick.rmf.pc1[2,2]*6.02+thick.rmf.pc1[7,2]*(3.32*6.02)
1.291638
thick.rmf.pc1.LL=thick.rmf.pc1[1,2]+thick.rmf.pc1[4,2]*-5.53+thick.rmf.pc1[2,2]*-6.06+thick.rmf.pc1[7,2]*(-5.53*-6.06)
3.091157
thick.rmf.pc1.HL=thick.rmf.pc1[1,2]+thick.rmf.pc1[4,2]*3.32+thick.rmf.pc1[2,2]*-6.06+thick.rmf.pc1[7,2]*(3.32*-6.06)
2.53052
thick.rmf.pc1.LH=thick.rmf.pc1[1,2]+thick.rmf.pc1[4,2]*-5.53+thick.rmf.pc1[2,2]*6.02+thick.rmf.pc1[7,2]*(-5.53*6.02)
2.976948

# Plot

xtick=seq(-2,3,1)
xlabels=c("0.14", "0.37", "1.00", "2.72", "7.39", "20.09")
ytick=seq(-1,4,1)
ylabels=c("0.37", "1.00", "2.72", "7.39", "20.09", "54.60")

plot(final.data$log.mean.thick, final.data$log.rgr, type="n", xlab="Mean Thickness", ylab="Relative Growth Rate", xaxt="n", yaxt="n", cex.lab=1.25)
axis(side=1, at=xtick, labels=xlabels, cex.axis=1.10)
axis(side=2, at=ytick, labels=ylabels, cex.axis=1.10)
abline(1.291638,0.7411207, col="black", lwd=4)
abline(2.976948,-1.428547, col="blue", lwd=4)
abline(3.091157,1.007932, col="coral", lwd=4)
abline(2.53052,-0.6721572, col="forest green", lwd=4)
legend("bottomleft", legend=c("RMF  PC1","High  High", "Low   High", "Low   Low", "High  Low"),
col=c("NA","black", "blue", "coral", "forest green"), lty=1, lwd=4)

dotchart(thick.rmf.pc1[2:9,2], labels=thick.rmf.pc1[2:9,7], xlim=c(-.2,.1), cex=1.5)
abline(v=0, lty=2)
lines(x=c(thick.rmf.pc1[2,4], thick.rmf.pc1[2,6]), y=c(1,1), lwd=2)
lines(x=c(thick.rmf.pc1[3,4], thick.rmf.pc1[3,6]), y=c(2,2), lwd=2)
lines(x=c(thick.rmf.pc1[4,4], thick.rmf.pc1[4,6]), y=c(3,3), lwd=2)
lines(x=c(thick.rmf.pc1[5,4], thick.rmf.pc1[5,6]), y=c(4,4), lwd=2)
lines(x=c(thick.rmf.pc1[6,4], thick.rmf.pc1[6,6]), y=c(5,5), lwd=2)
lines(x=c(thick.rmf.pc1[7,4], thick.rmf.pc1[7,6]), y=c(6,6), lwd=2)
lines(x=c(thick.rmf.pc1[8,4], thick.rmf.pc1[8,6]), y=c(7,7), lwd=2)
lines(x=c(thick.rmf.pc1[9,4], thick.rmf.pc1[9,6]), y=c(8,8), lwd=2)
points(x=thick.rmf.pc1[4,2], y=3, pch=19)
points(x=thick.rmf.pc1[8,2], y=7, pch=19)
points(x=thick.rmf.pc1[9,2], y=8, pch=19)

# thick.smf.pc1

thick.smf.pc1.HH=thick.smf.pc1[3,2]+thick.smf.pc1[5,2]*2.95+thick.smf.pc1[6,2]*6.02+thick.smf.pc1[8,2]*(2.95*6.02)
-0.3758739
thick.smf.pc1.LL=thick.smf.pc1[3,2]+thick.smf.pc1[5,2]*-6.51+thick.smf.pc1[6,2]*-6.06+thick.smf.pc1[8,2]*(-6.51*-6.06)
-1.102249
thick.smf.pc1.HL=thick.smf.pc1[3,2]+thick.smf.pc1[5,2]*2.95+thick.smf.pc1[6,2]*-6.06+thick.smf.pc1[8,2]*(2.95*-6.06)
0.4699858
thick.smf.pc1.LH=thick.smf.pc1[3,2]+thick.smf.pc1[5,2]*-6.51+thick.smf.pc1[6,2]*6.02+thick.smf.pc1[8,2]*(-6.51*6.02)
0.505414

# Intercept

thick.smf.pc1.HH=thick.smf.pc1[1,2]+thick.smf.pc1[4,2]*2.95+thick.smf.pc1[2,2]*6.02+thick.smf.pc1[7,2]*(2.95*6.02)
-0.8197453
thick.smf.pc1.LL=thick.smf.pc1[1,2]+thick.smf.pc1[4,2]*-6.51+thick.smf.pc1[2,2]*-6.06+thick.smf.pc1[7,2]*(-6.51*-6.06)
1.927117
thick.smf.pc1.HL=thick.smf.pc1[1,2]+thick.smf.pc1[4,2]*2.95+thick.smf.pc1[2,2]*-6.06+thick.smf.pc1[7,2]*(2.95*-6.06)
-0.5580381
thick.smf.pc1.LH=thick.smf.pc1[1,2]+thick.smf.pc1[4,2]*-6.51+thick.smf.pc1[2,2]*6.02+thick.smf.pc1[7,2]*(-6.51*6.02)
-0.3561466

#Plot

xtick=seq(-2,3,1)
xlabels=c("0.14", "0.37", "1.00", "2.72", "7.39", "20.09")
ytick=seq(-1,4,1)
ylabels=c("0.37", "1.00", "2.72", "7.39", "20.09", "54.60")

plot(final.data$log.mean.thick, final.data$log.rgr, type="n", xlab="Mean Thickness", ylab="Relative Growth Rate", xaxt="n", yaxt="n", cex.lab=1.25)
axis(side=1, at=xtick, labels=xlabels, cex.axis=1.10)
axis(side=2, at=ytick, labels=ylabels, cex.axis=1.10)
abline(-0.8197453,-0.3758739, col="black", lwd=4)
abline(-0.3561466,0.505414, col="blue", lwd=4)
abline(1.927117,-1.102249, col="coral", lwd=4)
abline(-0.5580381,0.4699858, col="forest green", lwd=4)
legend("topright", legend=c("SMF  PC1","High  High", "Low   High", "Low   Low", "High  Low"),
col=c("NA","black", "blue", "coral", "forest green"), lty=1, lwd=4)

dotchart(thick.smf.pc1[2:9,2], labels=thick.smf.pc1[2:9,7], xlim=c(-.25,.1), cex=1.5)
abline(v=0, lty=2)
lines(x=c(thick.smf.pc1[2,4], thick.smf.pc1[2,6]), y=c(1,1), lwd=2)
lines(x=c(thick.smf.pc1[3,4], thick.smf.pc1[3,6]), y=c(2,2), lwd=2)
lines(x=c(thick.smf.pc1[4,4], thick.smf.pc1[4,6]), y=c(3,3), lwd=2)
lines(x=c(thick.smf.pc1[5,4], thick.smf.pc1[5,6]), y=c(4,4), lwd=2)
lines(x=c(thick.smf.pc1[6,4], thick.smf.pc1[6,6]), y=c(5,5), lwd=2)
lines(x=c(thick.smf.pc1[7,4], thick.smf.pc1[7,6]), y=c(6,6), lwd=2)
lines(x=c(thick.smf.pc1[8,4], thick.smf.pc1[8,6]), y=c(7,7), lwd=2)
lines(x=c(thick.smf.pc1[9,4], thick.smf.pc1[9,6]), y=c(8,8), lwd=2)
points(x=thick.smf.pc1[2,2], y=1, pch=19)
points(x=thick.smf.pc1[4,2], y=3, pch=19)
points(x=thick.smf.pc1[8,2], y=7, pch=19)
points(x=thick.smf.pc1[9,2], y=8, pch=19)

# thick.ssl.pc2

thick.ssl.pc2.HH=thick.ssl.pc2[3,2]+thick.ssl.pc2[5,2]*2.84+thick.ssl.pc2[6,2]*3.90+thick.ssl.pc2[8,2]*(2.84*3.90)
-0.6532705
thick.ssl.pc2.LL=thick.ssl.pc2[3,2]+thick.ssl.pc2[5,2]*-3.25+thick.ssl.pc2[6,2]*-5.41+thick.ssl.pc2[8,2]*(-3.25*-5.41)
-0.7650798
thick.ssl.pc2.HL=thick.ssl.pc2[3,2]+thick.ssl.pc2[5,2]*2.84+thick.ssl.pc2[6,2]*-5.41+thick.ssl.pc2[8,2]*(2.84*-5.41)
0.7147707
thick.ssl.pc2.LH=thick.ssl.pc2[3,2]+thick.ssl.pc2[5,2]*-3.25+thick.ssl.pc2[6,2]*3.90+thick.ssl.pc2[8,2]*(-3.25*3.90)
0.4454995

# Intercepts

thick.ssl.pc2.HH=thick.ssl.pc2[1,2]+thick.ssl.pc2[4,2]*2.84+thick.ssl.pc2[2,2]*3.90+thick.ssl.pc2[7,2]*(2.84*3.90)
3.953492
thick.ssl.pc2.LL=thick.ssl.pc2[1,2]+thick.ssl.pc2[4,2]*-3.25+thick.ssl.pc2[2,2]*-5.41+thick.ssl.pc2[7,2]*(-3.25*-5.41)
3.81898
thick.ssl.pc2.HL=thick.ssl.pc2[1,2]+thick.ssl.pc2[4,2]*2.84+thick.ssl.pc2[2,2]*-5.41+thick.ssl.pc2[7,2]*(2.84*-5.41)
3.923716
thick.ssl.pc2.LH=thick.ssl.pc2[1,2]+thick.ssl.pc2[4,2]*-3.25+thick.ssl.pc2[2,2]*3.90+thick.ssl.pc2[7,2]*(-3.25*3.90)
3.394209

#Plot

xtick=seq(-2,3,1)
xlabels=c("0.14", "0.37", "1.00", "2.72", "7.39", "20.09")
ytick=seq(-1,4,1)
ylabels=c("0.37", "1.00", "2.72", "7.39", "20.09", "54.60")

plot(final.data$log.mean.thick, final.data$log.rgr, type="n", xlab="Mean Thickness", ylab="Relative Growth Rate", xaxt="n", yaxt="n", cex.lab=1.25)
axis(side=1, at=xtick, labels=xlabels, cex.axis=1.10)
axis(side=2, at=ytick, labels=ylabels, cex.axis=1.10)
abline(3.953492,-0.6532705, col="black", lwd=4)
abline(3.394209,0.4454995, col="blue", lwd=4)
abline(3.81898,-0.7650798, col="coral", lwd=4)
abline(3.394209,0.7147707, col="forest green", lwd=4)
legend("bottomleft", legend=c("SSL  PC2","High  High", "Low   High", "Low   Low", "High  Low"),
col=c("NA","black", "blue", "coral", "forest green"), lty=1, lwd=4)

dotchart(thick.ssl.pc2[2:9,2], labels=thick.ssl.pc2[2:9,7], xlim=c(-.15,.15), cex=1.5)
abline(v=0, lty=2)
lines(x=c(thick.ssl.pc2[2,4], thick.ssl.pc2[2,6]), y=c(1,1), lwd=2)
lines(x=c(thick.ssl.pc2[3,4], thick.ssl.pc2[3,6]), y=c(2,2), lwd=2)
lines(x=c(thick.ssl.pc2[4,4], thick.ssl.pc2[4,6]), y=c(3,3), lwd=2)
lines(x=c(thick.ssl.pc2[5,4], thick.ssl.pc2[5,6]), y=c(4,4), lwd=2)
lines(x=c(thick.ssl.pc2[6,4], thick.ssl.pc2[6,6]), y=c(5,5), lwd=2)
lines(x=c(thick.ssl.pc2[7,4], thick.ssl.pc2[7,6]), y=c(6,6), lwd=2)
lines(x=c(thick.ssl.pc2[8,4], thick.ssl.pc2[8,6]), y=c(7,7), lwd=2)
lines(x=c(thick.ssl.pc2[9,4], thick.ssl.pc2[9,6]), y=c(8,8), lwd=2)
points(x=thick.ssl.pc2[8,2], y=7, pch=19)
points(x=thick.ssl.pc2[9,2], y=8, pch=19)