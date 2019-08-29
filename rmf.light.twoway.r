# 2-way

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
plot(exp(final.data$log.rmf), final.data$log.rgr, type="n", xlab="Root Mass Fraction", ylab="Relative Growth Rate", ylim=c(-3,3))
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


results <- readRDS("rmf.light.ini.size.jags_results.RData")
rmf.light.results=do.call(rbind,results) # combine output from 6 MCMC chains
rmf.light.sample.df=rmf.light.results[,c(2,5)] # extract columns of parameters needed

rmf.light.Low=rmf.light[2,2]+rmf.light[5,2]*-2.69
0.0298893

rmf.light.slopes=matrix(data=NA, nrow=1559, ncol=1000)

for(i in 1:1000){
row.sample=rmf.light.sample.df[sample(nrow(rmf.light.sample.df), 1),]
rmf.light.output=row.sample[1]+row.sample[2]*final.data$log.light
rmf.light.slopes[,i]=rmf.light.output
}

# regression to get slope and intercept value for plotting

rmf.light.mean=rmf.light[2,2]+rmf.light[5,2]*final.data$log.light

rmf.light.mean.lm=lm(rmf.light.mean~final.data$log.light)

lm.1000=matrix(data=NA, ncol=2, nrow=1000)

for(i in 1:1000){
lin.mod=lm(rmf.light.slopes[,i]~final.data$log.light)
lm.1000[i,1]=lin.mod$coefficients[1] # intercept
lm.1000[i,2]=lin.mod$coefficients[2] # slope
}
# Determine CIs for intercept and slope distributions

quant.intercept=quantile(lm.1000[,1], probs=c(0.025, 0.975))
quant.slope=quantile(lm.1000[,2], probs=c(0.025, 0.975))

CI = -0.1099652 -0.0152051 

plot((final.data$log.rmf*final.data$log.light),lma.rmf.light.mean, type="n")
abline(-0.07649, 0.05305, col="black") # mean value
abline(-0.13564211,0.01492716, col="gray") # lower limit
abline(-0.01790535,0.09167430, col="gray") # upper limit

