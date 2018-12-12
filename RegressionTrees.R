# Making regression trees of the two-way interactions models

library(rpart)
library(rpart.plot)

mod.lma=rpart(log.rgr~log.lma+log.light, rgdata)
rpart.plot(mod.lma, type=3, digits=3, fallen.leaves=TRUE)
printcp(mod.lma)
new.mod.lma=prune.rpart(mod.lma, 0.010000)
rpart.plot(new.mod.lma, type=3, digits=3, fallen.leaves=TRUE)
1-sum(residuals(new.mod.lma)^2)/sum((rgdata$log.rgr-mean(rgdata$log.rgr))^2)
R2 = 0.03

mod.rmf=rpart(log.rgr~log.rmf+log.light, rgdata)
rpart.plot(mod.rmf, type=3, digits=3, fallen.leaves=TRUE)
printcp(mod.rmf)
new.mod.rmf=prune.rpart(mod.rmf, 0.010000)
rpart.plot(new.mod.rmf, type=3, digits=3, fallen.leaves=TRUE)
1-sum(residuals(new.mod.rmf)^2)/sum((rgdata$log.rgr-mean(rgdata$log.rgr))^2)
R2 = 0.05


# Making regression trees of the three-way interactions models

library(rpart)
library(rpart.plot)

mod.1=rpart(log.rgr~log.lma+log.rmf+log.light, data=rgdata)
rpart.plot(mod.1, type=3, digits=3, fallen.leaves=TRUE)
printcp(mod.1)
new.mod.1=prune.rpart(mod.1, 0.010000)
rpart.plot(new.mod.1, type=3, digits=3, fallen.leaves=TRUE)
1-sum(residuals(new.mod.1)^2)/sum((rgdata$log.rgr-mean(rgdata$log.rgr))^2)
R2 = 0.07

mod.2=rpart(log.rgr~log.lma+log.rmf+Comp.1, data=rgdata)
rpart.plot(mod.2, type=3, digits=3, fallen.leaves=TRUE)
printcp(mod.2)
new.mod.2=prune.rpart(mod.2, 0.010000)
rpart.plot(new.mod.2, type=3, digits=3, fallen.leaves=TRUE)
1-sum(residuals(new.mod.2)^2)/sum((rgdata$log.rgr-mean(rgdata$log.rgr))^2)
R2 = 0.06

mod.3=rpart(log.rgr~log.lma+log.rmf+Comp.2, data=rgdata)
rpart.plot(mod.3, type=3, digits=3, fallen.leaves=TRUE)
printcp(mod.3)
new.mod.3=prune.rpart(mod.3, 0.010000)
rpart.plot(new.mod.3, type=3, digits=3, fallen.leaves=TRUE)
1-sum(residuals(new.mod.3)^2)/sum((rgdata$log.rgr-mean(rgdata$log.rgr))^2)
R2 = 0.05

mod.4=rpart(log.rgr~log.lma+log.mean.thick+Comp.1, data=rgdata)
rpart.plot(mod.4, type=3, digits=3, fallen.leaves=TRUE)
printcp(mod.4)
new.mod.4=prune.rpart(mod.4, 0.013902)
rpart.plot(new.mod.4, type=3, digits=3, fallen.leaves=TRUE)
1-sum(residuals(new.mod.4)^2)/sum((rgdata$log.rgr-mean(rgdata$log.rgr))^2)
R2 = 0.05

mod.5=rpart(log.rgr~log.mean.thick+log.ssl+Comp.2, data=rgdata)
rpart.plot(mod.5, type=3, digits=3, fallen.leaves=TRUE)
printcp(mod.5)
new.mod.5=prune.rpart(mod.5, 0.011868)
rpart.plot(new.mod.5, type=3, digits=3, fallen.leaves=TRUE)
1-sum(residuals(new.mod.5)^2)/sum((rgdata$log.rgr-mean(rgdata$log.rgr))^2)
R2 = 0.06

mod.6=rpart(log.rgr~log.mean.thick+log.lar1+Comp.2, data=rgdata)
rpart.plot(mod.6, type=3, digits=3, fallen.leaves=TRUE)
printcp(mod.6)
new.mod.6=prune.rpart(mod.6, 0.010000)
rpart.plot(new.mod.6, type=3, digits=3, fallen.leaves=TRUE)
1-sum(residuals(new.mod.6)^2)/sum((rgdata$log.rgr-mean(rgdata$log.rgr))^2)
R2 = 0.13


mod.7=rpart(log.rgr~log.mean.thick+log.smf+Comp.1, data=rgdata)
rpart.plot(mod.7, type=3, digits=3, fallen.leaves=TRUE)
printcp(mod.7)
new.mod.7=prune.rpart(mod.7, 0.010000)
rpart.plot(new.mod.7, type=3, digits=3, fallen.leaves=TRUE)
1-sum(residuals(new.mod.7)^2)/sum((rgdata$log.rgr-mean(rgdata$log.rgr))^2)
R2 = 0.07


mod.8=rpart(log.rgr~log.mean.thick+log.rmf+Comp.1, data=rgdata)
rpart.plot(mod.8, type=3, digits=3, fallen.leaves=TRUE)
printcp(mod.8)
new.mod.8=prune.rpart(mod.8, 0.010000)
rpart.plot(new.mod.8, type=3, digits=3, fallen.leaves=TRUE)
1-sum(residuals(new.mod.8)^2)/sum((rgdata$log.rgr-mean(rgdata$log.rgr))^2)
R2 = 0.06

mod.9=rpart(log.rgr~log.lar1+log.lmf+Comp.2, data=rgdata)
rpart.plot(mod.9, type=3, digits=3, fallen.leaves=TRUE)
printcp(mod.9)
new.mod.9=prune.rpart(mod.9, 0.010000)
rpart.plot(new.mod.9, type=3, digits=3, fallen.leaves=TRUE)
1-sum(residuals(new.mod.9)^2)/sum((rgdata$log.rgr-mean(rgdata$log.rgr))^2)
R2 = 0.13

mod.10=rpart(log.rgr~log.ssl+log.smf+log.light, data=rgdata)
rpart.plot(mod.10, type=3, digits=3, fallen.leaves=TRUE)
printcp(mod.10)
new.mod.10=prune.rpart(mod.10, 0.010000)
rpart.plot(new.mod.10, type=3, digits=3, fallen.leaves=TRUE)
1-sum(residuals(new.mod.10)^2)/sum((rgdata$log.rgr-mean(rgdata$log.rgr))^2)
R2 = 0.08

mod.10=rpart(log.rgr~log.ssl+log.smf+log.light, data=rgdata)
rpart.plot(mod.10, type=3, digits=3, fallen.leaves=TRUE)
printcp(mod.10)
new.mod.10=prune.rpart(mod.10, 0.010000)
rpart.plot(new.mod.10, type=3, digits=3, fallen.leaves=TRUE)
1-sum(residuals(new.mod.10)^2)/sum((rgdata$log.rgr-mean(rgdata$log.rgr))^2)
R2 = 0.08