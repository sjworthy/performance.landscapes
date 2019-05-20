## Input for R

install.packages("dplyr")
install.packages("rjags")
install.packages("runjags")
install.packages("parallel")

library(dplyr)
library(rjags)
library(runjags)
library(parallel)

setwd("")

final.data = read.csv("Final.Data.csv", header = T, row.names = 1)

data <- list(nf = nrow(final.data), ns = 122, np = 215,
             pg = final.data$log.rgr, light = final.data$log.light, sp = final.data$sp.number,
            plot = final.data$Plot.2, trait1 = final.data$scale.log.lma, trait2 = final.data$log.rmf, ini.size = final.data$size)
             
n.adapt=10000
n.update=10000
n.iter=50000

mod1 <- jags.model("rgr.JAGS.r", data=data, n.chains=6, n.adapt=n.adapt)
update(mod1, n.update)
results <- coda.samples(mod1, c("mubeta", "pval","tau", "alpha", "beta2", "beta3", "beta4", "beta5", "beta6", "beta7", "beta8", "beta9", "beta10"), n.iter=n.iter)
results_dic <- dic.samples(mod1, n.iter, "pD") 

results
saveRDS(results, "LMA.RMF.Light.jags_results.RData")

results_dic
saveRDS(results_dic, "LMA.RMF.Light.jags_results_dic.RData")

summary(results)
gelman.diag(results)
