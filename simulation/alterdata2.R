library(dplyr)
library(plyr)
library(netmeta)
load("../real_data/louiscounty.RData")
load("../real_data/graphDist.RData")
source("./genData.R")


set.seed(2007)
samplesize <- sample(56:64, 64 * 1000, replace = TRUE)
betaind <- rep(1:64, 1000)



beta <- c(0.7, 0.5, -0.8)

truebetas <- matrix(rep(c(0.7, 0.5, -0.8), 64), nrow = 64, byrow = TRUE)
truebetas <- truebetas + 0.08 * (graphDist[45, ])

datalist <- purrr::map(1:64000, function(x) 
    genData(0.02, truebetas[betaind[x],], samplesize[x], 60, 0.1))

df <- ldply(datalist, data.frame)
df$block <- rep(rep(1:64, 1000), samplesize)
repind <- c(cumsum(samplesize)[64*(1:1000)][1], 
            diff(cumsum(samplesize)[64*(1:1000)]))
df$replicate <- rep(1:1000, repind)

save(df, truebetas, file = "./alterdata2.RData")
