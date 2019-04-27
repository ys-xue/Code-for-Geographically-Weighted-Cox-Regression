library(dplyr)
library(plyr)
source("./genData.R")

set.seed(942)
samplesize <- sample(56:64, 64*1000, replace = TRUE)

set.seed(2022)
datalist <- purrr::map(1:64000, function(x) 
    genData(0.02, c(0.7, 0.5, -0.8), samplesize[x], 0.1))

df <- ldply(datalist, data.frame)

df$block <- rep(rep(1:64, 1000), samplesize)

repind <- c(cumsum(samplesize)[64*(1:1000)][1],
            diff(cumsum(samplesize)[64*(1:1000)]))

df$replicate <- rep(1:1000, repind)

save(df, file = "./nulldata.RData")