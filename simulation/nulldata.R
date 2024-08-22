library(dplyr)
library(plyr)
source("./genData.R")

n_reps <- 1000
n_counties <- 64
sample_size_upper <- 64
sample_size_lower <- 56

set.seed(942)
samplesize <-
  sample(sample_size_lower:sample_size_upper,
         n_counties * n_reps,
         replace = TRUE)

set.seed(2022)

beta <- c(0.7, 0.5, -0.8)

datalist <- purrr::map(1:n_counties * n_reps, function(x)
  genData(0.02, beta, samplesize[x], 60, 0.1))

df <- ldply(datalist, data.frame)

df$block <- rep(rep(1:n_counties, n_reps), samplesize)

repind <- c(cumsum(samplesize)[n_counties * (1:n_reps)][1],
            diff(cumsum(samplesize)[n_counties * (1:n_reps)]))

df$replicate <- rep(1:n_reps, repind)

save(df, file = "./nulldata.RData")