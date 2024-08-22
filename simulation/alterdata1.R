library(dplyr)
library(plyr)
library(netmeta)
load("../real_data/louiscounty.RData")
source("./genData.R")

n_reps <- 1000
n_counties <- 64
sample_size_upper <- 64
sample_size_lower <- 56

set.seed(2133)
samplesize <-
  sample(sample_size_lower:sample_size_upper,
         n_counties * n_reps,
         replace = TRUE)
betaind <- rep(1:n_counties, n_reps)

beta <- c(0.7, 0.5, -0.8)

truebetas <-
  matrix(rep(beta, n_counties), nrow = n_counties, byrow = TRUE)
truebetas <-
  truebetas + 0.1 * (louiscounty[[2]][, 1] + louiscounty[[2]][, 2]
                     - mean(louiscounty[[2]][, 1]) -
                       mean(louiscounty[[2]][, 2]))

set.seed(2110)

datalist <- purrr::map(1:n_counties * n_reps, function(x)
  genData(0.02, truebetas[betaind[x],], samplesize[x], 60, 0.1))

df <- ldply(datalist, data.frame)
df$block <- rep(rep(1:n_counties, n_reps), samplesize)
repind <- c(cumsum(samplesize)[n_counties * (1:n_reps)][1],
            diff(cumsum(samplesize)[n_counties * (1:n_reps)]))
df$replicate <- rep(1:n_reps, repind)

save(df, truebetas, file = "./alterdata1.RData")
