library(survival)
library(parallel)
library(netmeta)
source("./process.R")
load("../real_data/louiscounty.Rdata")
load("./nulldata.RData")
load("../circdistance.RData")


nsims <- 1000
bw <- seq(0.5, 50, by = 0.5)

distMatbackup <- distMat
distMatbackup[distMatbackup == 1] <- 0

## null simulation, truebetas is a matrix with 64 rows, each row having the
## same value
truebetas <- matrix(rep(c(0.7, 0.5, -0.8), each = 64), nrow = 64)


tics <- matrix(0, nrow = length(bw), ncol = nsims)
p1locmat <- p2locmat <- p3locmat <-
    sd1locmat <- sd2locmat <- sd3locmat <-
    matrix(0, nrow = 64, ncol = nsims)

sim <- function(k){
    p1gwmat <- p2gwmat <- p3gwmat <-  
        sd1gwmat <- sd2gwmat <- sd3gwmat <- 
        matrix(0, nrow = 64, ncol = nsims)
    
    allloccoefs <- allcoefs <- matrix(0, nrow = 3 * nsims, ncol = 64)
    
    pb <- txtProgressBar(min=0, max = nsims, style=3)
    for (j in 1:nsims){
        repi <- subset(df, replicate == j)
        for (i in 1:64){
            distMat <- distMatbackup
            distMat <- exp(- distMat / k)
            wvec <- rep(distMat[i, ], times = as.numeric(table(repi$block)))
            gwfit <- coxph(Surv(survtime, status) ~ age + Black + Married, 
                           data = repi, weights = wvec)
            gwbeta <- gwfit$coefficients
            T <- subset(repi, block == i)$survtime
            X <- subset(repi, block == i, c("age", "Black", "Married"))
            delta <- subset(repi, block == i)$status
            X <- X[order(T), ]
            delta <- delta[order(T)]
            X <- as.matrix(X[dim(X)[1]:1, ])
            localfit <- coxph(Surv(survtime, status) ~ age + Black + Married, 
                              data = subset(repi, block == i),
                              init = gwbeta,
                              control = coxph.control(iter.max = 0))
            localJinv <- localfit$var
            localScore <- colSums(coxph.detail(localfit)$score)
            localK <- localScore %*% t(localScore)
            loglike <- localfit$loglik[1]
            penalty <- localJinv %*% localK
            tics[which(bw == k), j] <- tics[which(bw == k), j] - 2 * loglike +
                2 * sum(diag(gwfit$var %*% localK))
            
            if (k == bw[1]){
                locfit <- coxph(Surv(survtime, status) ~ age + Black + Married, 
                                data = subset(repi, block == i))
                loglike2 <- sum(rev(delta) * 
                                (X %*% locfit$coefficients -
                                log(cumsum(exp(X %*% locfit$coefficients)))))
                aiclocal[j] <- aiclocal[j] - 2 * loglike2
                
                p1locmat[i, j] <- locfit$coefficients[1]
                p2locmat[i, j] <- locfit$coefficients[2]
                p3locmat[i, j] <- locfit$coefficients[3]
                
                sd1locmat[i, j] <- sqrt(diag(locfit$var))[1]
                sd2locmat[i, j] <- sqrt(diag(locfit$var))[2]
                sd3locmat[i, j] <- sqrt(diag(locfit$var))[3]
            }
            p1gwmat[i, j] <- gwfit$coefficients[1]
            p2gwmat[i, j] <- gwfit$coefficients[2]
            p3gwmat[i, j] <- gwfit$coefficients[3]
            
            sd1gwmat[i, j] <- sqrt(diag(gwfit$var))[1]
            sd2gwmat[i, j] <- sqrt(diag(gwfit$var))[2]
            sd3gwmat[i,] <- sqrt(diag(gwfit$var))[3]
            
        }
        setTxtProgressBar(pb, j) 
    }
    a1 <- process(p1locmat, sd1locmat, trueBetas = truebetas[,1])
    a2 <- process(p2locmat, sd2locmat, trueBetas = truebetas[,2])
    a3 <- process(p3locmat, sd3locmat, trueBetas = truebetas[,3])
    
    a4 <- process(p1gwmat, sd1gwmat, trueBetas = truebetas[,1])
    a5 <- process(p2gwmat, sd2gwmat, trueBetas = truebetas[,2])
    a6 <- process(p3gwmat, sd3gwmat, trueBetas = truebetas[,3])
    
    
    result <- rbind(a4, a5, a6, a1, a2, a3)
    row.names(result) <- c("p1gw", "p2gw", "p3gw", "p1loc", "p2loc", "p3loc")
    return(list(result = result,
                tic = tics[which(bw == k),],
                p1gwmat = p1gwmat, p2gwmat = p2gwmat, p3gwmat = p3gwmat,
                p1locmat = p1locmat, p2locmat = p2locmat, p3locmat = p3locmat))
}
