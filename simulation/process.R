##' Summarize simulation results for gw_cox paper
##'
##' for each measure, compute over all replicates, and average across all
##' counties.
##'
##' @usage process(parMat, sdMat, trueBetas)
##'
##' @param parMat the matrix of estimated parameters for all replicates
##' @param sdMat the matrix of estimated sds for all replicates
##' @param trueBetas the true vector of underlying parameters for one covariate
##' across all counties
##' @import matrixStats
##' @export

process <- function(parMat, sdMat,  trueBetas){
    ## mean absolute bias
    MAB <- mean(abs(sweep(parMat, 1, trueBetas)))
    
    ## mean sd
    MSD <- mean(rowSds(parMat))
    
    ## mean MSE
    MSE <- mean(sweep(parMat, 1, trueBetas)^2)
    
    ## Mean coverage probability
    MCP <- mean((parMat + 1.96 * sdMat > trueBetas) *
                    (parMat - 1.96 * sdMat < trueBetas))
    
    return(list(MAB = MAB,
                MSD = MSD,
                MSE = MSE,
                MCP = MCP))
}
