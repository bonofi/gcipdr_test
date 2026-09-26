## corr_funcs.R contains simple modules to compute various correlation indexes. Correlation between a categorical and numerical variable is allowed.
## Copyright (C) 2018 Federico Bonofiglio

## This file is part of gcipdr.

    ## gcipdr is free software: you can redistribute it and/or modify
    ## it under the terms of the GNU General Public License as published by
    ## the Free Software Foundation, either version 3 of the License, or
    ## (at your option) any later version.

    ## gcipdr is distributed in the hope that it will be useful,
    ## but WITHOUT ANY WARRANTY; without even the implied warranty of
    ## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    ## GNU General Public License for more details.

    ## You should have received a copy of the GNU General Public License
    ## along with gcipdr.  If not, see <https://www.gnu.org/licenses/>.



## various correlation functions.....

# rank joint expectation

rjm <- function(x,y, ties=c("min","average", "max"),
                order=1)
{
    md <- model.frame(~x+y)
    ties<- match.arg(ties)
    rx <- rank(md[,1], ties=ties)
    ry <- rank(md[,2], ties=ties)
    prods <- (rx*ry)^order 
    res <- sum(prods)
    res  # note: since rank values > 0 always, all joint exps are > 0 always too. Thus, big values can always be managed with logarithm.

} 


# spearman rank correlation with unbroken ties: unbroken ties shall be more precise about each element order. Average and "max" is also unbroken but makes the ranking index value higher.

spcor <- function(x,y, ties=c("min","average", "max") )
{
    md <- model.frame(~x+y)
    ties<- match.arg(ties)
    rx <- rank(md[,1], ties=ties)
    ry <- rank(md[,2], ties=ties)
    res <- cor(rx,ry)
    res
}


rankcor <- function(X, y = NULL,
                    ties=c("min","average", "max") )
{
    ties<- match.arg(ties)
    if (!is.null(y) & is.numeric(y) & is.numeric(X) )
        res <- spcor(X, y)
    else
    {
        p <- dim(X)[2]
        combos <- combn(p, 2)
        J <- p*(p-1)/2
        lowtri <- unlist(
            lapply(1:J, function(j)
            {
                i <- combos[1, j]
                k <- combos[2, j]
                spcor(X[,i], X[,k])
            }
            )
        )
        res <- diag(nrow=p)
        res[lower.tri(res)] <- lowtri
        res[upper.tri(res)] <- t(res)[upper.tri(res)]
    }
    return(res)
}


#

pearsoncor <- function(x,y)
{
    z <- model.frame(~x+y)
    res <- cor(z[ ,1], z[ ,2])
    return(res)
}
#

momcor <- function(X, y = NULL )
{
    if (!is.null(y) & is.numeric(y) & is.numeric(X) )
        res <- pearsoncor(X,y)
    else
    {
        p <- dim(X)[2]
        combos <- combn(p, 2)
        J <- p*(p-1)/2
        lowtri <- unlist(
            lapply(1:J, function(j)
            {
                i <- combos[1, j]
                k <- combos[2, j]
                pearsoncor(X[,i], X[,k])
            }
            )
        )
        res <- diag(nrow=p)
        res[lower.tri(res)] <- lowtri
        res[upper.tri(res)] <- t(res)[upper.tri(res)]
    }
    return(res)
}


#  sample normal correlation score (van der waerden): ONLY FOR CONTINUOUS VARIABLES !!! CAN U USE BETA DIST FOR BERNOULLI IN ORDER TO GET SCORE ??

# smoothen up binary covariate

smoothen.binary <- function(bx)
{
    N <- length(bx)
    if (is.binary(bx))
        bx <- bx + rtruncnorm(N, a=0, sd = 0.001)    # add zero-sum continuos noise

    return(bx)        

}

#
smoothen.array <- function(X)
{
    res <- apply(X, 2, function(x)
        smoothen.binary(x)
        )
    return(res)
}

w.rankZ <- function(x, y,
                    ties=c("min","average", "max", "random"),
                    smooth = FALSE)
{
    ties<- match.arg(ties)  # should always be min ..
    md <- model.frame(~x+y)
    N <- dim(md)[1]
    if (smooth)
        md <-  smoothen.array(md)   # smoothens non-continuos variables
    rx <- rank(md[,1], ties=ties)
    ry <- rank(md[,2], ties=ties)
    ux <- rx/(N+1)  # 0-1 transformation
    uy <- ry/(N+1)   # 0-1
    zx <- qnorm(ux)  # mapped in normal standard space
    zy <- qnorm(uy) # mapped in normal standard space
    res <- cor(zx, zy) # sample normal score correlation
    res

}
# transformation based on CL, maybe slower convergence but also covers discrete rvs.

CL.trans <- function(x,y)
{
    md <- model.frame(~x+y)
    N <- dim(md)[1]
    zx <- (md[ ,1] - mean(md[ ,1], na.rm = T))/sd(md[ ,1], na.rm = T)
    zy <- (md[ ,2] - mean(md[ ,2], na.rm = T))/sd(md[ ,2], na.rm = T)
    res <- cor(zx, zy) # sample normal score correlation
    res
}

# smooth argument is deprecated !!!
Zcor <- function(x, y,
                 transf = c("central.limit", "weigted.rank"),
                 ties=c("min","average", "max", "random"),
                 smooth = FALSE )
{
    res <- switch(transf,
                  central.limit = CL.trans(x,y),
                  weighted.rank = w.rankZ(x,y,ties, smooth)
                  )
    return(res)
}


#

waerdencor <- function(X, y = NULL,
                       transf = c("central.limit", "weigted.rank"),
                       ties=c("min","average", "max", "random"),
                       smooth = FALSE )
{
    transf <- match.arg(transf)
    ties<- match.arg(ties)
    if (!is.null(y) & is.numeric(y) & is.numeric(X) )
        res <- Zcor(X, y, transf, ties, smooth)
    else
    {
        p <- dim(X)[2]
        combos <- combn(p, 2)
        J <- p*(p-1)/2
        lowtri <- unlist(
            lapply(1:J, function(j)
            {
                i <- combos[1, j]
                k <- combos[2, j]
                Zcor(X[,i], X[,k], transf, ties, smooth)
            }
            )
        )
        res <- diag(nrow=p)
        res[lower.tri(res)] <- lowtri
        res[upper.tri(res)] <- t(res)[upper.tri(res)]
    }
    return(res)
}



# mixed correlations in SN    
mixed_corr <- function(dat){
    
    N <- dim(dat)[1]
    wherenot <- apply(dat, 2, is.binary)
    cont <- which(!wherenot)
    
    dat[, cont] <- apply(
        dat[, cont], 2, 
        \(x) {
            # normal-ranks
            r_cont <- rank(x)
            qnorm(r_cont / (N + 1))
        }
    )
    
    mixed_matrix <- psych::mixedCor(
        dat, 
        c=cont, 
        d=which(wherenot)
    )
    
    return(
        mixed_matrix$rho)
    
}





