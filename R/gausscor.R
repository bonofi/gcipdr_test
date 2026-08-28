## gausscorr.R contains modules needed to perform NORTA transformation.
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




## Gaussian copula miscellanea (NORTA transform): Netwon-Raphson objective functions, first derivatives, MC integrals, additional steps to ensure SPD conditoin (see section A.1 of accompanyiyng reference and its supplementary material).
## PROJECTION OF AN ARBITRARY MATRIX IN STANDARD NORMAL SPACE (NORTA method)



### algorithm to convert an arbitrary correlation matrix to its corresponding one in standard normal space ....


## conversion of Rx into standard normal space corresponding matrix, Rz
# FINISH: introduce stoch argument and integral option everywhere

First.attempt.Rx_Rz.conversion <- function(Rx, marginals,
                                           means, sds, stoch=FALSE,
                                           lows=c(-5,-5), ups=c(5,5),
                                           pNorm=NULL, K=1000, NI_tol = 1e-05,
                                           NI_maxEval = 20)
{
                                        # Rx : corr matrix of X; marginals list of marginal inverse CDFs (quantile yelders)
    p <- length(marginals)    

    if(is.null(pNorm))
        pNorm <-rep(TRUE, p)
          
    combos <- combn(p, 2)

    J <- p*(p-1)/2 # number of separate root problems, equals dim(combos)[2]
    
    Rut <- lapply(1:J, function(j)
    {
        row <- combos[1, j]
        col <- combos[2, j]
        rxj <- Rx[row, col]

        g1 <- marginals[[row]]
        g2 <- marginals[[col]]

        mxj <- c(means[row], means[col] )
        sxj <- c(sds[row], sds[col] )
        pnrm <- pNorm[c(row,col)]

        # ── Closures capturing all arguments for this pair ────────────────
        # Give newtrap_one_cpp a clean single-argument (rz) interface
        
        fdist_closure <- function(rz)
          Compute.double.expectation.INS(
            rz,
            rx       = rxj,
            Gx1      = g1,
            Gx2      = g2,
            meanx    = mxj,
            sdx      = sxj,
            stoch    = stoch,
            lowlims  = lows,
            uplims   = ups,
            pNorm    = pnrm,
            K        = K,
            NI_tol   = NI_tol,
            NI_maxEval = NI_maxEval
          )
        
        fprime_closure <- function(rz)
          Compute.double.expectation.prime.INS(
            rz,
            rx       = rxj,
            Gx1      = g1,
            Gx2      = g2,
            sdx      = sxj,
            stoch    = stoch,
            lowlims  = lows,
            uplims   = ups,
            pNorm    = pnrm,
            K        = K,
            NI_tol   = NI_tol,
            NI_maxEval = NI_maxEval
          )
        
        # ── NEW: Combined wrapper — returns both objective and derivative in one call ──
        fdist_and_fprime_wrapper <- function(rz)
        {
          list(
            objective  = fdist_closure(rz),
            derivative = fprime_closure(rz)
          )
        }
        
        safecheck_closure <- function(rz)
          check.rz.bounds2(
            rz,
            rx       = rxj,
            Gx1      = g1,
            Gx2      = g2,
            meanx    = mxj,
            sdx      = sxj,
            stoch    = stoch,
            lowlims  = lows,
            uplims   = ups,
            pNorm    = pnrm,
            K        = K,
            NI_tol   = NI_tol,
            NI_maxEval = NI_maxEval
          )
        
        # ── Call Rcpp Newton-Raphson ─────────
        result <- newtrap_one_cpp(
          fdist_and_fprime = fdist_and_fprime_wrapper,  # ← NEW: combined function
          safecheck_arg    = safecheck_closure,
          start            = rxj,
          tol              = 0.01,
          maxit            = 50
        )
        
        out <- result$value
        
        if(sign(out)!=sign(rxj))
            warning("sign of copula pair correlation different from entry value !!!", call.=F)
        out
    }
    )
    
   res <- make.square.matrix(unlist( Rut ), p )
    res.bool <- make.square.matrix( unlist(
        lapply(Rut, function(x)
        attr(x, "adjusted", T)
        )
    ),
    p,
    diag.val = T
    )
    return( list( "matrix" = res, "flag" = res.bool  )
           )

}

                                        # Rz. FAILURE WARNING. if Rz is not semi-positive definite (e.g. is indefinite) the successive NORTA method may (will) fail (Gosh 01). In such case Gosh proposes to find the closest matrix S wich is semi-definite positive. This is a minimization task f(S-Rz) relative to S, where f is a linear operator. Approaches to solve these systems are "semidefinite programs" (see wiki: https://en.wikipedia.org/wiki/Semidefinite_programming). There is a R package here : https://cran.r-project.org/web/packages/Rcsdp/Rcsdp.pdf.


Rx.to.Rz.conv <- function(Rx, marginals,
                          means, sds,
                          stoch=FALSE, lows=c(-5,-5),
                          ups=c(5,5), pNorm=NULL, K=1000,
                          NI_tol = 1e-05, NI_maxEval = 20)
{
    first.try <- First.attempt.Rx_Rz.conversion(Rx,
                                                marginals,
                                                means, sds,
                                                stoch, lows,
                                                ups, pNorm, K,
                                                NI_tol, NI_maxEval
                                                )
    if ( is.SPD.matrix(first.try$matrix) )
        return(first.try$matrix)
    else
        res <- make.matrix.SPD( first.try$matrix, first.try$flag)
    return(res)
}

##

kruskalconv <- function(X) 2*sin(pi/6*X)



# this function goes into file 'ipd_rec.R'

convertRx <- function(Rx, marginals=NULL,
                      means=NULL, sds=NULL, stoch=FALSE, 
                      lows=c(-5,-5), ups=c(5,5),
                      pNorm=NULL, corrtype=c("moment", "rank"),
                      K=1000, NI_tol = 1e-05, NI_maxEval = 20)
{
    corrtype <- match.arg(corrtype)

    mc.stoch <- stoch & corrtype=="moment"
    Rz <- switch(corrtype,
                 moment= Rx.to.Rz.conv(Rx,
                                       marginals,
                                       means, sds,
                                       mc.stoch,
                                       lows, ups,
                                       pNorm, K,
                                       NI_tol, NI_maxEval
                                       ),                                        # also use mom in case rank need numerical search ...**
                 rank= kruskalconv(Rx)
                 )
    Rz

}


# COMPUTATION OF DOUBLE EXPECTATION 

## montecarlo integration approach


mccovx1x2 <- function(rz, Gx1, Gx2,..., rx= 0,
                      meanx=c(0,0), sdx=c(1,1),
                      pNorm=rep(TRUE,2), K=1000 )
{

    Sigmaz <- diag(2); Sigmaz[1,2] <- Sigmaz[2,1] <- rz
    z <- rmvnorm(K, sigma=Sigmaz)   # uniform sampling make little sense here ...
    p1 <- ifelse(rep(pNorm[1], K), pnorm(z[ ,1]), z[ ,1] )
    p2 <- ifelse(rep(pNorm[2], K), pnorm(z[ ,2]), z[ ,2] )
    x1 <- Gx1( p1 )
    x2 <- Gx2( p2 )
    covx1x2 <- sum(x1*x2)/K     # integral
    res <- (  covx1x2 -  prod(meanx)  ) / prod(sdx)  

    return(res - rx)
}

## left term, h, of first derivative of dmvnorm, with respect of rz. The original factorization is h*bnorm (see below)   

hLeftX <- function(z, rz)
{
    g <- 2*pi*sqrt(1-rz^2)
    gprime <- 2*pi*0.5*(-2*rz)*(1-rz^2)^(-0.5)
    l <- (-1)*( z[1]^2 - 2*rz*z[1]*z[2] + z[2]^2 )
    h <- 2*(1-rz^2) 
    f <- exp( l / h )   
    hsquare <- h^2
    hprime <- -4*rz
    lprime <- 2*z[1]*z[2]
    numint <- lprime*h - hprime*l
    fprime <- f*( numint/hsquare )
    numout <- ( numint/hsquare )*g - gprime
    derivfactor <- (numout/g)
    derivfactor
}


## mccovx1x2 first derivative wrt rz

mccovx1x2prime <- function(rz, Gx1, Gx2,..., sdx=c(1,1),  
                           pNorm=rep(TRUE,2), K=1000)
{
    Sigmaz <- diag(2); Sigmaz[1,2] <- Sigmaz[2,1] <- rz
    z <- rmvnorm(K,sigma=Sigmaz)
    p1 <- ifelse(rep(pNorm[1], K), pnorm(z[ ,1]), z[ ,1] )
    p2 <- ifelse(rep(pNorm[2], K), pnorm(z[ ,2]), z[ ,2] )

                                        # X-marginals ...
    x1 <- Gx1(p1)
    x2 <- Gx2(p2)
    hx <- apply(z, 1, function(x) hLeftX(x, rz) )
    covprime <- sum(x1*x2*hx)/K    # integral
    const <- 1/prod(sdx)
           
    return( const*covprime )
}



                                        # rewrite bivariate density without matrix notation
# TODO: see vignette("cubature"), example 'Double Gaussian', for possible performance improvements via vectorization ? 

bndens <- function(z, rz)
{
    const <- 1/(2*pi*sqrt(1-rz^2))
    err <- z[1]^2 - 2*rz*z[1]*z[2] + z[2]^2
    point <- const * exp( (-1)*err / (2*(1-rz^2)) )
    point
}



## factorized derivative with bivariate density as term !!

bnprime <- function(z, rz)
{
    g <- 2*pi*sqrt(1-rz^2)

                                        #     gsquare <- g^2  # this is not actually needed because it cancels out

    gprime <- 2*pi*0.5*(-2*rz)*(1-rz^2)^(-0.5)

    l <- (-1)*( z[1]^2 - 2*rz*z[1]*z[2] + z[2]^2 )

    h <- 2*(1-rz^2) 

    f <- exp( l / h )

    hsquare <- h^2

    hprime <- -4*rz

    lprime <- 2*z[1]*z[2]

    numint <- lprime*h - hprime*l

                                        #    fprime <- f*( numint/hsquare )  #not actually used term

    numout <- ( numint/hsquare )*g - gprime 

    deriv <- (numout/g)*bndens(z, rz)
    deriv

}


### distance function to minimize


### deterministic (numerical) double expectation

covx1x2 <- function(z, Gx1, Gx2,..., rz, pNorm)
{
          
    p1 <- ifelse(pNorm[1],pnorm(z[1]), z[1] )
    p2 <- ifelse(pNorm[2],pnorm(z[2]), z[2] )
    
    x1 <-  Gx1(p1, ...) 
    x2 <-  Gx2(p2, ...)

    cov <- x1*x2*bndens(z, rz)

    cov


}

#

covx1x2prime <- function(z, Gx1, Gx2,..., rz, pNorm)
{
    p1 <- ifelse(pNorm[1],pnorm(z[1]), z[1] )
    p2 <- ifelse(pNorm[2],pnorm(z[2]), z[2] )

    x1 <-  Gx1(p1, ...) 
    x2 <-  Gx2(p2, ...)

    deriv <- x1*x2*bnprime(z, rz)
    deriv
}



#  distance of covx1x2 from fiven rx


covx1x2dist <- function(rz, rx, Gx1, Gx2,...,
                        meanx=c(0,0), sdx=c(1,1),
                        lowlims= c(-5,-5), uplims=c(5,5),
                        pNorm=rep(TRUE,2), NI_tol = 1e-05,
                        NI_maxEval = 20)
{
                                        # K is defunct here only needed to to avoind arg conflicts with mc integration
    if ( abs(rz)==1 )
        rz <- sign(rz)*0.99  # at rz=pm 1 integral is degenerate and loopspins
    res <- adaptIntegrate(covx1x2,
                          lowlims,
                          uplims,
                          tol = NI_tol,
                          maxEval = NI_maxEval,
                          Gx1=Gx1,
                          Gx2=Gx2,
                          rz=rz,
                          pNorm=pNorm
                          )
    cor <- ( res$integral - prod(meanx) )/prod(sdx)

    cor - rx

}

# first derivative of dist func ...

covx1x2distprime <- function(rz, rx, Gx1, Gx2,...,
                             sdx=c(1,1),
                             lowlims= c(-5,-5),
                             uplims=c(5,5),
                             pNorm=rep(TRUE,2),
                             NI_tol = 1e-05, NI_maxEval = 20)
{
    if(abs(rz)==1)
        rz <- sign(rz)*0.99  # at rz=pm 1 integral is degenerate and loopspins
    res <- adaptIntegrate(covx1x2prime, lowlims,
                          uplims, tol = NI_tol,
                          maxEval = NI_maxEval, Gx1=Gx1,
                          Gx2=Gx2, rz=rz, pNorm=pNorm
                          )
    
    res$integral*(1/prod(sdx))

}



###  compute double expectation projected in normal space (INS), and norm between observed/projected correlation, numerically or stochastically .... ACTUALLY this is correlation no tdouble expectatoin ...

Compute.double.expectation.INS <- function( rz, rx, Gx1, Gx2,...,
                                           meanx=c(0,0), sdx=c(1,1),
                                           stoch = FALSE,
                                           lowlims= c(-5,-5),
                                           uplims=c(5,5),
                                           pNorm=rep(TRUE,2), K=NULL,
                                           NI_tol = 1e-05, NI_maxEval = 20 )
{
    if ( stoch )
        out <- mccovx1x2(rz, Gx1, Gx2,
                         rx = rx, meanx = meanx,
                         sdx = sdx, pNorm = pNorm,
                         K = K
                         )
    else
        out <- covx1x2dist(rz, rx,
                           Gx1, Gx2,
                           meanx = meanx,
                           sdx = sdx,
                           lowlims = lowlims,
                           uplims = uplims,
                           pNorm = pNorm,
                           NI_tol = NI_tol,
                           NI_maxEval = NI_maxEval
                           )
    return(out)
}

#

Compute.double.expectation.prime.INS <- function(rz, rx,
                                                 Gx1, Gx2,...,
                                                 sdx=c(1,1), stoch = FALSE,
                                                 lowlims= c(-5,-5),
                                                 uplims=c(5,5),
                                                 pNorm=rep(TRUE,2), K=NULL,
                                                 NI_tol = 1e-05, NI_maxEval = 20 )
{
    if ( stoch )
        out <- mccovx1x2prime(rz,
                              Gx1, Gx2,
                              rx = rx,
                              sdx = sdx,
                              pNorm = pNorm,
                              K = K
                              )
    else
        out <- covx1x2distprime(rz, rx,
                                Gx1, Gx2,
                                sdx = sdx,
                                lowlims = lowlims,
                                uplims = uplims,
                                pNorm = pNorm,
                                NI_tol = NI_tol,
                                NI_maxEval = NI_maxEval
                                )
    return(out)
}



# safecheck function in cor search ...

check.rz.bounds <- function(rz, ...)
{
    modified <- FALSE    
    if ( abs(rz) >= 1 )
    {
        rz <- sign(rz)*0.99
        modified <- TRUE
    }
    attributes(rz)$modified <- modified
    attributes(rz)$break.loop <- FALSE  # always false here 
    return(rz)
}

# AGAIN: problem of spilling search likely arises when, for given marginals, correlation in standard space is already close to the extremal. In that case, we may observe spilling, or singularities.

check.rz.bounds2 <- function(rz, rx,
                             Gx1, Gx2,...,
                             meanx, sdx,
                             stoch,
                             lowlims, uplims, pNorm, K, NI_tol, NI_maxEval )
{
                                        # if corr spill over natural boundary, space boundary is likely a cul de sac already. E.g. derivative may already quickly converge to 0 for 0.99, sending newtrap into a loop. By safechecking, better offer a small array of options on boundary, say 0.7 to 0.99 by 10, and pick the boundary value yielding minimal distance (this should also protect from the derivative freezing on a stationary point)

    modified <- FALSE
    break.loop <- FALSE
    
    if ( abs(rz) >= 1 )
    {
        if (abs(rz) == 1)
            rz <- sign(rz)*0.99  # BEWARE!! TODO(me): to use sign(rx) is much safer !!! think to change
        else
        {
            rz.grid <- sign(rz)*seq(0.735, 0.99, length = 10)  # BEWARE sign
            inits <- sapply(rz.grid, function(x)
                Compute.double.expectation.INS(x, rx,
                                               Gx1, Gx2,
                                               meanx = meanx,
                                               sdx = sdx,
                                               stoch = stoch,
                                               lowlims = lowlims,
                                               uplims = uplims,
                                               pNorm = pNorm, K = K,
                                               NI_tol = NI_tol, NI_maxEval = NI_maxEval
                                               )
                )
            init <- min( abs(inits), na.rm = T)
            rz <- rz.grid[abs(inits) == init]
        }
        modified <- TRUE  # value has been modified ....  ifelse( rz == 0.99, 0.9, rz - ( rz - round(rz,2)) )
                                        # if derivative already approaching saddle point on unity bordet (pm 0.9) order breaking loop and return border value with smallest norm ...

        flat.point <- Compute.double.expectation.prime.INS(
            sign(rz)*0.9,
            rx, Gx1, Gx2,
            sdx = sdx,
            stoch = stoch,
            lowlims = lowlims,
            uplims = uplims,
            pNorm = pNorm, K = K,
            NI_tol = NI_tol,
            NI_maxEval = NI_maxEval
        )
        if (  abs(flat.point) < 0.01 ) # if yes, flag out with breaking command for external loop in newton-raphson search (newtrap) ....
            break.loop <- TRUE
    }
    attributes(rz)$modified <- modified
    attributes(rz)$break.loop <- break.loop # bug: attributes must be always set out of loops or logical blocks (better before returning), otherwise may turn NULL
    
    return(rz)
}

#

is.square.matrix <- function( mat )
{
    if (!is.matrix(mat))
        stop("argument must be a matrix")

    p <- dim(mat)[1]
    k <- dim(mat)[2]

    return( p == k )

}

#

make.square.matrix <- function( vec, p,
                               from.triangular = T,
                               diag.val = F )
{
    if ( from.triangular )
    {
        res <- diag(nrow=p, ncol=p) - diag.val  # if logical matrix diagonal must be zero
        res[lower.tri(res)] <- vec
        res[upper.tri(res)] <- t(res)[upper.tri(res)]
    }
    else
        res <- matrix( vec, nrow = p)
    return(res)
}



is.SPD.matrix <- function(mat)
    all( eigen(mat, T,T)$values > 0 )  

# make matrix SPD if flagged: if matrix is flagged some matrix items are artificially set to pm0.99. Likely there exist no such given pearson correlation for simulated variables (problem likely to arise if one variable is binary and other continuous). Then, likely solution is to let correlation diverge from extremal. In the given settings we can at most reduce the flagged item of an amount 0.005*50 = 0.25. Change these values if u thing a greater discount is needed.

make.matrix.SPD.if.flagged <- function( num.vec, bool.vec, p )
{
    mat <- make.square.matrix( num.vec, p, FALSE)
    i <- 0
    while( !is.SPD.matrix(mat) )
    {
        num.vec[bool.vec] <- num.vec[bool.vec] - sapply(num.vec[bool.vec], function(i) sign(i)*0.005 ) # tweak flagged items an amount opposite to their sign
        mat <- make.square.matrix( num.vec, p, FALSE)
        i <- i+1
        if ( i > 200 )
        {
                                        # 50  # increase iteration treshold !!!
            warning("no SPD matrix found")
            break
        }
    }

    attributes(mat)$SPD.adjustment <- "targeted"  # only if flagged terms are targeted ...
    return(mat)
}

                                        # tweak sequentially to make matrix SPD: first tweak higher corr values which may be the ones on boundary space. Tweak sequentially by decreasing rank ...

tweak.to.SPD.sequentially <- function( num.vec, p )
{
                                        # num.vec is lower triangular vector here
    ref.sign <- sign(num.vec)
    inv.rank <- match(rev(sort(abs(num.vec))), abs(num.vec))
    n <- length(inv.rank)
    mat <- make.square.matrix( num.vec, p )
    j <- 1  # inner loop
    i <- 0  # outer loop
    while( !is.SPD.matrix(mat) )
    {
        num.vec[inv.rank[j]] <- num.vec[inv.rank[j]] - sign( num.vec[inv.rank[j]] )*runif(1, max = 0.01)
        mat <- make.square.matrix( num.vec, p )  # update
        j <- j +1 # descend along matrix elements rank
        if ( j > n )
        {
            j <- 1  # reset j if all inv.rank elements have been exauhsted
        }
        i <- i+1

        if ( i > 200 )
        {
            warning("no SPD matrix found")
            break
        }
    }

    if (any(sign(num.vec) != ref.sign))
        warning("make.matrix.SPD: some sign has changed")

    attributes(mat)$SPD.adjustment <- "sequential"
                      
    return(mat)

}

# make matrix SPD in general case: For some reason matrix may not be SPD even if conversion succeded. In these case we need some sort of general tweaking in the hope matrix will turn SPD. The principle is again that some matrix items may lie on the boundary of existence for the give nmarginals. Since the identity matrix is always SPD, by letting converging values to zero we hope to make matrix SPD. Thus, we discount each items by randomly reducing its value by a certain amount. U may have to put in place some safacheck if for instanc ethe procedure end up inverting the sign (not allowed).

tweak.to.SPD.not.sequentially <- function( num.vec, p )
{
                                        # num.vec is lower triangular vector here
    ref.sign <- sign(num.vec)
    mat <- make.square.matrix( num.vec, p )
    i <- 0
    while( !is.SPD.matrix(mat) )
    {
        num.vec <- num.vec - sapply(num.vec, function(i) sign(i)*runif(1, max = 0.01) ) # tweak flagged items an amount opposite to their sign
        mat <- make.square.matrix( num.vec, p )
        i <- i+1
        if ( i > 50 )
        {
            warning("no SPD matrix found")
            break
        }
    }

    if (any(sign(num.vec) != ref.sign))
        warning("make.matrix.SPD: some sign has changed")
    attributes(mat)$SPD.adjustment <- "not.sequential"
    return(mat)
}

#

make.matrix.SPD.if.NOT.flagged <- function( num.vec, p )
{
    mat <- tweak.to.SPD.sequentially( num.vec, p)
    if ( !is.SPD.matrix(mat) )
        mat <- tweak.to.SPD.not.sequentially( num.vec, p)

    if ( !is.SPD.matrix(mat) )
        warning("sequential/not-sequential SPD search: no SPD matrix found")
    return(mat)
}




make.matrix.SPD <- function( mat, flag )
{

    if (!is.square.matrix(mat) | !is.square.matrix(flag))
        stop( "arguments must be square matrices" )

    warning("make.matrix.SPD was invoked: correlation matrix conversion failed. Suggestions: if you used numerical integration try decreasing 'NI_tol' (just a little) and/or increasing 'NI_maxEval'. If you used stochastic integration, try increasing SI_k")
    p <- dim(mat)[1]
    vec <- as.numeric(mat)
    bool <- as.logical(flag)
    if ( any(bool) )
    {
                                        # if any matrix element is flagged ....
        if ( is.SPD.matrix(mat) )
            return(mat)
        else
            res <- make.matrix.SPD.if.flagged( vec, bool, p)
    }
    else
    {

        if ( is.SPD.matrix(mat) )  # TODO(me) : is missing the case where mat is flaged but not flagged elements need also be tweaked ...
            return(mat)
        else
            res <- make.matrix.SPD.if.NOT.flagged( mat[lower.tri(mat)], p )  # feed lower trinagular only
    }
    return(res)
}






################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################

#### redundant functions DEPRECATED....

## first derivative of standard bivariate normal (active)

bndensprime <- function(z, rz){  ## redundant function !!!

 g <- 2*pi*sqrt(1-rz^2)
    
     gsquare <- g^2

 gprime <- 2*pi*0.5*(-2*rz)*(1-rz^2)^(-0.5)
    
 l <- (-1)*( z[1]^2 - 2*rz*z[1]*z[2] + z[2]^2 )
    
        h <- 2*(1-rz^2) 

    f <- exp( l / h )

    hsquare <- h^2

    hprime <- -4*rz

    lprime <- 2*z[1]*z[2]

     numint <- lprime*h - hprime*l
    
    fprime <- f*( numint/hsquare )

    numout <- fprime*g - gprime*f

    deriv <- numout/gsquare

    deriv
}




# +
# +
# +
# +
# +
# +
# +
# +
# +
# +
# +
# +



# +
# +
# +
# +
# +
# +
# +
# +
# +
# +
# +
# +

## guess init (DEPRECATED)

guess.rz.init <- function(rxj, g1, g2,
                          mxj, sxj, pnrm,
                          lowlims = c(-5, -5), uplims = c(5, 5),
                          stoch=FALSE, K = 1000, NI_tol = 1e-05, NI_maxEval = 20 )
{
    grid <- seq( sign(rxj)*0.01, sign(rxj)*0.99, length = 100 )
    inits <- sapply( grid, function(x)
        Compute.double.expectation.INS(x, rx = rxj,
                                       Gx1 = g1, Gx2 = g2,
                                       meanx = mxj, sdx = sxj,
                                       stoch = stoch,
                                       lowlims = lowlims,
                                       uplims = uplims, pNorm = pnrm,
                                       K = K, NI_tol = NI_tol, NI_maxEval = NI_maxEval
                                       )
        )
    init <- min( abs(inits), na.rm = T)

    out <- grid[abs(inits) == init]
    
    return(out)
}




##############  obsolete code ( deprecated, not deactivated )

 distrXrZ <- function(rz, X1, X2, Zgrid, rx, dz=0.05, code=c("C++","R")){

     code <- match.arg(code)
  
    dist <- jointE(X1, X2, Zgrid, dz, code=code) - rz 

     dist

    }



### u can try to directly minimize the above function with oprim(), or optimize()


## in the case the first derivative of distfunc is (first derivative of binormal standanrd)...

 distrXrZprime <- function(rz, X1, X2, Zgrid, dz=0.05, code=c("C++","R")){

     code <- match.arg(code)

 BNder <- switch(code,

       "C++"= primebiGaussC(Zgrid, rz),

       "R"=standBivNprime(Zgrid, rz)

                 )

  deriv <- sum( X1*X2*BNder*dz^2 )
     
     deriv

  }


## joint expectation Cairo, Nelson 97 Equation 1 . jointE DEPRECATED


jointE <- function(X1, X2, Zgrid, rz, dz, code=c("C++","R")){

  code <- match.arg(code)
        
  ZbivarNorm <- switch(code,
                       
       "C++"= densbiGaussC(Zgrid, rz) ,
                       "R"={
    Sigmaz <- diag(2)
      Sigmaz[1,2] <- Sigmaz[2,1] <- rz
                  dmvnorm(Zgrid, sigma=Sigmaz) 
               }
    )

                         
    COVx1x2 <- sum( X1*X2*ZbivarNorm*dz^2 )  # double integral 

    COVx1x2
    
    }


##


#

standBivNorm <- function(z, rz){ # z is n x 2 grid

   dens <- apply(z, 1, function(x) bndens(x, rz ) )  ### this step possibly very slow ..code in C ++

   dens

   } 

#
standBivNprime <- function(z, rz){

 deriv <- apply(z, 1, function(x) bnprime(x, rz ) )  ### this step possibly very slow ..code in C ++

   deriv

    }



 
