## ipd_rec.R contains modules to reconstruct IPD from IPD summaries only.
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



## IPD reconstruction using NORTA.
## Comparisons with original IPD: IPD summary statistics are directly computed on design matrix of (in principle unavailable) IPD example. 
 ## IPD summaries can otherwise be directly piped to function 'DataRebuild' that is core module of IPD reconstruction.


 

################ LOWER LEVEL OF DATA MANIPULATION: DATA RECONSTRUCTION ################


##################################################################
########### LEVEL OF LOOPING : SINGLE DATA-FRAME (raw data source) ###########
##################################################################

## RAW DATA HOLDER: IPD layer (HIGHER HIERARCHY)

create.dm.formula <- function(names)
{
    K <- length(names) # number of variables
    dmfrm <- paste(names, "+") # design matrix formula
    dmfrm[K] <- substr(dmfrm[K], 1, nchar(names[K]) ) # drop last plus sign (under convention that last var is not outcome ..., but it may be if K == 2, ADJUST later)
    # dmformula <- as.formula(c("~", dmfrm) )  # fixing late deprecation
    dmformula <- as.formula(paste0("~", dmfrm))
    attributes(dmformula)$formula.blocks <- dmfrm
    return(dmformula)
}


#' @title Conversion of original IPD into suitable format for further processing.
#'
#' @description At the level of a raw data holder this function converts the original IPD into a design-matrix-like format onto which marginal moments and pairwise correlations, needed by `DataRebuild()`, are computed.
#'
#' @param data an IPD.
#' @param fill.missing logical: if TRUE, application of `Return.key.IPD.summaries()` on the output of this function will instruct to compute IPD summaries such that `DataRebuild()` also generates imputations of the missing IPD records. 
#' @details See Details of `DataRebuild()` for more info on the required IPD format. 
#'
#' @seealso [Return.key.IPD.summaries()] for IPD summary computation, [DataRebuild()]


# data layer : return IPD desing matrix

Return.IPD.design.matrix <- function(data, fill.missing = F)
{
    N <- dim(data)[1] # raw number of observation (with missin values)
    K <- dim(data)[2]
    if ( is.null( names(data)  ) )
        names(data) <- paste("V",1:K, sep="")
    xnms <- names(data) # use variables names for formula symbols !!
### DESIGN MATRIX (ideally the raw data holder is working on this object and report its syntheses)
    dmformula <- create.dm.formula(xnms)
                                        # model.frame (option to keep missing values here !!!)
    mframe <- model.frame(dmformula, data, na.action= ifelse( fill.missing, na.pass, na.omit  ) )
                                        # design matrix contains ALL variables (including outcome as third column)
    md <- model.matrix( dmformula, mframe)  # NA removed automatically. md[,3]=y
    attributes(md)$is.sample.size.complete <- fill.missing # comunicates complete sample size withouth NA observation deletion (useful for missing data imputation later !!!!!)
    return(md)    
    
}

#

## INTERFACE between the raw- and reduced- data holder

### Raw data holder computes key summaries from the IPD design matrix (it must follow column order convention !!!),....and transfers them to the reduced data holder !!! 


## compute different typed of correlation matrixes  

Return.correlation.matrix <- function(x,
                                      corrtype = c("rank.corr", "moment.corr", "normal.corr"),
                                      y = NULL  )
{
    corrtype <- match.arg(corrtype)     
    res <- switch( corrtype,
             "rank.corr" =    rankcor(x, y), # Spearman
             "moment.corr" =    momcor(x, y),   # Pearson
             "normal.corr" = waerdencor(x, y) # van der Waerden 
             )
    attributes(res)$corr.type <- corrtype
    return( res )
}

# function extracts key IPD summaries (next data hierarchy resolution)- Typically complete.N arg not used, since it is inherited from md in my application. Must also comunicate data variables names ...

 
## NOTE: cannot know in advantage what model will be used, thus must offer widest possible variable combination for model induced reduction

compute.model.induced.stats.for.complete.data <- function( outcome, covariates)
{  # no NAs
    covariate <- as.matrix(covariates)
    sy <- t(covariates)%*%outcome  # X'y term (numerator): (model induced) "partial" sufficient statistics
    sysq <- sum(outcome^2) # sum of squares for gaussian model (MUST include all vars, and select later the one needed by name ...)
    mirs <- rbind( sy, sysq ) # (m)odel (i)nduced (r)educed (s)tats 
    attributes(mirs)$model.suff.stat <- TRUE  # ensure ordering convention is applied to this object
    return(mirs)
} 

#

compute.model.induced.stats.for.incomplete.data <- function( outcome, covariates)
{ # with possible NAs
    covariate <- as.matrix(covariates)
    sy <- apply(covariates, 2, function(x) sum(outcome*x, na.rm = T) )
    sysq <- sum(outcome^2, na.rm = T) # sum of squares for gaussian model (MUST include all vars, and select later the one needed by name ...)
    mirs <- matrix(c( sy, sysq ),ncol=1) # (m)odel (i)nduced (r)educed (s)tats 
    attributes(mirs)$model.suff.stat <- TRUE  # ensure ordering convention is applied to this object
    return(mirs)

} 


Return.model.induced.statistics <- function(outcome, covariates, is.data.missing)
{
    if ( is.data.missing )
        out <- compute.model.induced.stats.for.incomplete.data(outcome, covariates)
    else
        out <- compute.model.induced.stats.for.complete.data(outcome, covariates)
    return(out)
}


# looped across variables

# out.col = possible outcome columns (by number)

Return.model.induced.statistics.looped <- function(md, is.data.missing,
                                                   out.col = 2:dim(md)[2])
{
    # in my application outcome-column may be only one between number 2 or 3 (in more general applications, any)
    K <- dim(md)[2]
    if ( any(out.col > K ) )
        out.col <- out.col[out.col <= K]  # in case data has only one/two colums
    out <- lapply(out.col,
                  function(i)
                      Return.model.induced.statistics(
                          md[ ,i],
                          md[ ,-i],
                          is.data.missing )
                  )
    names(out) <- colnames(md)[out.col]
    return(out)

}


#' @title Computation of input IPD summaries from the original IPD
#'
#' @description This function computes marginal moments up to fourth degree, all pairwise correlation of an available IPD and Johnson system parameters for all IPD non-binary marginals. The IPD is supposed to be in a specific format (see `Return.IPD.design.matrix()` and Details of `DataRebuild()`).
#'
#' @details Ideally the the IPD holder uses this function output as a template for the IPD summaries format to be disclosed and needed in `DataRebuild()` by a IPD summary holder. If IPD has missing data, the module can compute IPD summaries as if missing records were imputed. This must be decided before via argument fill.missing in `Return.IPD.design.matrix()`. 
#'
#'
#' @seealso [Return.IPD.design.matrix()] for currently requested IPD format, [DataRebuild()]


# un can either have na.omit(md) or na.pass(md) .... handle accordingly ..

Return.key.IPD.summaries <- function( md,
                                     corrtype=c("moment.corr", "rank.corr", "normal.corr"),
                                     only.fixed.stats = FALSE, na.erase = TRUE,
                                     compute.jsp = FALSE )
{ # set na.erase = F if NAs may insert bias
    corrtype <- match.arg(corrtype)
    is.data.missing <- attr(md, "is.sample.size.complete", T)  # complete = including NAs
    if (is.null(is.data.missing) & na.erase) # should not: if fact it can be if md = Xspace (function: compute expected fixedstat)
        is.data.missing <- any(is.na(md))
    else
    {
        if (is.null(is.data.missing))
            is.data.missing <- FALSE   # keep NAs as last option
    }
    n <- sum(md[,1]) # number of independent observations
    variable.names <- attr(md, "dimnames")[[2]][-1]  # drop intercept !! THESE NAMES ARE INHERETED from fixed.stat !
    
############## check storage mode (important for data simulation later)
    x.mode <- apply(md[,-1], 2, is.binary)  # boolean: is x_j binary or not ?

#### (presumely reported statistics) : *compression phase* ### SIGNAL
    mirs <- Return.model.induced.statistics.looped( md, is.data.missing, 3:4 ) # BEWARE: now mirs is a list !!!!!
    if (only.fixed.stats)
        return( mirs )  # escape
    if ( is.data.missing ) # *data marginals*: ancillary stats (u must include outcome for norta simulation, which is redundant with sy[1], but it makes it easier to handle with other marginals)
        sx <- apply(md[ ,-1], 2, sum, na.rm = T) # marginals for data imputation
    else
        sx <- diag(diag(nrow=n))%*%md[ ,-1] # marginals not suited for imputation
    mx <- as.vector(sx/n) # firs moment
    sdx <- apply(md[, -1], 2, sd, na.rm=TRUE) # second moment

### optional moments
    skx <- apply(md[, -1], 2, skewness, na.rm=TRUE) # third
    ktx <- apply(md[, -1], 2, kurtosis, na.rm=TRUE) # fourth
#####################################
    moms <- cbind(mx, sdx, skx, ktx) # collect sample moments for NonLinearEquationSystem used in Johnson method
    if (compute.jsp)
        jsp <- apply( cbind( moms, x.mode ), 1, function(x)
        {
            if ( !x[5] )
                FitJohnsonDistribution(x[1], x[2], x[3], x[4])
            else
                rep(NA, 6)
        }
        )
    else
        jsp <- NULL
    Rx <- Return.correlation.matrix( md[,-1], corrtype)  #  sample correlation matrix
    rownames(Rx) <- colnames(Rx) <- variable.names # adding names to corr pairs
    out <- list(sample.size = n ,
                is.sample.size.complete = is.data.missing , # if NAs are kept, sample size includes all records (is complete)
                variable.names = variable.names ,
                first.four.moments = moms ,  # as an array (each column speaks for a data variable marginals)
                correlation.matrix = Rx , # data correlation matrix
                model.sufficient.statistics = mirs, # model-induced data reduction, X'y term in GLM and Cox models
                johnson.parameters = jsp , # solution to Johnson system
                is.binary.variable = x.mode  # boolean (is variable binary or continuous ? )
                )
    class(out) <- "key.ipd.summary"
    return( out )
}



##################################################################
########### LEVEL OF LOOPING : SINGLE DATA-FRAME AND ITS *H* REPETITIONS (simulations) ###########
##################################################################


# LOWER HIERARCHY RESOLUTION.  
## REDUCED DATA HOLDER is handed over key IPD summaries, which must be now exploited to try to recover much of the original raw data information.
# recovery happens via a method of simulating the raw data from its summaries ....

### DATA SIMULATION METHODS (permutation or norta)

# incomplete correlation reproduction (permutation search)

# function samples gamma or bernoulli distributed variables
# H : number of independnet R^n samples; n = number of independent unit observations, ...

# LOOPING LEVEL here is across data.frame variables and variables repetitions (H) .... no data level looping (across several data-sets)

gamma.method <- function( H, n,  mx, sdx, x.mode )
{
    if (x.mode)
        out <- replicate(H, rbinom(n, 1, mx) ) # COMMENT: U CAN MAKE IT FIXED HERE ? mx = 1, 1-mx = 0 !!!! (CHECK method: hint: it may introduce a systematic pattern, thus more bias ....)
    else
        out <- replicate(H, rgmom(n, mx, sdx) )
    return( out )
}

# looped over number of variables. So many variables as lenght of mx (sdx, and x.mode)

gamma.method.looped <- function( H, n,  mx, sdx, x.mode )
{
    K <- length(mx) # number of (compressed) variables
    if ( K != length(sdx) | K != length(x.mode) )
        stop( "mx, sdx, x.mode must have all equal length" )
    out <- lapply( 1:K,
                  function(i)
                      gamma.method( H, n, mx[i], sdx[i], x.mode[i] )
                  )
    return( out )
}


# function samples johnson or bernoulli distributed variables
# COMMENT : sometimes Johnsonsampler returns negative values tohugh Johnsonfit yielded type SB = 3. In tihs case introduce adjusting argument for quintile sampler

# rJohn is a wrapper for yJohnsonDistribution, allowing for control of unwanted negative values if IPYPE = SB

rJohn <- function(z, jp1, jp2, jp3, jp4, jp5, SBjohn.correction = F)
{  # also works for quantile computation
    out <- yJohnsonDistribution( z, jp1, jp2, jp3, jp4, jp5 )
    if ( SBjohn.correction & jp1 == 3 & any(out < 0))
    {
        N <- length(out[out < 0])
        if ( length(out[out > 0]) > 0 )
        {
            minpos <- min(out[out > 0]) # minimum hook. Alternative: out >= 0 ?
            replace <- seq(0, minpos, length = N + 1)[1:N]  # untie
            out[out < 0] <- replace  # replace
        }
        else   # if all values are negative
            out <- abs(out)
    }
                                        #  out[out < 0] <- 0  # easier version TODO(me) : check above is indeed a better approach then here
    return(out)
}

##


johnson.method <- function( H, n, mx, jsp , x.mode, SBjohn.correction = F)
{
    if (x.mode)
        replicate(H, rbinom(n, 1, mx) )
    else
        do.call(
            "cbind",
            mclapply(1:H, function(x)
                rJohn( rnorm(n),
                      jsp[1],
                      jsp[2],
                      jsp[3],
                      jsp[4],
                      jsp[5],
                      SBjohn.correction
                      )
                )
        )
}
  

# looped over number of variables. So many variables as lenght of mx (sdx, and x.mode) 

johnson.method.looped <- function(H, n, mx, jsp,  
                                  x.mode, SBjohn.correction = F)
{
    K <- length(mx)
    if ( K != length(x.mode) )
        stop( "mx and x.mode must have all equal length" )
    out <- lapply( 1:K,
                  function(j)
                      johnson.method(H, n,
                                     mx[j],
                                     jsp[ ,j],
                                     x.mode[j],  # jsp is a matrix now
                                     SBjohn.correction )
                  )
    out

}


# write function returning an expected extremal correlation to use in rsearch2 ...

Compute.extremal.correlation <- function(x, y, corr, corrtype=c("rank.corr", "moment.corr") )
{
    corrtype <- match.arg( corrtype )
    x.sorted <- sort(x)
    y.sorted <- sort(y)
    if(corr < 0 )
        y.sorted <- rev(y.sorted)
    out <-  Return.correlation.matrix(x.sorted, corrtype, y.sorted)
    return(out)
}

 # it operates over H variable repetitions

Compute.expected.extremal.correlation <- function(ref, target,
                                                  corr, corrtype = c("rank.corr", "moment.corr") )
{
    corrtype <- match.arg( corrtype )    
    if ( all(dim(ref) != dim(target)) )
        stop("dimension of ref and target do not match")
    H <- dim(ref)[2]
    res <- unlist( lapply(1:H, function(i)
        Compute.extremal.correlation(ref[ ,i], target[ ,i], corr, corrtype)
        )
        )
    out <- list("mc.average" = mean(res, na.rm = T),
                "mc.median" = median(res, na.rm = T),
                "mc.extremal" = ifelse(corr < 0, min(res, na.rm = T), max(res, na.rm = T) )
                )
    return(out)

}


# functions loops pairwise Rüschendorf (random-bernoulli) permutation search over H repetitions
# .repeated.obj = a matrix with columns as simulations of the *same* variable
# note: sorting is not directly necessary in rsearch2, but only in final reordering phase (Disentagle cluster). Sorting is directly needed if gsearch is applied (not actually)

#
Ruesch.permutation.search.looped <- function(H, sorted.repeated.reference, # sorted in ascending order
                                             repeated.variable, reference.corr, 
                                             corrtype = c("rank.corr", "moment.corr"),
                                             compute.eec = F )
{
    corrtype <- match.arg( corrtype )
    mc.extremal.corr <- NULL
    if (compute.eec)
        mc.extremal.corr <- Compute.expected.extremal.correlation(sorted.repeated.reference, repeated.variable,
                                                                  reference.corr, corrtype)
    out <- do.call(
        "cbind",
        mclapply(1:H, function(i)
        {
            mat <- rsearch2( sorted.repeated.reference[,i],
                            repeated.variable[,i], reference.corr, mc.extremal.corr,
                            d=0.01, rf.type = corrtype)
            mat[order(mat[ ,1]), 2]
        }
        )
    )
    return( out ) # return *only* permutated variable, not reference
}


## function disentagles data cluster generated by Rüsch search: if there are K clusters of dim n*H, the function reassemble column-wise to yield H clusters of dim n*(K + 1) 

Disentagle.data.clusters <- function(H, kcluster, reference, variable.names)
{
    K <- length(kcluster)  # number of variables clusters; H variable repetitions
    res <- mclapply(1:H, function(i)
    {  # list of H data-sets
        out <- cbind(reference[ ,i],
                     do.call(
                         "cbind",
                         lapply(1:K, function(j)
                             kcluster[[j]][ ,i]
                             )
                     )
                     )
        colnames(out) <- variable.names
        return( out )
    }
    )
    return( res )
}


# function applies Rüsch search over a *list* of repeated variables ( add additional search methods maybe ...)
# variables.list = list of repeated variables, slot 1 contains reference.
# reference.corr.vec = vector with reference correlations.
# output: disentagled list of H data-sets

permutation.search <- function(H, variables.list, reference.corr.vec,
                               variable.names = NULL, 
                               corrtype = c("rank.corr", "moment.corr"),
                               compute.eec = F, quick.check = FALSE )
{
    corrtype <- match.arg( corrtype )
    K <- length(variables.list) - 1  # minus reference (first list slot)
    if ( K != length(reference.corr.vec) )
        stop( "not enough reference correlations given" )
    if (quick.check)
    {
        h.cluster.size <- unlist( lapply( variables.list, function(x) dim(x)[2] ) )
        if ( any( h.cluster.size != H ) )
            stop( "some clusters have not H repetitions" )
    }
    sorted.reap.ref <- apply( variables.list[[1]], 2, function(x) sort(x) ) # sort reference in ascending order
    target.variables <- variables.list[-1] # drop reference
    kcluster <- lapply( 1:K, function(j)
        Ruesch.permutation.search.looped(H, sorted.reap.ref,
                                         target.variables[[j]],
                                         reference.corr.vec[j], corrtype,
                                         compute.eec
                                         )
        )
    list.of.H.permutated.datasets <- Disentagle.data.clusters( H, kcluster, sorted.reap.ref, variable.names )
    return( list.of.H.permutated.datasets )
    
   }



# function applies a method of permutaiton search based on a stochastic mechanism (extremal correlation structures) conditonal on only few correlation matrix entries (typically first row). Thus the obtained correlation structure is incomplete. OUTPUT: a list of of lenght H of (random) datasets (the space of simiar data X*)


Generate.with.Incomplete.correlation <- function(H, n, correlation.matrix,
                                                 moments, johnson.parameters, x.mode,
                                                 variable.names, SBjohn.correction = F,
                                                 compute.eec = F, corrtype=c("rank.corr", "moment.corr"),
                                                 marg.model=c("gamma", "johnson") )
{
    corrtype <- match.arg(corrtype)
    marg.model <- match.arg(marg.model)
    if (corrtype != "rank.corr" & corrtype != "moment.corr")
        stop("Permutation ordering currently admits only moment or rank correlations: other correlation type detected.")
    mx <- moments[ ,1] # first moment
    sdx <- moments[ ,2] # second moment     
    K <- length(mx) # number of variables 
    if( K != length(variable.names))
        stop("variable labels do not match the number of variables")
    Xstar <- switch(marg.model,  # simulate marginal columns relative to chosen method
                    gamma = gamma.method.looped( H, n,  mx, sdx, x.mode ),
                    johnson = johnson.method.looped( H, n, mx, johnson.parameters, x.mode, SBjohn.correction )
                    )
    if (K > 1)  # if there is more than 1 variable apply permutation search (always the case under convention) ...
    { 
        if (is.null(correlation.matrix))
        {  # but exit without permutation search if correlation matrix is null  ## ALLOW FOR NULL CORRELATION MATRIX !!!!
            unordered <- Disentagle.data.clusters( H, Xstar[2], Xstar[[1]], variable.names )
            return(unordered)
        }
        detected.corrtype <- attr(correlation.matrix, "corr.type", T) # TODO(me) : must allow for NULL corr matrix if data has only one column !!!
        if (detected.corrtype != corrtype)
            stop("Correlation-matrix type and argument corrtype do not match")
        lrx <- correlation.matrix[lower.tri(correlation.matrix)] # extraxt lower triangular
        refstat <- lrx[1:K-1]   # reference stats: first row of correlation matrix
        res <- permutation.search( H, Xstar, refstat, variable.names, corrtype, compute.eec ) # reorder variables according to given correlation !!!!
    }
    else
    { # else return variable's simulations unordered bundle (wil never be used in practice under my convention) ALLOW FOR SINGLE VARIABLE vector BE SIMULATED ....
        res <- Xstar[[1]]
                                        #  colnames(res) <- variable.names # TODO(me): dimension issues with colnames here .... fix later
    }
    return( list( artificial.data = res, copula.params = NULL ) )
}


#### norta approach: complete correlation reproduction


# return inverse marginal distribution

gamma.quintiles <- function( mx, sdx, n, x.mode,
                            corrtype=c("rank.corr", "moment.corr", "normal.corr") )
{
    corrtype <- match.arg(corrtype)
    
    if (x.mode)
        switch(corrtype,       # BEWARE: replacement of bernoulli with beta could bias second moment
               rank.corr =  function(x) qbmom(x, mx, n) ,   # beta quantiles (better for rank.norta method)
               moment.corr = function(x) qbinom(x, 1, mx),
               normal.corr = function(x) qbinom(x, 1, mx)  )
    else
        function(x) qgmom(x, mx, sdx)
}

# looped over different dataset variables .....
gamma.quintiles.looped <- function( mx, sdx, n, x.mode,
                                   corrtype=c("rank.corr", "moment.corr", "normal.corr") )
{
    corrtype <- match.arg(corrtype)
    K <- length(mx)
    if ( K != length(sdx) | K != length(x.mode) )
        stop( "mx, sdx, x.mode must have all equal length" )
    inverse.functions <- lapply( 1:K, function(j)
        gamma.quintiles( mx[j], sdx[j], n,  x.mode[j], corrtype )
        )
    return( inverse.functions )
}

#

johnson.quintiles <- function( mx, sdx, jsp, n, x.mode, 
                              corrtype=c("rank.corr", "moment.corr", "normal.corr"),
                              SBjohn.correction = F )
{
    corrtype <- match.arg(corrtype)     
    if (x.mode)
        switch(corrtype,  # BEWARE: replacement of bernoulli with beta could bias second moment
               rank.corr =  function(x) qbmom(x, mx, n) ,   # beta quantiles (better for rank.norta method)
               moment.corr = function(x) qbinom(x, 1, mx),
               normal.corr = function(x) qbinom(x, 1, mx)  )
    else
        function(x) rJohn(x,
                          jsp[1],
                          jsp[2],
                          jsp[3],
                          jsp[4],
                          jsp[5],
                          SBjohn.correction
                          )
}

 # looped over different dataset variables .....
johnson.quintiles.looped <- function( mx, sdx, jsp, n, x.mode,     
                                     corrtype=c("rank.corr", "moment.corr", "normal.corr"),
                                     SBjohn.correction = F )
{
    corrtype <- match.arg(corrtype)
    K <- length(mx)
    if ( K != length(sdx) | K != length(x.mode) )
        stop( "mx, sdx, x.mode must have all equal length" )
    inverse.functions <- lapply( 1:K, function(j)
        johnson.quintiles(mx[j],
                          sdx[j],
                          jsp[ ,j],
                          n,
                          x.mode[j],
                          corrtype,
                          SBjohn.correction
                          )
        )
    return( inverse.functions )
}


# norta :  convert correlation matrix

NORTAconvert.correlation.matrix <- function(correlation.matrix , marginal.inverse.distributions,
                                            first.moments, second.moments, stochastic.integration,
                                            yes.normal.percentile,
                                            corr.type=c("rank.corr", "moment.corr", "normal.corr"),
                                            SI_k = 8000, NI_tol = 1e-05, NI_maxEval = 20)
{
    corr.type <- match.arg(corr.type)
    if (corr.type == "normal.corr")
        return(correlation.matrix)  # already in standard normal space
    else
    {
        K <- length(marginal.inverse.distributions)
        if (dim(correlation.matrix)[1] != K |
            length(first.moments) != K |
            length(second.moments) != K |
            length(yes.normal.percentile) != K
            )
            stop("arguments must have same length")
        if (K > 1)
                                        # solves NORTA problem (finds corresponding corr matrix in standard normal space, CAIRO & NELSON 97)
            correlation.in.standard.normal.space <- switch(corr.type,   # Rz
                                                           "moment.corr" =
                                                               {
                                                                   convertRx(correlation.matrix,
                                                                             marginal.inverse.distributions,
                                                                             first.moments, second.moments,
                                                                             stoch = stochastic.integration,
                                                                             pNorm = yes.normal.percentile,
                                                                             K = SI_k,
                                                                             NI_tol = NI_tol,
                                                                             NI_maxEval = NI_maxEval
                                                                             )  # TODO(me): set K as formal argument or tweak here ....
                                                               },
                                                           "rank.corr" =  convertRx( correlation.matrix , corrtype="rank")  )
        else
            correlation.in.standard.normal.space <- diag(nrow = K)  # Rz
        if ( any( eigen(correlation.in.standard.normal.space, T,T)$values < 0 )  )
            stop("solution to the NORTA problem is not definite semipositive ")
        return(correlation.in.standard.normal.space)
    }
     
}




#  resample from NORTA multivariate distribution (one data-frame)

rNORTA <- function( sample.size, L, yes.normal.percentile,
                   marginal.inverse.distributions, variable.names )
{
    K <- length(marginal.inverse.distributions)
    if (length(yes.normal.percentile) != K | length(variable.names) != K)
        stop("arguments 'marginal.' to 'variable.names': must have same length")
    W <- replicate(K, rnorm( sample.size ))
    Z <- L%*%t(W)  # correlated standard multinormals (cholesky transformation)
    Z <- t(Z)
    datamatrix <- do.call(
        "cbind",
        lapply(1:K, function(j)
        {
            if (yes.normal.percentile[j])  # Phi(z) or z
                PZ <- pnorm(Z[ ,j])  # normal percentiles (P = pnorm)
            else
                PZ <- Z[ ,j]    # or normal quantiles (P = identity)
            marginal.inverse.distributions[[j]]( PZ )  # G(Phi(z)) or G(z)
        }
        )
    )
    
    colnames(datamatrix) <- variable.names
    return(datamatrix)
}

# reiterate NORTA sampling H times (several data-frames)

rNORTA.looped <- function(simulation.size, sample.size, L, yes.normal.percentile,
                          marginal.inverse.distributions, variable.names )
{
    out <- mclapply(1:simulation.size, function(i)
        rNORTA( sample.size, L, yes.normal.percentile,
               marginal.inverse.distributions, variable.names)
        )
    return(out)

}


# norta method

norta.method <- function( simulation.size, sample.size, correlation.matrix,
                         first.moments, second.moments, stochastic.integration, yes.normal.percentile,
                         marginal.inverse.distributions, variable.names,
                         corr.type=c("rank.corr", "moment.corr", "normal.corr"),
                         SI_k = 8000, NI_tol = 1e-05, NI_maxEval = 20, input.sn.corr = NULL)
{
    corr.type <- match.arg(corr.type)
    if ( any( eigen(correlation.matrix, T,T)$values < 0 )  )   # empirically observed correlation matrix
        stop("input correlation matrix is not definite semipositive")
    K <- length( marginal.inverse.distributions ) # nmbr of inverses
    if ( K != length(first.moments) |
         K != length(second.moments) |
         K != length(yes.normal.percentile) |
         K != dim(correlation.matrix)[1] )
        stop( "marginal inverses, .moments, percentiles, and correlation must have same length" )
    if (is.null(input.sn.corr))
    {
        if (corr.type == "normal.corr")
            correlation.in.standard.normal.space <- correlation.matrix # already is ...
        else
            correlation.in.standard.normal.space <- NORTAconvert.correlation.matrix(correlation.matrix,
                                                                                    marginal.inverse.distributions,
                                                                                    first.moments, second.moments,
                                                                                    stochastic.integration,
                                                                                    yes.normal.percentile,
                                                                                    corr.type, SI_k, NI_tol, NI_maxEval
                                                                                    )
    }
    else
        correlation.in.standard.normal.space <- input.sn.corr # adding option to externally tweak copula paramenter
                                        #
    rownames(correlation.in.standard.normal.space) <- colnames(correlation.in.standard.normal.space) <- variable.names
    U <- chol(correlation.in.standard.normal.space)  # Cholesky decomposition
    L <- t(U)
                                        # data-sets list
    out <- rNORTA.looped(simulation.size, sample.size, L, yes.normal.percentile,
                         marginal.inverse.distributions, variable.names)
    return( list(copula.inverse = out,
                 copula.params = correlation.in.standard.normal.space
                 )
           ) # list of H datasets of dimension sample.size*K each
}


## function directly sample from an arbitrary multivariate random vector with given marginal distributions and given correlation structure, Rx. The method exploits inverse probability sampling to map a multivatriate standard normal vector with correlation Rz (unknown) to the desired arbitrary multivariatre random vector. To succed a root problem must be solved for Rz, such that a solution for Rz ensures the properties of Rx are transfered into the desired multivariate vector.         *OUTPUT*: a list of (random) datasets of lenght H (the space of simiar data X*)                               

Generate.with.Complete.correlation <- function(H, n, correlation.matrix, moments, user_defined_marginals,
                                               johnson.parameters, stochastic.integration, x.mode,
                                               variable.names, SBjohn.correction = F,
                                               corrtype=c("moment.corr", "rank.corr" , "normal.corr"),
                                               marg.model=c("gamma", "johnson", "user_defined"),
                                               SI_k = 8000, NI_tol = 1e-05, NI_maxEval = 20,
                                               input.sn.corr = NULL, rescale.smoothed.binary = FALSE)
{
    corrtype <- match.arg(corrtype)
    marg.model <- match.arg(marg.model)
    detected.corrtype <- attr(correlation.matrix, "corr.type", T)
    if (detected.corrtype != corrtype & is.null(input.sn.corr))
        stop("Correlation-matrix type and argument corrtype do not match")
    mx <- moments[ ,1] # first moment
    sdx <- moments[ ,2] # second moment
    K <- length(mx)
    if( K != length(variable.names))
        stop("variable labels do not match the number of variables")
    marginals <- switch(marg.model, # marginal inverse distributions
                        gamma= gamma.quintiles.looped( mx, sdx, n, x.mode, corrtype ),
                        johnson = johnson.quintiles.looped(mx,
                                                           sdx,
                                                           johnson.parameters,
                                                           n,
                                                           x.mode, corrtype, SBjohn.correction
                                                           ),
                        user_defined = user_defined_marginals     # must provide a list of inverse distrib funcs (quintiles) for each marginal variable
                        )
                                        # if x is binary and corrtype not rank : norta second moment is bernoulli
    Ind <- x.mode & corrtype == "moment.corr"
    nortasd <- sqrt( Ind*( mx*(1-mx) ) ) + (1-Ind)*sdx
                                        # set boolean: use normal percentiles or not
    if (marg.model == "johnson")
        norm.perc <- x.mode  # (only Johnson meth don't use normal percentiles)
    else
        norm.perc <- x.mode | (1-x.mode) # else always true 

    res <- norta.method(H, n, correlation.matrix, mx, nortasd,
                        stochastic.integration, norm.perc,
                        marginals, variable.names,
                        corrtype, SI_k, NI_tol, NI_maxEval, input.sn.corr
                        )
    if ( ( corrtype == "rank.corr" & any(x.mode) & rescale.smoothed.binary ) ) # convert smoothed binary to integer bits
        res$copula.inverse <- lapply(res$copula.inverse, function(x)
        {
            x[, x.mode] <- apply( x[, x.mode], 2, map2bit )  ## TODO: x.mode will include categorical (Dir -> mult)
            return(x)
        }
        )
    return( list( artificial.data = res$copula.inverse,
                 copula.params = res$copula.params
                 )
           )

} 



#' @title IPD reconstruction from IPD summaries only.
#'
#' @description `DataRebuild()` generates artificial data, that is stochastic copies of the original IPD, by taking empirical IPD distributional summaries as input data only.
#'
#' @param H integer number of independent pseudodata replications to be generated.
#' @param n integer number of independent IPD observations.
#' @param correlation.matrix pairwise IPD correlation matrix. See 'Return.key.IPD.summaries'.
#' @param moments numeric array of IPD marginal moments. See 'Return.key.IPD.summaries'.
#' @param x.mode logical vector: is IPD marginal variable binary (TRUE) or not ?
#' @param user_defined_marginals list of inverse distributions (quantile functions) moment-defined for each marginal variable
#' @param johnson.parameters array of Johnson parameters for each IPD marginal variable. Depends on CRAN archived 'JohnsonDistribution' package. If NULL it is computed on given 'moments'.
#' @param stochastic.integration logical: should Monte Carlo integration be used to resolve Gaussian copula inversion (NORTA transformation)? Default to FALSE, that is numerical integration relying on package 'cubature' is used first. 
#' @param data.rearrange method of IPD dependence reconstruction based on all pairwise IPD correlations (norta), or on first degree correlations only (incomplete).
#' @param corrtype what type of IPD correlation matrix are you feeding in ? Spearman (rank.corr), Pearson (moment.corr), or Waerden (normal.corr). see Deatails.
#'
#' @param marg.model either "gamma", "johnson" (continuous marginals) or "user_defined". If not user defined, all binary marginals are modeled via a Bernoulli distribution or a Beta distribution if Kruskal analytic conversion is used (see below).
#'
#' @param variable.names names of IPD marginal variables. If NULL (Default) automatic labels are generated.
#'
#' @param SBjohn.correction logical. Should be Johnson marginal values corrected ? Default to FALSE. If TRUE, wrongly sampled negative values are set to the minimum positive sampled value
#'
#' @param compute.eec currently deprecated. Do not edit default value.
#'
#' @param checkdata logical: if TRUE it compares the IPD summary (marginal moments and pairwise correlations) averages over the H IPD reconstructions against the original IPD summary input values. 
#'
#' @param tabulate.similar.data if TRUE and also checkdata = TRUE it returns the full tabular comparison between the reconstructed and original IPD summaries.
#' 
#' 
#' @param SI_k resampling size of stochastic integration approach. Default to 8000.
#'
#' @param NI_tol error tolerance for numerical integration. Default 1e-02, do not decrease too much. As reasonable max value use 1e-05.
#'
#' @param NI_maxEval max number of evaluations during numerical integration. Default 500 (instead 0 implies infinite number of evaluations).
#'
#' @param input.sn.corr solution of 'correlation.matrix' into standard normal space (the Gaussian copula parameter, see Details). Default is NULL and solution is found internally. If matrix solution is instead given, it overrides internal computations and it is directly used to generate artificial data. This can be useful as a post hoc data generation tuning procedure. See Details.
#'
#'@param cp.finetune logical. If NORTA method is used and x.mode = TRUE, it iteratively fine-tunes Kruskal analytic solution (corrtype = rank.corr) of copula parameter, until the correlation bias of the generated artificial data is reduced. It can also be used along with argument 'assume.all.smooth = TRUE' (see below). Default FALSE.
#'
#'@param rescale.smoothed.binary if Kruskal analytic conversion was used and x.mode = T, it rescales smoothed binary variables into integer format (typically needed). Default FALSE.
#'
#'@param assume.all.smooth logical. If NORTA method is used, it pretends an input Pearson correlation matrix is already a valid Kruskal solution, which falsely assumes all variables are continuous, when some are actually discrete. This is biased but it can yield quick (fine-tunable -- see 'cp.finetune'). Default FALSE.
#'  
#' @return An object of class 'similar.data'.
#'
#' @details `DataRebuild()` is based on a Gaussian Copula inversion technique also known as NORmal To Anything (NORTA) transformation. Inversion occurs upon conversion of an input empirical matrix into standard normal space (copula parameter solution). If data.rearrange = "norta", conversion (optimization) expects a Pearson correlation matrix as input (corrtype = "moment.corr" is chosen automatically default). Using "norta" and "rank.corr" performs Kruskal analytic conversion (theoretically valid if all marginals are continous), whereas "normal.corr" simply returns the input matrix as it is. If optimization fails with numerical integration (default), try stochastic integration (stochastic.integration = TRUE) instead. 
#' 
#' @section Note:  this program currently assumes that previous to calculation of the input IPD summaries every IPD categorical variable with \eqn{m} levels was first converted to \eqn{m-1} dummy (binary) variables. As an alternative one can, in the future, allow for categorical marginals as well and use a Multinomial distribution modeling. This program relies on archived package 'JohnsonDistribution'.


#'
#'@seealso [Return.key.IPD.summaries()] for allowed input IPD summary format, [FitJohnsonDistribution()] from archived package JohnsonDistribution, [adaptIntegrate()] from package cubature
#' 
#' @references 
#' 
#' @examples
#' \dontrun{
#' DataRebuild( H = 100, n = 1000 )
#' }
#'
#'help("gcipdr")
#' 


DataRebuild <- function(H, n, correlation.matrix, moments, x.mode,
                        user_defined_marginals = NULL, johnson.parameters = NULL,
                        stochastic.integration = FALSE, data.rearrange = c("norta", "incomplete"),
                        corrtype = c("moment.corr", "rank.corr", "normal.corr"),
                        marg.model = c("gamma", "johnson", "user_defined"), variable.names = NULL,
                        SBjohn.correction = FALSE, compute.eec = FALSE, checkdata = FALSE,
                        tabulate.similar.data =  FALSE, SI_k = 8000, NI_tol = 1e-02, NI_maxEval = 500,
                        input.sn.corr = NULL, cp.finetune = FALSE,
                        rescale.smoothed.binary = FALSE, assume.all.smooth = FALSE)
{
    data.rearrange <- match.arg(data.rearrange)
    corrtype <- match.arg(corrtype)
    marg.model <- match.arg(marg.model)
    if (is.null(variable.names))
    {
        warning("variables are not named. Fictious names assigned")
        variable.names <- paste("V",1:dim(moments)[1], sep="") # moments rows = number of variables

    }

### get Johnson parameters from empirical moments using Hill AS 99 algorithm
    if (is.null(johnson.parameters) & marg.model == "johnson")
        johnson.parameters <- apply( cbind( moments, x.mode ), 1, function(x)
        {
            if ( !x[5] )
                FitJohnsonDistribution(x[1], x[2], x[3], x[4])
            else
                rep(NA, 6)
        }
        )

                                        
    if (is.null(correlation.matrix))
        data.rearrange <- "incomplete"  # always if corr matrix is null ALLOW FOR NULL CORR MATRIX (METHOD SET TO INCOMPLETE ALWAYS)
    rkcbool0 <- (corrtype == "moment.corr" & assume.all.smooth & is.null(input.sn.corr) & data.rearrange == "norta")
    if(rkcbool0)
        input.sn.corr <- correlation.matrix # assumed as Kruskal solution here
    Xlist <- switch(data.rearrange, # (sample) space of similar data-sets
                    incomplete= Generate.with.Incomplete.correlation(H, n,  correlation.matrix, moments,
                                                                     johnson.parameters, x.mode,
                                                                     variable.names, SBjohn.correction,
                                                                     compute.eec, corrtype, marg.model
                                                                     ),
                    norta= Generate.with.Complete.correlation(H, n, correlation.matrix, moments, user_defined_marginals,
                                                              johnson.parameters, stochastic.integration,
                                                              x.mode, variable.names, SBjohn.correction,
                                                              corrtype, marg.model, SI_k, NI_tol, NI_maxEval,
                                                              input.sn.corr, rescale.smoothed.binary
                                                              )
                    )
    rkcbool <- (cp.finetune & corrtype == "rank.corr" & any(x.mode) & rescale.smoothed.binary & data.rearrange == "norta")
    rkcbool2 <- (cp.finetune & rkcbool0)
    mex <- NA
    if (checkdata | rkcbool | rkcbool2)
        mex <- is.data.similar(Xlist$artificial.data,
                               correlation.matrix,
                               moments,
                               corrtype,
                               ifelse(rkcbool | rkcbool2, TRUE, tabulate.similar.data)
                               ) # boolean
    if ( rkcbool | rkcbool2)
    {  # new option: fine tune analytic solution (corrtype = rank.corr) of copula parameter
        ftr <- fine.tune.cp.looped(Xlist$copula.params,
                                   mex[[5]]$diff,
                                   Xlist$artificial.data,
                                   H, n,
                                   correlation.matrix,
                                   moments,
                                   johnson.parameters,
                                   x.mode, variable.names,
                                   SBjohn.correction,
                                   marg.model, corrtype,
                                   rescale.smoothed.binary
                                   )  # iterative fine tuning of copula parameter on newly generated artif data
        Xlist$artificial.data <- ftr$artificial.data  # replace with finetuned artif data
        Xlist$copula.params <- ftr$copula.params  # replace with finetuned cp
        mex <- ftr$is.similar # replace with updated data-check
    }
    out <- list(Xspace = Xlist$artificial.data,
                is.similar = mex,
                johns.params = johnson.parameters,
                copula.parameters = Xlist$copula.params
                )
    class(out) <- "similar.data"
    return(  out  )
}


### NEW: fine tuning options for analytic solution (corrtype = "rank.corr") of copula parameter

fine.tune.cp <- function(cp, diff)
{  # it assumes beta->bern conversion always underestimates original correlation
    x <- cp[lower.tri(cp)]
    except <- which(abs(cp[lower.tri(cp)]) >= 0.99)
    if (length(except) > 0)
        y <- which(abs(diff) > 0.05)[-except]
    y <- which(abs(diff) > 0.05)
    x[y] <- sign(x[y]) * ( abs(x[y]) + (1-abs(x[y]))*abs(diff[y]) + replicate(length(y),
                                                                              runif(1, 0.001, 0.005)
                                                                              )
    )  # fine tune based on observed bias + random noise
    out <- make.square.matrix(x, dim(cp)[1])
    if (!is.SPD.matrix(out)) # check spd condition
        out <- make.matrix.SPD.if.NOT.flagged( out[lower.tri(out)], dim(cp)[1] )
    return( out )
}

#  iterative fine tuning on updated artificial data generations

fine.tune.cp.looped <- function(cp, diff, ad,
                                H, n,
                                correlation.matrix, moments,
                                johnson.parameters, x.mode,
                                variable.names, SBjohn.correction,
                                marg.model, corrtype, rescale = FALSE)
{
    i <- 0
    checkad <- NULL
    while (any(abs(diff) > 0.055 ) )
    {
        cp <- fine.tune.cp(cp, diff) # fine tuned analytic solution of copula parameter
        print(paste("fine tuning iteration", i))
        updated <- Generate.with.Complete.correlation(H, n,
                                                      correlation.matrix, moments,
                                                      johnson.parameters, x.mode = x.mode,
                                                      variable.names = variable.names,
                                                      SBjohn.correction = SBjohn.correction,
                                                      marg.model = marg.model, input.sn.corr = cp,
                                                      rescale.smoothed.binary = rescale
                                                      ) # regenerate data with fine tuned copula parameter
        ad <- updated$artificial.data
        checkad <- is.data.similar(ad, correlation.matrix, moments, corrtype, tabulate = TRUE )
        diff <- checkad[[5]]$diff # updated correlation bias on new artif data. If bias is still too big, fine tunes cp again.
        i <- i+1
        if (i > 10)
                                        #   if (any(abs(cp[lower.tri(cp)]) >= 0.99) | i > 10) # tuning can enter loop where diff becomes almost constant
            break
    }
    return( list( copula.params = cp, artificial.data = ad, is.similar = checkad ) )
}
 

#### QUICK CHECK on generated artificial data
## data check: is.data.similar (Yes/No). data is similar in the sense of Definition xx of Thesis

is.data.similar <- function(Xspace, correlation.matrix, moments,
                            corrtype, tabulate = FALSE, split = FALSE)
{  # for corrfunc use
    if ( !is.null(correlation.matrix) )
    {
        detected.corrtype <- attr(correlation.matrix, "corr.type", T)
        if (detected.corrtype != corrtype)
            stop("Correlation-matrix type and argument corrtype do not match")
        ltRx <- lower.tri(correlation.matrix)
        lowtri <- do.call(
            "cbind",  # Monte Carlo raw correlations
            mclapply(Xspace, function(x)
                Return.correlation.matrix( x, corrtype)[ltRx]
                )
        )  # iteratively retrieve correlation from Xspace
        lowtri.mc.average <- apply(lowtri, 1, mean, na.rm = T) # monte carlo average
        norm0 <- lowtri.mc.average - correlation.matrix[ltRx]  # norm : mc average - true estimates
        norm0b <- sign(lowtri.mc.average) == sign(correlation.matrix[ltRx])
    }
    else
    {
        norm0 <- norm0b <- NA
        lowtri.mc.average <- ltRx <- NULL
    }
    first.moment <- do.call(
        "cbind",   # extract mc first moment
        mclapply(Xspace, function(x) apply(x, 2, mean, na.rm = T)  )
    )
    second.moment <- do.call(
        "cbind",  # extract mc second moment
        mclapply(Xspace, function(x) apply(x, 2, sd, na.rm = T)  )
    )
    third.moment <- do.call(
        "cbind",
        mclapply(Xspace, function(x) apply(x, 2, skewness, na.rm = T)  )
    )
    fourth.moment <- do.call(
        "cbind",
        mclapply(Xspace, function(x) apply(x, 2, kurtosis,na.rm = T)  )
    )

    first.moment.mc.average <- apply(first.moment, 1, mean, na.rm = T)  
    second.moment.mc.average <- apply(second.moment, 1, mean, na.rm = T)
    third.moment.mc.average <- apply(third.moment, 1, mean, na.rm = T)
    fourth.moment.mc.average <- apply(fourth.moment, 1, mean, na.rm = T)
    if ( any(is.na(first.moment.mc.average)) )  # valid for all moments ...
        warning("some variables are NaN in all MC runs: MC average is also NaN")
    norm1 <- first.moment.mc.average - moments[, 1]
    norm2 <- second.moment.mc.average - moments[, 2]
    norm3 <- third.moment.mc.average - moments[, 3]
    norm4 <- fourth.moment.mc.average - moments[, 4]
    mex0 <- any(abs(norm0) > 0.05, na.rm = T) | !all(norm0b, na.rm = T)
    mex1 <- any(abs(norm1) > 1, na.rm = T )
    mex2 <- any(abs(norm2) > 1, na.rm = T)
    mex3 <- any(abs(norm3) > 1, na.rm = T)
    mex4 <- any(abs(norm4) > 1, na.rm = T)
    ltRxdat <- data.frame(ipd = correlation.matrix[ltRx],
                          mc.mean = lowtri.mc.average,
                          diff = norm0
                          )
    rownames(ltRxdat) <- corrPairByName(rownames(correlation.matrix)) # corr pair names

    moms1 <- data.frame( ipd = moments[ ,1] , mc.mean = first.moment.mc.average , diff = norm1 )
    moms2 <- data.frame( ipd = moments[ ,2] , mc.mean = second.moment.mc.average, diff = norm2  )
    moms3 <- data.frame( ipd = moments[ ,3] , mc.mean = third.moment.mc.average , diff = norm3 )
    moms4 <- data.frame( ipd = moments[ ,4] , mc.mean = fourth.moment.mc.average, diff = norm4  )
    if ( mex1 )
    {
        print("Badly generated 1st moments:", quote = F)
        print(moms1[abs(moms1$diff) > 1, ], digits =  4, quote = F) }
    if ( mex2 )
    {
        print("Badly generated 2nd moments:", quote = F )
        print( moms2[abs(moms2$diff) > 1, ], digits = 4, quote = F)
    }
    if ( mex3 )
    {
        print("Badly generated 3rd moments:", quote = F)
        print(moms3[abs(moms3$diff) > 1, ], digits = 4, quote = F)
    }
    if ( mex4 )
    {
        print("Badly generated 4th moments:", quote = F )
        print(moms4[abs(moms4$diff) > 1, ], digits = 4, quote = F)
    }
    if ( mex0 )
    {
        print("Badly generated correlations:", quote = F )
        print(ltRxdat[abs(ltRxdat$diff) > 0.05 | !norm0b, ], digits = 4, quote = F)
    }
    if (split)
        bool <- data.frame(bool.marg1_2 = ifelse(mex1 | mex2, F, T ) ,
                           bool.marg3_4 = ifelse(mex3 | mex4, F, T ) ,
                           bool.corr = mex0 )
    else
        bool <- mex0 | mex1 | mex2 | mex3 | mex4  # does data match overall ? T | F
    if (tabulate)
        return(list( "first.moment" = moms1,
                    "second.moment" = moms2,
                    "third.moment" = moms3,
                    "fourth.moment" = moms4,
                    "lower.triangular.Rx" = ltRxdat,
                    "bool"= 1 - bool) )
    else
        return( 1 - bool) # is data similar ?
}





#' @title IPD reconstruction from IPD summaries directly computed on IPD (DEPRECATED).
#'
#' @description `Simulate.data.given.IPD()` is basically a wrapper for `DataRebuild()` but additionally assumes original IPD is available.
#'
#' @details This module is either meant for experimental purposes or when original IPD is actually available. For instance the module eases the analysis pipeline from IPD to IPD reconstruction when the user wants simply to check to utility of this package. Alternatively a network server having the IPD but not disclosure rights can think to either share the IPD summaries or a synthetic, anonymous, version of the original IPD.
#'
#'@section Note: the arguments of this function are mostly inherited by `DataRebuild()`. However fill.missing (logical, default to FALSE) is inherited by `Return.key.IPD.summaries()`. Also argument method (integer 1 to 4) sum up the available combinations from data.rearrange, corrtype, and marg.model into just four sensible options (see 'setting.comb.matrix()' in source code).
#'
#' @seealso [Return.key.IPD.summaries()] for currently allowed input IPD summary format, [DataRebuild()] 


## automatize data reconstruction from higher to lower resolution levels over one data set .. 

# similar data object, arg method = 1 to 8


Simulate.data.given.IPD <- function(data, H = NULL, method, fill.missing = F,  # TODO(me) : u may want to set fill.missing = TRUE later to ensure maximal generality
                                    SBjohn.correction = F, stochastic.integration = F,
                                    checkdata = F, compute.eec = F,
                                    tabulate.similar.data =  FALSE, print.message = TRUE,
                                    set.corr.matr2null = FALSE, SI_k = 8000, NI_tol = 1e-05,
                                    NI_maxEval = 20, input.sn.corr = NULL)
{
    method.settings.combo <-  setting.comb.matrix()  # see below ..
    data.rearrange <- method.settings.combo[1 ,method]
    corrtype <- method.settings.combo[3 ,method]
    marg.model <- method.settings.combo[2 ,method]
    md <- Return.IPD.design.matrix(data, fill.missing) # imputation option is decided here
    key.summaries <- Return.key.IPD.summaries( md, corrtype ) # also computed according to fill.missing option
    n <- key.summaries$sample.size
    moms <- key.summaries$first.four.moments
    Rx <- key.summaries$correlation.matrix
    sy <- key.summaries$model.sufficient.statistics # a list now. must be matched by variable name later
    jsp <- key.summaries$johnson.parameters
    x.mode <- key.summaries$is.binary.variable
    variable.names <- key.summaries$variable.names
    if (is.null(H))
        H <- ifelse( n < 500, 300, 100)
    if (set.corr.matr2null)
        Rx <- NULL
    data.simulation <- DataRebuild(H, n, Rx, moms,
                                   x.mode, jsp, stochastic.integration,
                                   data.rearrange, corrtype, marg.model,
                                   variable.names, SBjohn.correction,
                                   compute.eec, checkdata, tabulate.similar.data,
                                   SI_k = SI_k, NI_tol = NI_tol, NI_maxEval = NI_maxEval,
                                   input.sn.corr = input.sn.corr
                                   )
    similar.data.space <- data.simulation$Xspace
    is.data.statistically.similar <- data.simulation$is.similar  # can be NA is checkdata = FALSE
    out <- list(
        similar.data = similar.data.space,  # simulated data
        is.data.similar = is.data.statistically.similar, # is data similar to reference ?
        ipd.moments = moms,
        ipd.johnson.parameters =
            {
                if (is.null(jsp))
                    data.simulation$johns.params
                else
                    jsp
            },
        copula.parameters = data.simulation$copula.parameters,
        re.ordering.method = data.rearrange , # data simulating settings
        marginal.model.distr = marg.model ,               # data simulating settings
        corr.type = corrtype ,                # data simulating settings
        simulation.method = method, # 1 to 8 encoding of all combinations of previous three attributes..
        sample.size = n,
        is.data.imputed  = fill.missing,  # is similar data in missing data imputation mode ?
        model.sufficient.statistics = sy , # a list. slot of interest must be matched by variable names
        ipd.raw.data = data , # raw IPD-frame
        is.variable.binary = x.mode,
        data.name = attr( data, "data.name", T),  # data should possess this attribute or value will be NULL
        mc.runs = H
    )
    class(out) <- "similar.data.object"
    Print.simul.settings(method.settings.combo[, method],
                         1, H, stochastic.integration,
                         set.corr.matr2null, !print.message
                         )
    return( out )
}
                     
                    

# pick density

pick.density <- function( marg.model )
{
    density <- list(
        gamma = function(x, theta)  dgmom( x,
                                          theta[1],
                                          theta[2]
                                          ),
        johnson = function(x, theta) dJohnson(x,
                                              theta[1],
                                              theta[2],
                                              theta[3],
                                              theta[4],
                                              theta[5]
                                              )
    )
    out <- density[[marg.model]]
    return(out)
}

                                        #  graphically check convergence of marginal.model distr to empirical frequencies ...
check.convergence <- function(similar.data)
{
    if ( class(similar.data) != "similar.data.object")
        stop("function argument must of class 'similar.data.object' ")
    is.bin <- similar.data$is.variable.binary
    density.model <- similar.data$marginal.model.distr
                                        # density parameters
    if ( density.model == "gamma")
        param <- similar.data$ipd.moments[!is.bin, ] # similar data shall always have at least two columns ...
    else
        param <- t(similar.data$ipd.johnson.parameters)[!is.bin, ]
    raw.ipd <- similar.data$ipd.raw.data[ ,!is.bin]  # select non-binary variables only, may return only a vector ...
    vnames <- colnames(similar.data$ipd.raw.data)[!is.bin]
    K <- ifelse( is.null(dim(raw.ipd)[2]),
         ifelse( is.numeric(raw.ipd), 1, NA), dim(raw.ipd)[2]) # if dim null supposedly only one column present ...

                                        # roughly establish par window size
    cl <- 1
    rw <-  1     # K%/%cl + K%%cl
    i <- 0
    while( (rw*cl) < K )
    {
        cl <- cl + i
        rw <-  K%%cl # + K%/%cl
        i <- i + 1
    }
#### plot stuff
    par(mfrow=c(rw,cl))
    for(j in 1:K)
    {
                                        # recursively plot histogram (parmfrow) of non-binary variables (is.binary)
        if (K < 2)
        {
            toplot <- raw.ipd  # ...if K = 1
            dparam <- param
        }
        else
        {
            toplot <- raw.ipd[ ,j]
            dparam <- param[j, ]
        }
        density <- pick.density(density.model)
        hist( toplot, xlab = "observed IPD quantiles",
             main = paste("variable:", vnames[j]), freq = FALSE)
                                        # choose curve expression (write function switching throug densities) based on marginal.model distr arg
        curve( density(x, dparam), add = TRUE, lwd= 2, col = "red")
    }

}




##################################################################
########### LEVEL OF LOOPING : OVER SEVERAL DATA-FRAMES ###########
##################################################################

# data simulation parameters combination
setting.comb.matrix <- function()
    matrix( c(
        c( "incomplete", "gamma", "rank.corr" ),   # rank corr more natural in Rüsch meth (easier to find valid extremal correlation with binary variables ...)
        c( "incomplete", "johnson", "rank.corr" ),
        c( "norta", "gamma", "moment.corr" ),
        c( "norta", "johnson", "moment.corr" ),
        c( "incomplete", "gamma", "moment.corr" ),
        c( "incomplete", "johnson", "moment.corr" ),
        c( "norta", "gamma", "rank.corr" ),
        c( "norta", "johnson", "rank.corr" ) ),
        nrow = 3
        )   # eight combination for data simulation: create varaible 'simulation type: approach.1, approach.2, ...' DONE in bias.looped

Print.simul.settings <- function(settings, N, H,
                                 stoch, null.corr, silent = FALSE )
{
    if (!silent)
    {
        data.rearrange <- settings[1]
        corrtype <- settings[3]
        marg.model <- settings[2]
        cat("Number of data-set(s):", N,
            "\n simulated with following settings \n MC REPLICATES:",
            H, "\n CORRELATION GENERATION:",
            data.rearrange, "\n", "CORRELATION FUNCTION:",
            corrtype, "\n", "MARGINAL DISTRIBUTION:",
            marg.model, "(and bernoulli for binary variables)",
            "\n"
            )
        if (!null.corr)
        {
            if ( settings[1] == "norta")
                cat( "Stochastic integration for corr. matrix:",
                    stoch, "\n")
        }
        else
            cat("NULL corr. matrix:",
                null.corr,
                ". No  corr. structure induced: data variables are randomly ordered",
                "\n"
                )
    }
}


#' @title Reconstruction of several (unrelated) available IPDs (DEPRECATED).
#'
#' @description This function can be used as simple application of `Simulate.data.given.IPD()` on a list of IPDs, DATA_1, DATA_2, DATA_3...
#'
#' @seealso [Simulate.data.given.IPD()], [DataRebuild()]


##  loop data simulation over different data-sets, need named datalist
# (modeltype free); argument method = 1 to 8

Simulate.many.datasets <- function( datalist, H, method ,
                                   fill.missing = FALSE, # TODO(me) : later set fill TRUE
                                   SBjohn.correction = F, stochastic.integration  = F,
                                   checkdata = F, compute.eec = F,
                                   tabulate.similar.data =  FALSE,
                                   set.corr.matr2null = FALSE, SI_k = 8000,
                                   NI_tol = 1e-05, NI_maxEval = 20)
{
    K <- length( datalist )  # number of different data-sets
    if ( is.null(names(datalist)))
        datanames <- paste( "data.", 1:K, sep = "" )
    else
        datanames <- names(datalist)
    list.of.similar.data.objects <- mclapply(datalist, function(x)
        try(
            Simulate.data.given.IPD(x, H, method, fill.missing,
                                    SBjohn.correction, stochastic.integration,
                                    checkdata, compute.eec, tabulate.similar.data,
                                    print.message = FALSE, set.corr.matr2null = set.corr.matr2null,
                                    SI_k = SI_k, NI_tol = NI_tol, NI_maxEval = NI_maxEval
                                    ),
            TRUE  )  ## silent true
        )

    names(list.of.similar.data.objects) <- datanames
    Print.simul.settings(
        setting.comb.matrix()[ ,method],
        length(datalist),
        "check by data-set source",
        stochastic.integration,
        set.corr.matr2null
    )
    class(list.of.similar.data.objects) <- "list.of.similar.data.objects"
    return( list.of.similar.data.objects )

}

#

## write function to modify similar data dimensions (looped) based on predefinite features ...
# depends: Manage.datalist (initially conceived for list of IPDs) taken from datalist.all.covariates.R
# (previous to run MC_givenIPD u can select right type of data, e.g. continuous outcome for linear model)


post.trim.similar.data <- function( similar.data,
                                   action=c("only.continuous","only.binary", "mixed"),
                                   from.col=3, to.col=NULL, exact = TRUE)
{
    out <- lapply(similar.data, function(x)
        Manage.datalist(x, action, from.col, to.col, exact)
        )
    return(out)
}  


# function to check simulation approach (against given matrix above)

check.simul.approach <- function( a, b, c, ref)
    all( c(a,b,c) == ref)



#

label.simul.approach <- function(a, b, c, refmat)
{
    out <- apply(refmat, 2, function(x)
        check.simul.approach(a,b,c, x) )
    which( out )  # simulation approach code ..
    
}












 
