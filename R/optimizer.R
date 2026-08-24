## optimizer.R contains optimization modules needed to perform NORTA transformation.
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




### univariate and vectorized Newton-Raphson optimization

# start is scalar
# safeguard value against f/fprime explosion, in case fprime --> zero.
# ============================================================
# MODIFIED: newtrap.one with stochastic integration awareness
# Generates z once per iteration, reuses for both fdist/fprime
# ============================================================

newtrap.one <- function(fdist, fprime, start, ..., tol=0.01,
                        maxit=50, safecheck = NULL,
                        stoch_integration = FALSE,  # ← NEW: flag for stochastic mode
                        stoch_K = NULL)             # ← NEW: K for rmvnorm
{ # safecheck shall output an attribute break.loop and modified (boolean)
    adjusted <- FALSE
    i <- 0 # to exit main loop
    j <- 0 # to exit safecheck loop
    x <- start
    # Initial evaluation (no z reuse on first call)
    out <- fdist(x, ...) # must be scalar
    while(  abs( out ) > tol  )
    {  
        # ← NEW: Generate z ONCE per iteration if stochastic
        z_current <- NULL
        if (stoch_integration)
        {
            Sigmaz <- diag(2)
            Sigmaz[1,2] <- Sigmaz[2,1] <- x  # current rz
            z_current <- rmvnorm(stoch_K, sigma=Sigmaz)
        }
        # Evaluate derivative (reuse z_current)
        den <- fprime(x, ..., precomputed_z = z_current) #  must be scalar

############# safeguarding step : try to avoid bottlenecks ####################
        if ( abs( out/den ) >= 1e5 )
        {
            x <- x - 2*out
            i <- i + 1
            if (i > maxit)
            {
                warning("Search could not escape bottleneck: procedure has failed")
                break
            }
            next
        }
#############################################

        x <- x - ( out/den )  # update x
       
############# safeguarding step : try to avoid impossible values  ####################
        if ( !is.null(safecheck) )
        {
            x <- safecheck(x, ...)  # controls that some necessay conditions for x are fulfilled in order for newtrap to proceed safely ....
            if ( attr(x, "modified", T) )
            {
                                        # has safecheck produced a change in x ? if yes proceed with further booleans ...
                j <- j + 1
                                        #  if ( j > 2 | attr(x, "break.loop", T) ){ # if u have to safecheck more than 3 times module is probably run into loop: just exit.
                if ( attr(x, "break.loop", T) )
                {
                    adjusted <- TRUE
                    break
                }
            }
        }
#############################################
        # Evaluate objective at NEW x (reuse SAME z_current — still valid for this iteration)
        out <- fdist(x, ..., precomputed_z = z_current)

        i <- i + 1
        if (i > maxit)
        {
            warning(paste("no zero found after max iteration"))
            break
        }

    }
    attributes(x)$adjusted <- adjusted  # flag if value exit on safecheck or not
    if (adjusted)
        warning("newtrap module: solution needed be adjusted to escape loop")
    return(list( value = x, iter = i ))  
}




#

newtrap.multi <- function( fdist, fprime, start, ..., tol=0.01,
                          maxit=50, safecheck = NULL)
{
                                        # f must return a vector while fprime a square matrix 

    adjusted <- FALSE
    k <- length(start)
    zero <- matrix(rep(tol, k), nrow=1)
    x <- matrix(start, nrow=1)
    i <- 0
    j <- 0
    out <- matrix( fdist(x, ...), nrow=1) # make horizontal
    den <-  fprime(x, ...)   # is square matrix (Hessian, or Jacobian)
    while( any( abs( out ) > zero )  )
    {
        den <- fprime(x, ...)
        x <- x - out%*%solve(den)   # update x

        if ( !is.null(safecheck) )
        {
            x <- safecheck(x, ...)  # controls that some necessay conditions for x are fulfilled in order for newtrap to proceed safely ....
            j <- j + 1
            if ( j > 2 | attr(x, "break.loop", T)) # if u have to safecheck 10 times module is probably run into loop: just exit.
                adjusted <- TRUE
            break
        }

        out <- matrix( fdist(x, ...), nrow = 1)
        i <- i + 1
        if(i > maxit)
        {
            warning(paste("no zero found after max iteration"))
            break
        }

    }
    if (adjusted)
        warning("newtrap module: adjusted solution to escape loop")
    return( list(value= as.numeric(x), iter= i) )
}

                                        
#
# ============================================================
# MODIFIED: newtrap wrapper passes stochastic parameters
# ============================================================
newtrap <- function(fdist, fprime, start, ..., tol=0.01, maxit=50,
                    stoch_integration = FALSE, stoch_K = NULL)  # ← NEW parameters
{
    if (length(start) > 1 )
        out <- newtrap.multi(fdist,
                             fprime,
                             start, ...,
                             tol=tol,
                             maxit=maxit
                             ) # TODO(me) : must also safeguard against singular matrixes
    else
        out <- newtrap.one(fdist, fprime, start, ..., tol=tol, maxit=maxit,
                           stoch_integration = stoch_integration,  # ← PASS THROUGH
                           stoch_K = stoch_K)                      # ← PASS THROUGH
    return(out)
}
