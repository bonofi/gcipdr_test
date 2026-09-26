## This file contains miscellaneous functions and info on dependencies.
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



### libraries to be loaded

## url <- "https://cran.r-project.org/src/contrib/Archive/JohnsonDistribution/JohnsonDistribution_0.24.tar.gz"
## pkgFile <- "JohnsonDistribution_0.24.tar.gz"
## download.file(url = url, destfile = pkgFile)

## install.packages(pkgs=pkgFile, type="source", repos=NULL)

## unlink(pkgFile)


 warning("This package depends on CRAN archived package 'JohnsonDistribution'. Please manually download this dependency from: https://cran.r-project.org/src/contrib/Archive/JohnsonDistribution/JohnsonDistribution_0.24.tar.gz")


# warning("Other dependecies: moments, parallel, cubature, mvtnorm, ggplot2")



## ..... plus some miscellaneous functions


is.binary <- function(x, tol = .Machine$double.eps^0.5)
{
                                        # check if var is binary
    x <- na.omit(x)
    #all( abs(x - round(x)) < tol & x >= 0 & x <= 1  )
    all(x %in% c(0,1))
}  

# check if both variables in a pair are both continuous
# result: correlation matrix where all-continuous pairs have
# Kruskal analytic solution 
kruskal_init_matrix <- function(dat){
  
  p <- dim(dat)[2]
  combos <- combn(p, 2)
  is_pair_cont <- apply(
    combos, 2,
    \(x) !any(
      c(
        is.binary(dat[, x[1]]),
        is.binary(dat[, x[2]])
      )
      
    )
  )
  
  # indicator for binary variable presence
  I <- make.square.matrix(is_pair_cont, p)
  Rx <- momcor(dat)
  Kx <- kruskalconv(
    rankcor(dat)
  ) 
   
  # return matrix mixture of moment and Kruskal correlations
  
  out <- (1-I)*Rx + Kx*I 
  
  # flag all-continuous pairs as attribute
  attributes(out)$is_cont_flag <- I
  
  return(out)
}


reldist <- function(mean, true, range)
{
                                        # relative distance
    (mean-true)/range

}


                                        # function to check is a parameter is finite or available

is.inf <- function( vec )
{

   vec == Inf | vec == -Inf

}


is.inf.OR.na <- function( vec  )
{
    is.inf(vec) | is.na(vec)    

}


 ## is vector in continuous domain ??

is.continuous <- function(x)
{
    if ( is.numeric(x) | is.integer(x))
        !is.binary(x)
    else
    {
        if ( is.character(x) )
            FALSE
        else
        {
            if ( is.factor(x) )
                FALSE
        }
    }
}


# are vectors in continuous domain ?

are.continuous <- function( ds )
{
    if ( !is.null( dim(ds) ) )
        unlist( lapply( ds, function(x) is.continuous(x) ) )
    else
        is.continuous(ds)

}


## two simple utilities to map from a beta to a bit variable

sample.bit <- function(p) sapply(p, function(x) rbinom(1, 1, x) )


map2bit <- function(p) p >= quantile(p, 1 - mean(p)) # keeps both moments and approx correlations


### add rownames to lower triangular correlation matrix

corrPairByName <- function(stnames)
{
    cpnamesfull <- unlist(
        lapply(stnames, function(i)
            lapply(stnames, function(j)
            {
                paste( i, "-", j, sep = "" )

            }
            )
            )
    )
     
                                        
    cpnames <- cpnamesfull[as.vector(
        lower.tri(matrix(
            nrow = length(stnames),
            ncol=length(stnames)
        )
        )
    )
    ] 

    return( cpnames)
    
}




#### function to convert some data columns into multirow latex code to be used within xtable

# df = data.frame, cols = target columns, span = number of multirow span for each column ...

make.multirow <- function( df, cols,
                          span = rep(NA, length = dim(df)[2]),
                          em = rep(NA, length = dim(df)[2]),
                          rotate = rep(FALSE, length(cols)
                                       )
                          )
{
    newcols <- do.call(
        "cbind",
        lapply(1:length(cols), function(j)
        {
            label <- as.character(df[[ cols[j] ]]) # convert value into character
            if (is.na(span[j]))  # define row span
                span.row <- rle(label)$lengths
            else
                span.row <- span[j]
            if (is.na(em[j]))
                col.width <- "4em"
            else
                col.width <- em[j]

            first <- !duplicated(label)   # set duplicates to void
            label[!first] <- ""

                                        # define appearance of \multirow
            if( rotate[j] )
                label[first] <- paste0(
                    "\\parbox[t]{2mm}{\\multirow{",
                    span.row,
                    "}{*}{\\rotatebox[origin=c]{90}{",
                    label[first], "}}}"
                )
            else
                label[first] <- paste0(
                    "\\multirow{",
                    span.row,
                    "}{",
                    col.width,"}{",
                    label[first], "}"
                )

            return( label )
        }

        )
    )

    df[, cols] <- newcols
    return( df )
    
}



#####  ## scatterplot function with linear and loess | lowess regression



panel.fit <- function(x,y)
{
    b <- lm(y~x)$coef
    pred <- b[1] + b[2]*x
    points(x,y, pch = ".")
    lines(lowess(x,y), col= 'blue', lty = 4, lwd = 2)
    lines(x, pred , col='red', lwd = 2)       
}

                                        

panel.fit2 <- function(x,y)
{
    b <- lm(y~x)$coef
    pred <- b[1] + b[2]*x
    points(x,y, pch = ".")
    lines(x[order(x)],
          loess(y~x)$fitted[order(x)],
          col= 'blue', lty = 4, lwd = 2)
    lines(x, pred , col='red', lwd = 2)
       
}

#

panel.cor <- function(x,y)
{

    r <- cor(x,y)
    legend("center", paste(round(r,2)), bty = "n")
    
  }


                                        #

scatter <- function(dat2,
                    smooth = c("loess", "lowess"),
                    text = NULL)
{
    smooth <- match.arg(smooth)
    pairs(dat2,
          upper.panel = panel.cor,
          lower.panel =switch(smooth,
                              "lowess" = panel.fit,
                              "loess" = panel.fit2
                              ),
          main = text
          )
}

