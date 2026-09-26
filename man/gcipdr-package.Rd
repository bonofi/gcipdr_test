\name{gcipdr-package}
\alias{gcipdr-package}
\alias{gcipdr}
\docType{package}
\title{
\packageTitle{gcipdr}
}
\description{
\packageDescription{gcipdr}
}

\details{
This package is a beta version and extensive editing, documenting, and
amelioration shall be expected in the future.
}

\references{
 
}

\keyword{ package }
\keyword{ copula,  }
\seealso{
 \code{\link[<pkg>:<pkg>-package]{https://github.com/bonorico/gcipdr}} 
}
\examples{

# eventually set global parallelization options with 'options(mc.cores='your numebr of cores')'

##### LEVEL 1: RAW DATA HOLDER (RDH) COMPUTES IPD SUMMARY AND DISCLOSE IT


IPD <- airquality[, -6]  # drop redundat variable Day 

# RDH gets IPD in suitable format for analysis 

md <- Return.IPD.design.matrix(IPD)  # do you want to impute NA's later ?

attributes(md)

# COMPRESSION PHASE: RDH computes IPD summary 


IPDsummary <- Return.key.IPD.summaries( md, "moment.corr" )

class(IPDsummary)

names(IPDsummary)

## summaries to be disclosed

moms <- IPDsummary$first.four.moments   

corr <- IPDsummary$correlation.matrix

corr[lower.tri(corr)]  # inspect

n <- IPDsummary$sample.size

supp <- IPDsummary$is.binary.variable

names <- IPDsummary$variable.names

##### LEVEL 2: SUMMARY DATA HOLDER (SDH) RETRIEVE IPD SUMMARY ONLY AND SIMULATES IPD FROM IT



 H <- 300  ## number of artificial IPD copies to be generated


    set.seed(8736, "L'Ecuyer")

 ## Gaussian copula based IPD simulation

 
 system.time(
     
 IPDstar <- DataRebuild( H, n, correlation.matrix = corr, moments = moms,
      x.mode = supp, data.rearrange = "norta",
              corrtype = "moment.corr", marg.model = "gamma",
     variable.names = names, checkdata = TRUE, tabulate.similar.data = TRUE  )

 )



 class(IPDstar)
 names(IPDstar)

 # diagnostic: is input IPD summary well reproduced by artificial IPD ?

 IPDstar$is.similar



 ### REPRODUCE AN ORIGINAL IPD INFERENCE FROM ARTIFICIAL IPD



   formula <- Ozone ~ Solar.R + Wind + Temp + Month

 origIPDinfer <- glm( formula, data = IPD  )

 oom <- summary(origIPDinfer)
  origbeta <- oom$coef[, 1]
 origsd <- oom$coef[, 2]
origaic <- oom$aic

#

 stat <- function( x ){

     x <- as.data.frame(x)
     
 model <- glm( formula, data = x )

     out <- summary(model)
     
  return(
      list( coef= out$coefficient[, 1], sd = out$coefficient[, 2], aic = out$aic )
      )
    }



 ### APPLY STATISTICAL PROCEDURE ON EACH IPD ARTIFICIAL COPY


 IPDinferstar <- lapply( IPDstar$Xspace, function(x) stat(x) )



### extract wanted inferences

  IPDinferstar.out <- lapply(c("coef", "sd", "aic"), function(j){
do.call("rbind",
       lapply(1:H, function(i){

           IPDinferstar[[i]][[j]]
           
       }) )
    }
       )  ; names(IPDinferstar.out) <- c("coef", "sd", "aic")


 ### SUMMARY OPERATION ON SIMULATED INFERENCES

MCmeanINFER <- lapply(IPDinferstar.out, function(x) apply(x, 2, mean, na.rm = TRUE))


### MC standard deviation of coef only 

MCsdINFER <- apply(IPDinferstar.out$coef, 2, sd, na.rm = TRUE)



 #### compare inferences !!!!

data.frame( names= c(rep(names,2), "") ,
 est= c(rep("logbeta", length(origbeta)),
rep("logsdbeta", length(origbeta)), "aic") ,
originalIPD = c( origbeta, origsd, origaic),
reconstructed = unlist(MCmeanINFER) )



}
