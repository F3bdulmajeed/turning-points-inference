
################################################################################
#
# Purpose: Efficient statistical inference of turning points in animal movement data.
#
# Author: Abdulmajeed Alharbi
#
# Date created: January 11 2025
#
################################################################################



# The required packages
required_packages <- c("Rcpp", "Rfast", "zoo")

# Install missing packages
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

# Load the packages
library(Rcpp)
library(Rfast)
library(zoo)

# read the c++ functions
sourceCpp("The_algorithms.cpp")


# simulate an example 
set.seed(111)
h <- c( rvonmises(1000 , -pi/2, 10) , rvonmises(1000 , pi/2, 10) , rvonmises(1000 , -pi/2, 10) )


#############################################################################
# Stage 1: Calculate rho , kappa_rho and kappa_mu
#############################################################################

## Using secussive difference 
rho <- median(1-cos(diff(h,2))) / median(1-cos(diff(h,1))) -1 
kappa_rho <- 0.4592 / ( (1+rho)*(1-sqrt( median(sin(diff(h)))^2 + median(cos(diff(h)))^2 )))
kappa_mu <- (1-rho^2)*kappa_rho

# Using non-parameteric method
window <- 21 # Need to be adjusted.
med_head <- atan2(rollmedian(sin(h), window, na.pad = TRUE, align = "center"), 
                  rollmedian(cos(h), window, na.pad = TRUE, align = "center"))

# residuals of moving median filter
errors <- na.omit((h - med_head) %% (2 * pi))
sin_errors <- sin(errors)
cos_errors <- cos(errors)
n <- length(errors)

# Compute rho (autocorrelation)
rho <- sum(sin_errors[-1] * sin_errors[-n]) / 
        sqrt(sum(sin_errors[-1]^2) * sum(sin_errors[-n]^2))

# Calculate whitened errors
whitened_errors <- (errors[-1] - atan2(rho * sin_errors[-n], 1 - rho + rho * cos_errors[-n])) %% (2 * pi)

# Estimate kappa for whitened_h and errors
kappa_rho <- vm.mle(whitened_errors)$param[2]
kappa_mu  <- vm.mle(errors)$param[2]



#############################################################################
# Stage 2: Estimate the turning-points
#############################################################################

# fit the IID model 
penalty <-  -log(length(h))/kappa_mu   * (1+rho)/(1-rho)
tps <-  AFPOP(sin(h), cos(h), penalty) 
sort(tps)


# fit the AR1 model 
penalty <-  -log(length(h))/kappa_rho
mu <- seq(-pi,pi,  length = 720 )
cos_h <- cos(h) ; sin_h <- sin(h) ; cos_mu<- cos(mu) ; sin_mu<- sin(mu)
Msine <- Outer(sin_h , sin_mu , "*")
Mcosine <- Outer(cos_h , cos_mu , "*")
D <- cos(diff(h))
M <-  Msine + Mcosine
tps <- DFPOP(M,D,rho,penalty)
sort(tps)

