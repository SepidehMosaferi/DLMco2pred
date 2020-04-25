# ============================================================================================
# author: Sepideh Mosaferi
# Dynamic Linear Models
# code name: Dynamic Linear Models
# version: 0.1 April 2019
# ============================================================================================

# Required Libraries -------------------------------------------------------------------------

library("reshape2")
library("ggplot2")
library("eply")
library("fpp")
library("cowplot")
library("dplyr")
library("ggplot2")
library("dlm")
library("invgamma")
library("mcmcse")
library("coda")

setwd("/home/sepidehmosaferi/DLM")

FinalData <- read.csv("DLMdata.csv" , header = TRUE)

# Preliminary Data Analysis ------------------------------------------------------------------

par(mfrow = c(2 , 2))
ts.plot(FinalData$OBS , col = "red" , lwd = 1.5 , ylab = "CO2 Flux" , 
        ylim = c(min(FinalData$OBS) , max(FinalData$OBS)) ,
        main = "Time-Series Plot of Aggregated \n Values for CO2 Flux")
acf(FinalData$OBS , ylim = c(-1 , 1) , main = "Autocorrelation Plot")
pacf(FinalData$OBS , ylim = c(-1 , 1) , main = "Partial Autocorrelation Plot")
qqnorm(FinalData$OBS , ylim = c(-2 , 2))
qqline(FinalData$OBS , col = "blue")
shapiro.test(FinalData$OBS)

# DLM Project [Main Part Programming] --------------------------------------------------------

# DLM Analysis ------------------------------------------------------------------------------- 

buildModPoly1 <- function(v){
  dV <- exp(v[1]) 
  dW <- exp(v[2]) 
  dlmModPoly(order = 1 , dV = dV , dW = dW)}

y <- FinalData$OBS
varGuess <- var(y) 
parm <- c(log(varGuess) , log(varGuess)) 
mle <- dlmMLE(y , parm = parm , buildModPoly1)
if (mle$convergence != 0) stop(mle$message)
model <- buildModPoly1(mle$par) ; model

# MCMC1 -------------------------------------------------------------------------------------
# Explanation: 
# We need to use Filtering Backward Sampling (FFBS) Algorithm. 
# Therefore, "dlm::dlmBSample" and "dlm::dlmFilter" functions should be useful in R.

# Initial values
T <- nrow(FinalData)
sigma2 <- as.vector(model$V)
sigma2beta <- as.vector(model$W)
sh1 <- 1 + T / 2
sh2 <- 1 + T / 2
R <- 10000
FF <- as.vector(model$FF)
m0 <- as.vector(model$m0)
C0 <- as.vector(model$C0)
sigma2_save <- numeric(R)
sigma2beta_save <- numeric(R)
theta_save <- matrix(numeric(R * 107) , nrow = R , ncol = 107)

ptm <- proc.time( )
for (it in 1 : R){
  filt <- dlmFilter(y , model)
  level <- dlmBSample(filt)
  levelr <- level[-1 , ]
  rate1 <- 1 + 0.5 * sum((y - levelr)^2)
  sigma2 <- rinvgamma(1 , sh1 , rate1)
  rate2 <- 1 + 0.5 * sum(sapply(1:(nrow(level) - 1) , function(i){(level[i + 1] - level[i])^2}))
  sigma2beta <- rinvgamma(1 , sh2 , rate2)
  V(model) <- sigma2
  W(model) <- sigma2beta
  
  sigma2_save[it] <- sigma2
  sigma2beta_save[it] <- sigma2beta
  theta_save[it , ] <- levelr
}

OUTPUT <- cbind(sigma2_save , sigma2beta_save , theta_save)
proc.time( ) - ptm

library(xtable)
xtable(summary(OUTPUT[ , 1:2] %>% sqrt))
summary(OUTPUT[ , 1:2] %>% sqrt) # Important output

OUTPUT <- OUTPUT[500:10000 , ]

colnames(OUTPUT) <- c("Sigma2" , "Sigma2beta" , paste0("theta" , 1:107))
MCMC1 <- coda::mcmc(OUTPUT)
geweke.diag(MCMC1)
geweke.plot(MCMC1)

plot(MCMC1[ , 1:2] , col = "blue")
plot(MCMC1[ , c(3 , 4 , 108 , 109)] , col = "blue")
par(mfrow = c(2 , 2))
acf(OUTPUT[ , 1] , main = "Autocorrelation Plot for Sigma2")
pacf(OUTPUT[ , 1] , main = "Partial Autocorrelation Plot for Sigma2")
acf(OUTPUT[ , 2] , main = "Autocorrelation Plot for Sigma2beta")
pacf(OUTPUT[ , 2] , main = "Partial Autocorrelation Plot for Sigma2beta")

# Posterior Estimates (will be used in other parts)
MCMCMEAN <- mcmcMean(OUTPUT)
MCMCMEAN   #values with their point uncertainties
SIGMA2 <- MCMCMEAN[1 , 1]; SIGMA2sd <- MCMCMEAN[2 , 1]
SIGMA2BETA <- MCMCMEAN[1 , 2]; SIGMA2BETAsd <- MCMCMEAN[2 , 2]
THETA <- MCMCMEAN[1 , 3:109]
Uncertianty <- MCMCMEAN[2 , 3:109]

# Histograms ---------------------------------------------------------------------------------

OUTPUTa <- as.data.frame(OUTPUT)

dinvgamma1 <- function(x , a , b) dgamma(x , a , b) / x^2   #for sigma2
dsqrtinvgamma1 <- function(x , a , b) dinvgamma(x^2 , a , b) * 2 * x

dinvgamma2 <- function(x , a , b) dgamma(x , a , b) / x^2   #for sigma2beta
dsqrtinvgamma2 <- function(x , a , b) dinvgamma(x ^ 2 , a , b) * 2 * x

p1 <- ggplot(OUTPUTa , aes(Sigma2) , main = "MCMC 1") +
  geom_histogram(aes(y = ..count..) , color = "gray" , bins = 100) +
  stat_function(fun = dsqrtinvgamma1 , color = "red", geom = "area" ,
                alpha = 0.6 , args = list(a = 1 , b = 1)) +
  ggtitle("Histogram for the Posterior of Sigma2") + theme_bw( )

p2 <- ggplot(OUTPUTa , aes(Sigma2beta) , main = "MCMC 1") +
  geom_histogram(aes(y = ..density..) , color = "gray" , bins = 100) +
  stat_function(fun = dsqrtinvgamma2 , color = "red" ,
                args = list(a = 1 , b = 1)) +
  ggtitle("Histogram for the Posterior of Sigma2beta") + theme_bw( )

cowplot::plot_grid(p1 , p2 , labels = c("A" , "B") , align = "v")

# Filtering Study ----------------------------------------------------------------------------

buildModPoly1 <- function(v){
  dV <- exp(v[1]) 
  dW <- exp(v[2]) 
  dlmModPoly(order = 1 , dV = dV , dW = dW)}

y <- FinalData$OBS
varGuess <- var(y) 
parm <- c(log(varGuess) , log(varGuess)) 
mle <- dlmMLE(y , parm = parm , buildModPoly1)
if (mle$convergence != 0) stop(mle$message)
model <- buildModPoly1(mle$par) ; model

V1 <- as.vector(model$V) ; W1 <- as.vector(model$W)
V2 <- SIGMA2 ; W2 <- SIGMA2BETA

MODEL1 <- dlmModPoly(order = 1 , dV = V1 , dW = W1)
MODEL2 <- dlmModPoly(order = 1 , dV = V2 , dW = W2)
OBSFilt1 <- dlmFilter(FinalData$OBS , MODEL1)
OBSFilt2 <- dlmFilter(THETA , MODEL2)

dev.off( )
plot(FinalData$OBS , type = "o" , col = c("darkgray") , xlab = " " , ylab = "Level" ,
     main = "Plot of Filtered Values Based on MLE and Bayesian Inference")
lines(dropFirst(OBSFilt1$m) , lty = "longdash" , col = "red")  # m gives filtered values
lines(dropFirst(OBSFilt2$m) , lty = "dotdash" , col = "blue") # m gives filtered values
leg <- c("data" , "filtered values based on MLE [MSE: 0.045]" ,
         "filtered values based on Bayesian inference [MSE: 0.025]")
legend("bottomleft" , legend = leg ,
       col = c("darkgray" , "red" , "blue") ,
       lty = c("solid" , "longdash" , "dotdash") ,
       pch = c(1 , NA , NA) , bty = "n")

n <- 107 # number of obs
MSEfiltMLE <- sum((FinalData$OBS - dropFirst(OBSFilt1$m))^2) / n ; MSEfiltMLE
MSEfiltBayes <- sum((FinalData$OBS - dropFirst(OBSFilt2$m))^2) / n ; MSEfiltBayes


# Smoothing Study MLE and Bayesian -----------------------------------------------------------

V1 ; W1 # From MLE
V2 <- SIGMA2 ; W2 <- SIGMA2BETA  #From Bayesian
MODEL1 <- dlmModPoly(order = 1 , dV = V1 , dW = W1)
MODEL2 <- dlmModPoly(order = 1 , dV = V2 , dW = W2)

filtMLE1 <- dlmFilter(FinalData$OBS , MODEL1)
OBSSmoothMLE <- dlmSmooth(filtMLE1)
str(OBSSmoothMLE , 1)
attach(OBSSmoothMLE)
drop(dlmSvd2var(U.S[[n+1]] , D.S[n+1 , ]))
drop(dlmSvd2var(U.S[[n/2+1]] , D.S[n/2+1 , ]))
hwidMLE <- qnorm(0.025 , lower = FALSE) * sqrt(unlist(dlmSvd2var(U.S , D.S)))
SMOOTHMLE <- cbind(s , as.vector(s) + hwidMLE %o% c(-1 , 1)) #smoothed values

filtBayes2 <- dlmFilter(THETA , MODEL2)
OBSSmoothBayes <- dlmSmooth(filtBayes2)
str(OBSSmoothBayes , 1)
attach(OBSSmoothBayes)
drop(dlmSvd2var(U.S[[n+1]] , D.S[n+1 , ]))
drop(dlmSvd2var(U.S[[n/2+1]] , D.S[n/2+1 , ]))
hwidBayes <- qnorm(0.025 , lower = FALSE) * sqrt(unlist(dlmSvd2var(U.S , D.S)))
SMOOTHBayes <- cbind(s , as.vector(s) + hwidBayes %o% c(-1 , 1)) #smoothed values

plot(FinalData$OBS , type = "o" , col = c("darkgray") , xlab = " " , ylab = "Level" , 
     main = "Smoothed Values and CIs Based on MLE and Bayesian" , ylim = c(-0.5 , 1.5))
lines(SMOOTHMLE[ , 1] , lty = "longdash" , col = "black")
lines(SMOOTHMLE[ , 2] , lty = "dotdash" , col = "black")
lines(SMOOTHMLE[ , 3] , lty = "dotdash" , col = "black")

lines(SMOOTHBayes[ , 1] , lty = "longdash" , col = "blue")
lines(SMOOTHBayes[ , 2] , lty = "dotdash" , col = "blue")
lines(SMOOTHBayes[ , 3] , lty = "dotdash" , col = "blue")

legend("bottomleft" , col = c("darkgray" , "black" , "black" , "blue" , "blue") ,
       lty = c("solid" , "longdash" , "dotdash" , "longdash" , "dotdash") ,
       pch = c(1 , NA , NA , NA , NA) , bty = "n" , legend = c("data" , "smoothed level [MLE]" ,
                                             "95% probability limits MLE" ,
                                             "smoothed level [Bayesian]" ,
                                             "posterior 95% probability intervals"))

SmoothOBSMLE <- SMOOTHMLE[ , 1][-1]
MSE2 <- sum((FinalData$OBS - SmoothOBSMLE)^2) / n ; MSE2

# Ploting Both the Smoothed Values and Filtered Values [MLE] ---------------------------------

plot(FinalData$OBS , type = "o" , col = c("darkgray") , xlab = " " , ylab = "Level" ,
     main = "Comparison of Smoothed and Filtered Values Based on MLE")
lines(SMOOTHMLE[ , 1] , lty = "longdash" , col = "blue")
lines(dropFirst(OBSFilt1$m) , lty = "dotdash" , col = "red")

legend("bottomleft" , col = c("darkgray" , "blue" , "red") ,
       lty = c("solid" , "longdash" , "dotdash") ,
       pch = c(1 , NA , NA) , bty = "n" , legend = c("data" , "smoothed level (MSE = 0.049)" ,
                                       "filtered level (MSE = 0.006)"))

# Ploting Both the Smoothed Values and Filtered Values [Bayesian Inference] ------------------

plot(FinalData$OBS , type = "o" , col = c("darkgray") , xlab = " " , ylab = "Level" ,
     main = "Comparison of Smoothed and Filtered Values Based on Bayesian Inference")
lines(SMOOTHBayes[ , 1] , lty = "longdash" , col = "blue")
lines(dropFirst(filtBayes2$m) , lty = "dotdash" , col = "red")

legend("bottomleft" , col = c("darkgray" , "blue" , "red") ,
       lty = c("solid" , "longdash" , "dotdash") ,
       pch = c(1 , NA , NA) , bty = "n" , legend = c("data" , "smoothed level (MSE = 0.05)" ,
                                       "filtered level (MSE = 0.035)"))

SmoothOBSBayes <- SMOOTHBayes[ , 1][-1]
MSE3 <- sum((FinalData$OBS - SmoothOBSBayes)^2) / n ; round(MSE3 , 3)

# Forecasting (Ten Step Ahead) [MLE and Bayesian Approaches] ---------------------------------

SIGMA2 <- 0.07527344
SIGMA2BETA <- 0.06392886
THETA <- MCMCMEAN[1 , 3:109]

V2 <- SIGMA2 ; W2 <- SIGMA2BETA  #From Bayesian
MODEL1 <- dlmModPoly(order = 1 , dV = V1 , dW = W1)
MODEL2 <- dlmModPoly(order = 1 , dV = V2 , dW = W2)
aMLE <- window(cbind(FinalData$OBS , OBSFilt1$f) , start = 2 , end = 107)

filtBayes <- dlmFilter(THETA , MODEL2)
aBayes <- window(cbind(FinalData$OBS , filtBayes$f) , start = 2 , end = 107)
set.seed(136)
filtMLE <- dlmFilter(FinalData$OBS , MODEL1)

expdForeMLE <- dlmForecast(filtMLE , nAhead = 10 , sampleNew = 10)
expdForeBayes <- dlmForecast(filtBayes , nAhead = 10 , sampleNew = 10)  #10 step aheads Bayesian
expdForeBayes$f <- sapply(1:10 , function(i){mean(expdForeMLE$newObs[[i]])})

plot(aBayes[ , 1] , type = "o" , col = "darkgray" , ylab = "Level" , 
     xlim = c(0 , 120) , ylim = c( 0 , 1) ,
     main = "Ten Step Ahead Forcast Based on Bayesian Inference" ,
     xlab = "List of Variation of Future OBS: (0.184, 0.247, 0.311, 0.374, 0.438, 0.502, 0.565, 0.629, 0.692, 0.756)")
FUTUREBayes1 <- c(as.vector(aBayes[ , 2]))

FUTUREBayes2 <- c(as.vector(expdForeBayes$f))
VarBayes2 <- c(as.vector(expdForeBayes$Q)); unlist(VarBayes2)

lines(FUTUREBayes1 , lty = "longdash" , col = "red" , type = "o")
lines(108:117 , FUTUREBayes2 , col = "blue" , type = "o")

leg <- c("data" , "ten-step-ahead forcast")
legend("bottomleft" , legend = leg ,
       col = c("darkgrey" , "blue") ,
       lty = c("solid" , "longdash") ,
       pch = c(1 , NA) , bty = "n")

# The Innovation Process and Model Checking MLE and Bayesian--------------------------------

sum(residuals(filtBayes , sd = FALSE)^2)
shapiro.test(as.vector(residuals(filtBayes)$res))

# END ========================================================================================

