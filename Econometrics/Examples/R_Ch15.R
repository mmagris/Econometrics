library(AER)
library(dynlm)
library(forecast)
library(readxl)
library(stargazer)
library(scales)
library(quantmod)
library(urca)

rm(list=ls())

# Stuff for importing and making the data to use
USMacroSWQ <- read_xlsx("C:/Users/Martin/Desktop/us_macro_quarterly.xlsx", sheet = 1, col_types = c("text", rep("numeric", 9)))
USMacroSWQ$...1 <- as.yearqtr(USMacroSWQ$...1, format = "%Y:0%q")
colnames(USMacroSWQ) <- c("Date", "GDPC96", "JAPAN_IP", "PCECTPI", "GS10", "GS1", "TB3MS", "UNRATE", "EXUSUK", "CPIAUCSL")
GDP <- xts(USMacroSWQ$GDPC96, USMacroSWQ$Date)["1960::2013"]
GDPGrowth <- xts(400 * log(GDP/lag(GDP)))
GDPGRSub <- GDPGrowth["1962::2012"]
GDPGR_level <- as.numeric(GDPGRSub[-1])
# 3-month Treasury bills interest rate
TB3MS <- xts(USMacroSWQ$TB3MS, USMacroSWQ$Date)["1960::2012"]
# 10-year Treasury bonds interest rate
TB10YS <- xts(USMacroSWQ$GS10, USMacroSWQ$Date)["1960::2012"]
# term spread
TSpread <- TB10YS - TB3MS
# convert growth and spread series to ts objects
GDPGrowth_ts <- ts(GDPGrowth, start = c(1960, 1), end = c(2013, 4), frequency = 4)
TSpread_ts <- ts(TSpread, start = c(1960, 1), end = c(2012, 4), frequency = 4)
# join both ts objects





# BIC
# compute BIC for AR model objects of class 'dynlm'
BIC <- function(model) {
  
  ssr <- sum(model$residuals^2)
  t <- length(model$residuals)
  npar <- length(model$coef)
  
  return(
    round(c("p" = npar - 1,
            "BIC" = log(ssr/t) + npar * log(t)/t,
            "Adj.R2" = summary(model)$adj.r.squared), 4)
  )
}


# apply the BIC() to an intercept-only model of GDP growth
BIC(dynlm(ts(GDPGR_level) ~ 1))

#loop BIC over models of different orders
order <- 1:6

BICs <- sapply(order, function(x) BIC(dynlm(ts(GDPGR_level) ~ L(ts(GDPGR_level), 1:x))))
BICs
BICs[, which.min(BICs[2, ])]



# loop 'BIC()' over multiple ADL models 
order <- 1:12

BICs <- sapply(order, function(x) BIC(dynlm(GDPGrowth_ts ~ L(GDPGrowth_ts, 1:x) + L(TSpread_ts, 1:x),
                                            start = c(1962, 1), end = c(2012, 4))))
BICs
BICs[, which.min(BICs[2, ])]













## ADF
# repetitions
N <- 1000

# observations
n <- 1000

# define constant, trend and rho
drift <- 0.5
trend <- 1:n
rho <- 1

# function which simulates an AR(1) process
AR1 <- function(rho) {
  out <- numeric(n)
  for(i in 2:n) {
    out[i] <- rho * out[i-1] + rnorm(1)
  }
  return(out)
}

# simulate from DGP with constant 
RWD <- ts(replicate(n = N, drift + AR1(rho)))

# compute ADF test statistics and store them in 'ADFD'
ADFD <- numeric(N)

for(i in 1:ncol(RWD)) {
  ADFD[i] <- summary(
    dynlm(diff(RWD[, i], 1) ~ L(RWD[, i], 1)))$coef[2, 3]
}

# simulate from DGP with constant and trend
RWDT <- ts(replicate(n = N, drift + trend + AR1(rho)))

# compute ADF test statistics and store them in 'ADFDT'
ADFDT <- numeric(N)

for(i in 1:ncol(RWDT)) {
  ADFDT[i] <- summary(
    dynlm(diff(RWDT[, i], 1) ~ L(RWDT[, i], 1) + trend(RWDT[, i]))
  )$coef[2, 3]
}

# estimate quantiles for ADF regression with a drift
round(quantile(ADFD, c(0.1, 0.05, 0.01)), 2)

# estimate quantiles for ADF regression with drift and trend
round(quantile(ADFDT, c(0.1, 0.05, 0.01)), 2)

# plot standard normal density
curve(dnorm(x), 
      from = -6, to = 3, 
      ylim = c(0, 0.6), 
      lty = 2,
      ylab = "Density",
      xlab = "t-Statistic",
      main = "Distributions of ADF Test Statistics",
      col = "darkred", 
      lwd = 2)

# plot density estimates of both Dickey-Fuller distributions
lines(density(ADFD), lwd = 2, col = "darkgreen")
lines(density(ADFDT), lwd = 2, col = "blue")

# add a legend
legend("topleft", 
       c("N(0,1)", "Drift", "Drift+Trend"),
       col = c("darkred", "darkgreen", "blue"),
       lty = c(2, 1, 1),
       lwd = 2)


## DOES GDP HAVE A UNIT ROOT?
# generate log GDP series
LogGDP <- ts(log(GDP["1962::2012"]))

# estimate the model
coeftest(
  dynlm(diff(LogGDP) ~ trend(LogGDP, scale = F) + L(LogGDP) 
        + diff(L(LogGDP)) + diff(L(LogGDP), 2)))

#Compare -2.3119  against
round(quantile(ADFDT, c(0.1, 0.05, 0.01)), 2)

# test for unit root in GDP using 'ur.df()' from the package 'urca'
summary(ur.df(LogGDP, type = "trend", lags = 2, selectlags = "Fixed"))






# QLR

#set up a range of possible break dates
tau <- seq(1970, 2005, 0.25)

# initialize vector of F-statistics
Fstats <- numeric(length(tau))

# estimation loop over break dates
for(i in 1:length(tau)) {
  
  # set up dummy variable
  D <- time(GDPGrowth_ts) > tau[i]
  
  # estimate ADL(2,2) model with intercations
  test <- dynlm(GDPGrowth_ts ~ L(GDPGrowth_ts) + L(GDPGrowth_ts, 2) + 
                  D*L(TSpread_ts) + D*L(TSpread_ts, 2),
                start = c(1962, 1), 
                end = c(2012, 4))
  
  # compute and save the F-statistic
  Fstats[i] <- linearHypothesis(test, 
                                c("DTRUE=0", "DTRUE:L(TSpread_ts)", 
                                  "DTRUE:L(TSpread_ts, 2)"),
                                vcov. = sandwich)$F[2]
  
}

# identify QLR statistic
QLR <- max(Fstats)
QLR

# identify the time period where the QLR-statistic is observed
as.yearqtr(tau[which.max(Fstats)])

# series of F-statistics
Fstatsseries <- ts(Fstats, 
                   start = tau[1], 
                   end = tau[length(tau)], 
                   frequency = 4)

# plot the F-statistics 
plot(Fstatsseries, 
     xlim = c(1960, 2015),
     ylim = c(1, 7.5),
     lwd = 2,
     col = "steelblue",
     ylab = "F-Statistic",
     xlab = "Break Date",
     main = "Testing for a Break in GDP ADL(2,2) Regression at Different Dates",
     cex.main=0.8)

# dashed horizontal lines for critical values and QLR statistic
abline(h = 4.71, lty = 2)
abline(h = 6.02, lty = 2)
segments(0, QLR, 1980.75, QLR, col = "darkred")
text(2010, 6.2, "1% Critical Value",cex=0.8)
text(2010, 4.9, "5% Critical Value",cex=0.8)
text(1980.75, QLR+0.2, "QLR Statistic",cex=0.8)
