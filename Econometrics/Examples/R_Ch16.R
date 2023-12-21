
install.packages(c("AER","quantmod","dynlm","orcutt","nlme","stargazer"))

library(AER)
library(quantmod)
library(dynlm)
library(orcutt)
library(nlme)
library(stargazer)

# load the frozen orange juice data set
data("FrozenJuice")

# compute the price index for frozen concentrated juice
FOJCPI <- FrozenJuice[, "price"]/FrozenJuice[, "ppi"]
FOJC_pctc <- 100 * diff(log(FOJCPI))
FDD <- FrozenJuice[, "fdd"]

# convert series to xts objects
FOJCPI_xts <- as.xts(FOJCPI)
FDD_xts <- as.xts(FrozenJuice[, 3])

# Plot orange juice price index
plot(as.zoo(FOJCPI),
     col = "steelblue", 
     lwd = 2,
     xlab = "Date",
     ylab = "Price index", 
     main = "Frozen Concentrated Orange Juice")

# divide plotting area
par(mfrow = c(2, 1))

# Plot percentage changes in prices
plot(as.zoo(FOJC_pctc),
     col = "steelblue", 
     lwd = 2,
     xlab = "Date",
     ylab = "Percent",
     cex.main=0.8,
     main = "Monthly Changes in the Price of Frozen Conentrated Orange Juice")

# plot freezing degree days
plot(as.zoo(FDD),
     col = "steelblue", 
     lwd = 2,
     xlab = "Date",
     ylab = "Freezing degree days",
     cex.main=0.8,
     main = "Monthly Freezing Degree Days in Orlando, FL")


# simple regression of percentage changes on freezing degree days
orange_SR <- dynlm(FOJC_pctc ~ FDD)
coeftest(orange_SR, vcov. = vcovHAC)


# distributed lag model with 6 lags of freezing degree days
orange_DLM <- dynlm(FOJC_pctc ~ FDD + L(FDD, 1:6))
coeftest(orange_DLM, vcov. = vcovHAC)



# compute cumulative multipliers
cum_mult <-cumsum(orange_DLM$coefficients[-1])

# rename entries
names(cum_mult) <- paste(0:6, sep = "-", "period CDM")

cum_mult


# estimate cumulative dynamic multipliers using the modified regression
cum_mult_reg <-dynlm(FOJC_pctc ~ d(FDD) + d(L(FDD,1:5)) + L(FDD,6))
coef(cum_mult_reg)[-1]

# obtain coefficient summary that reports HAC standard errors
coeftest(cum_mult_reg, vcov. = vcovHAC)









# HAC

# function that computes rho tilde
acf_c <- function(x, j) {
  return(
    t(x[-c(1:j)]) %*% na.omit(Lag(x, j)) / t(x) %*% x
  )
}

# simulate time series with serially correlated errors
set.seed(1)

N <- 100

eps <- arima.sim(n = N, model = list(ma = 0.5))
X <- runif(N, 1, 10)
Y <- 0.5 * X + eps

# compute OLS residuals
res <- lm(Y ~ X)$res

# compute v
v <- (X - mean(X)) * res

# compute robust estimate of beta_1 variance
var_beta_hat <- 1/N * (1/(N-2) * sum((X - mean(X))^2 * res^2) ) / 
  (1/N * sum((X - mean(X))^2))^2

# rule of thumb truncation parameter
m <- floor(0.75 * N^(1/3))

# compute correction factor
f_hat_T <- 1 + 2 * sum(
  (m - 1:(m-1))/m * sapply(1:(m - 1), function(i) acf_c(x = v, j = i))
) 

# compute Newey-West HAC estimate of the standard error 
sqrt(var_beta_hat * f_hat_T)


# Using NeweyWest():
NW_VCOV <- NeweyWest(lm(Y ~ X), 
                     lag = m - 1, prewhite = F, 
                     adjust = T)

# compute standard error
sqrt(diag(NW_VCOV))[2]






#### DIFFERENT REPRESENTATIONS
## DL MODEL

# set seed for reproducibility
set.seed(1)

# simulate a time series with serially correlated errors
obs <- 501
eps <- arima.sim(n = obs-1 , model = list(ar = 0.5))
X <- arima.sim(n = obs, model = list(ar = 0.25))
Y <- 0.1 * X[-1] + 0.25 * X[-obs] + eps
X <- ts(X[-1])

# estimate the distributed lag model
dlm <- dynlm(Y ~ X + L(X))

# check that the residuals are serially correlated
acf(residuals(dlm))

# coefficient summary using the Newey-West SE estimates
coeftest(dlm, vcov = NeweyWest, prewhite = F, adjust = T)

## ADL MODEL

# estimate the ADL(2,1) representation of the distributed lag model
adl21_dynamic <- dynlm(Y ~ L(Y) + X + L(X, 1:2))

# plot the sample autocorrelaltions of residuals
acf(adl21_dynamic$residuals)

# compute estimated dynamic effects using coefficient restrictions
# in the ADL(2,1) representation
t <- adl21_dynamic$coefficients

c("hat_beta_1" = t[3],
  "hat_beta_2" = t[4] + t[3] * t[2])

