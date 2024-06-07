# Install and load required packages
install.packages(c("urca", "vars", "tseries", "ggplot2", "dplyr", "tidyr", "lmtest", "bptest"))
library(urca)
library(vars)
library(tseries)
library(ggplot2)
library(dplyr)
library(tidyr)
library(lmtest)
library(sandwich)

# Read and preprocess data
BTC <- read.csv("Data/BTC-USD.csv") %>% select(-Adj.Close)
ETH <- read.csv("Data/ETH-USD.csv") %>% select(-Adj.Close)
SP500 <- read.csv("Data/^SPX.csv") %>% select(-Adj.Close)
GOLD <- read.csv("Data/Gold Futures Historical Data.csv") %>%
  mutate(Close = as.numeric(gsub(",", "", Price)),
         Volume = as.numeric(gsub("K", "", gsub(",", "", Vol.))) * 1000,
         Date = format(as.Date(Date, "%m/%d/%Y"), "%Y-%m-%d")) %>%
  select(-c(Price, Vol., Change..))

# Find common dates and filter datasets
common_dates <- Reduce(intersect, lapply(list(BTC, ETH, SP500, GOLD), `[[`, "Date"))
datasets <- lapply(list(BTC, ETH, SP500, GOLD), function(df) filter(df, Date %in% common_dates))

# Extract data frames
BTC <- datasets[[1]]
ETH <- datasets[[2]]
SP500 <- datasets[[3]]
GOLD <- datasets[[4]]

# Ensure Date is in Date format
BTC$Date <- as.Date(BTC$Date)
ETH$Date <- as.Date(ETH$Date)
SP500$Date <- as.Date(SP500$Date)
GOLD$Date <- as.Date(GOLD$Date)

# Summarize datasets
lapply(list(BTC, ETH, SP500, GOLD), summary)

# Prepare data for plots
df_close <- data.frame(Date = BTC$Date, BTC = BTC$Close, ETH = ETH$Close, GOLD = GOLD$Close, SP500 = SP500$Close) %>%
  gather(key = "Asset", value = "Close", -Date)
df_volume <- data.frame(Date = BTC$Date, BTC = BTC$Volume, ETH = ETH$Volume, GOLD = GOLD$Volume, SP500 = SP500$Volume) %>%
  gather(key = "Asset", value = "Volume", -Date)
df_volatility <- data.frame(Date = BTC$Date[-1],
                            BTC = diff(log(BTC$Close)),
                            ETH = diff(log(ETH$Close)),
                            GOLD = diff(log(GOLD$Close)),
                            SP500 = diff(log(SP500$Close))) %>%
  gather(key = "Asset", value = "Volatility", -Date)

# Plot data
ggplot(df_close, aes(x = Date, y = Close, color = Asset, group = Asset)) +
  geom_line() +
  labs(title = "Prices", x = "Date", y = "Price")

ggplot(df_volume, aes(x = Date, y = Volume, color = Asset, group = Asset)) +
  geom_line() +
  labs(title = "Volumes", x = "Date", y = "Volume")

ggplot(df_volatility, aes(x = Date, y = Volatility, color = Asset, group = Asset)) +
  geom_line() +
  labs(title = "Volatility (Daily Returns)", x = "Date", y = "Volatility")

# Remove rows with NA values for ADF test
df_close <- df_close %>% drop_na()

# Stationarity test
adf_results <- lapply(df_close %>% select(-Date) %>% split(., df_close$Asset), function(x) adf.test(x$Close))
print(adf_results)

# Differencing the data and retesting
differenced_data <- df_close %>% select(-Date) %>% split(., df_close$Asset) %>% lapply(function(x) diff(x$Close))
adf_diff_results <- lapply(differenced_data, adf.test)
print(adf_diff_results)

# Cointegration test
data_diff <- data.frame(BTC = differenced_data$BTC, ETH = differenced_data$ETH, GOLD = differenced_data$GOLD, SP500 = differenced_data$SP500)
johansen_test_diff <- ca.jo(data_diff, type = "trace", ecdet = "const", K = 2)
summary(johansen_test_diff)

# VECM and VAR model
vecm_initial <- vec2var(johansen_test_diff, r = 3)
var_initial <- VAR(data_diff, p = 2)

# Function to compute AIC
compute_aic <- function(var_model) AIC(var_model)

# Manual stepwise selection
manual_stepwise_selection <- function(var_model, data) {
  best_model <- var_model
  best_aic <- compute_aic(best_model)
  improved <- TRUE

  while (improved) {
    improved <- FALSE
    for (i in 1:ncol(data)) {
      reduced_data <- data[, -i]
      reduced_model <- VAR(reduced_data, p = 2)
      reduced_aic <- compute_aic(reduced_model)
      if (reduced_aic < best_aic) {
        best_model <- reduced_model
        best_aic <- reduced_aic
        improved <- TRUE
      }
    }
  }

  return(best_model)
}

best_var_model <- manual_stepwise_selection(var_initial, data_diff)
summary(best_var_model)

# Coefficients and diagnostics
coefficients <- coef(best_var_model)
print(coefficients)
serial_test <- serial.test(best_var_model, lags.pt = 10, type = "PT.asymptotic")
print(serial_test)

# Cointegration and adjustment coefficients
cointegration_matrix <- johansen_test_diff@V
adjustment_coefficients <- johansen_test_diff@W

cat("Long-term relationships (cointegration vectors):\n", cointegration_matrix)
cat("Short-term adjustment coefficients:\n", adjustment_coefficients)

# Ramsey-RESET test for model specification
residuals_data <- data.frame(residuals(best_var_model))

reset_tests <- lapply(names(residuals_data), function(var) {
  lm_model <- lm(residuals_data[[var]] ~ ., data = residuals_data)
  resettest(lm_model, power = 2:3, type = "fitted")
})

names(reset_tests) <- names(residuals_data)
cat("Ramsey RESET Tests:\n")
print(reset_tests)

# Breusch-Pagan and White tests for heteroscedasticity
bp_tests <- lapply(names(residuals_data), function(var) {
  lm_model <- lm(residuals_data[[var]] ~ ., data = residuals_data)
  list(bp_test = bptest(lm_model), white_test = bptest(lm_model, ~ fitted(lm_model) + I(fitted(lm_model)^2)))
})

names(bp_tests) <- names(residuals_data)
cat("Breusch-Pagan and White Tests:\n")
print(bp_tests)

# Breusch-Godfrey test for autocorrelation
bg_tests <- lapply(names(residuals_data), function(var) {
  lm_model <- lm(residuals_data[[var]] ~ ., data = residuals_data)
  bgtest(lm_model, order = 2)
})

names(bg_tests) <- names(residuals_data)
cat("Breusch-Godfrey Tests:\n")
print(bg_tests)
