month <- 12
lambda <- dat[,.N]^(-1/3)
# - Confirm prepared dataset is loaded into memory
if (!exists('datCredit_train_TFD')) unpack.ffdf(paste0(genPath,"creditdata_train_TFD"), tempPath);gc()

# Initialize cox model
cox <- coxph(Surv(Start, End, Default_Ind) ~ BalanceToTerm, data = datCredit_train_TFD); gc() # Build a cox model based on a single variable

# Initialize dataset
dat <- datCredit_train_TFD %>% subset(select=c(Date, LoanID, Start, End, Default_Ind))
dat[, Marker := round(predict(cox,type="risk"),1)]
thresholds <- dat$Marker %>% unique() %>% sort()
nthresh <- length(thresholds)

# House keeping
rm(cox, datCredit_train_TFD); gc()

# Sort dataset according to start time
setorder(dat,End)

# Create unique endpoints before given time
uEnd <- dat$End %>% unique() %>% sort()
nTimes <- sum(uEnd <= month) # Number of months before and including the final month

# Initialize Nearest Neighbor Estimation variables
uMarker <- dat$Marker %>%  unique() %>% sort() # Obtain unique markers
uDTimes <- unique(dat$End[dat$Default_Ind == 1]) %>% sort() # Unique defaulting times
DTimes <- uDTimes[uDTimes <= month] # Unique defaults before given time
S_t <- numeric(length(uMarker)) # Vector to contain survival estimates

# Calculate survival probability for each unique marker value using weights
gc();expWeights <- lapply(uMarker,function(um) exp(-(dat$Marker - um)^2 / lambda^2))# Compute exponential weights for all uMarker values

# Calculate pooled survival probability at times using Kaplan-Meier estimator
start_before_time <- outer(dat$Start, DTimes, "<=")# Ensure observation existed before time t
end_after_time <- outer(dat$Start, DTimes, ">=")# Ensure observation exists after time t
end_at_time <- outer(dat$End, DTimes, "==")# Ensure observation exists at time t

# Calculate boolean matrices
at_risk <- (start_before_time & end_after_time) # Populations at risk at each time point
events <- (start_before_time & end_at_time & (dat$Default_Ind == 1)) # Number of defaults at each time point

for (i in 1:length(expWeights)){
  # Calculate weighted sums of n and d across uTimes
  n_values <- colSums(expWeights[[i]] * at_risk) # Weighed populations at risk
  d_values <- colSums(expWeights[[i]] * events) # Weighed loans leaving the populations
  
  # Calculate survival factors and handle division by zero cases
  survival_factors <- 1 - (d_values / n_values)
  survival_factors[is.na(survival_factors)] <- 1  # Set NaN cases to 1
  S_t[i] <- prod(survival_factors)
}

S_all <- S_t[match(dat$Marker, uMarker)]
n <- NROW(dat)
S_Marg <- sum(S_all)/n

# Preallocate ROC matrix
roc.matrix <- matrix(NA, nthresh, 2)
roc.matrix[nthresh, ] <- c(0, 1)

for (c in 1:(nthresh - 1)) {
  p1 <- sum(dat$Marker > thresholds[c])/n
  Sx <- sum(S_all[dat$Marker > thresholds[c]])/n
  roc.matrix[c, 1] <- (p1 - Sx)/(1 - S_Marg)
  roc.matrix[c, 2] <- 1 - Sx/S_Marg
}
sensitivity = roc.matrix[, 1]
specificity = roc.matrix[, 2]
x <- 1 - c(0, specificity)
y <- c(1, sensitivity)
n <- length(x)
dx <- x[-n] - x[-1]
mid.y <- (y[-n] + y[-1])/2
area <- sum(dx * mid.y)

# Create a data frame for plotting
datGraph <- data.frame(x = x, y=y)
vCol <- brewer.pal(8,"Set1")[c(2)]

# Plot ROC curve
(gg <- ggplot(datGraph,aes(x=x,y=y)) + theme_minimal() + 
    theme(text = element_text(family="Cambria")) + labs(x = "FP", y = "TP") +
    geom_line(color=vCol) + geom_abline(linetype="dashed",color="grey"))
 