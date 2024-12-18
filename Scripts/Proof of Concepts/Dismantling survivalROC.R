# ========================= INVESTIGATING SURVIVALROC() ==========================
# As a poof of concept, we shall try to replicate the inner workings of the 
# survivalROC() function from Heagerty et al. (2000) in the package of the same name.
# We detected a minor error [index1 should be index2] in survivalROC(), which impacted
# the AUC-statistic materially.
# --------------------------------------------------------------------------------
# PROJECT TITLE: Default Survival Modelling
# SCRIPT AUTHOR(S): Dr Arno Botha (AB)
# VERSION: 1.0 (November-2024)
# --------------------------------------------------------------------------------
# -- Script dependencies:
#   - 0.Setup.R
#   - Comparison_tROC.R
# ================================================================================




# ----------------- 1. Preliminaries
# NOTE: Extract and prepare an example dataset, and fit a Cox regression model as baseline

# The cgd dataset from the survival package contains survival data from a clinical
# trial on patients with chronic granulomatous disease (CGD), a rare immune deficiency.
force(data(cgd,package="survival"))
data(cgd) # Load data set
# Lightly prepare data into a generic format that can span our eventual credit dataset as well
dat2 <- as.data.table(cgd)[, .(id, tstart, tstop, status, sex, age, height, weight, inherit, enum, steroids, treat)] %>% 
  rename(ID=id, Start=tstart,End=tstop,Event_Ind=status)


# --- Fit Cox Regression Model correctly, where observations are clustered around a given ID without assuming independence
coxExample2 <- coxph(Surv(Start,End,Event_Ind) ~ weight + age + enum + steroids + treat, data=dat2, id=ID)
summary(coxExample2)




# ----------------- 2. Investigate survivalROC()-function from survivalROC-package
# NOTE: Parts hereof must interactively executed after running the corresponding part 
# in the tROC()-function in script 0b(iii).FunkySurv_tROCkit.R.

# --- Assign explicitly the inner quantities found in survivalROC() function
Stime <- dat2$End; status <- dat2$Event_Ind; entry <- dat2$Start
method <- "NNE"; cut.values <- NULL; span <- 0.05; window <- "symmetric"
lambda <- NULL; predict.time <- 203; marker <- round(predict(coxExample2, type="lp"),2)
describe(marker); hist(marker, breaks="FD")


  # AB: ----------- Relevant starting point of the survivalROC()-function's logic
  times = Stime # AB: endpoints
  x <- marker
  if (is.null(entry)) 
    entry <- rep(0, length(times))
  # AB: Checking for bad records and then removing them from all quantities
  bad <- is.na(times) | is.na(status) | is.na(x) | is.na(entry)
  entry <- entry[!bad] #AB: entry times
  times <- times[!bad]
  status <- status[!bad] # event flags
  x <- x[!bad]
  # AB: posting about badness
  if (sum(bad) > 0) 
    cat(paste("\n", sum(bad), "records with missing values dropped. \n"))
  # AB: if not given, then get unique markers and use them as thresholds, and order them
  if (is.null(cut.values)) 
    cut.values <- unique(x)
  cut.values <- cut.values[order(cut.values)]
  ncuts <- length(cut.values)
  # AB: Obtain the rank order of raw end points, should they be sorted ascendantly
  ooo <- order(times)
  # AB: Now re-sort the raw end points according to these particular  rank-order indices, whereafter
  # the same is performed to the vectors of event indicators and marker-values
  times <- times[ooo]
  status <- status[ooo]
  x <- x[ooo]
  
  
  ## AB: UNNECESSARY ----------------------------------------------------------------
  # AB: Implements the classical product-limit Kaplan-Meier estimator of S(t)
  s0 <- 1
  unique.t0 <- unique(times)
  unique.t0 <- unique.t0[order(unique.t0)]
  n.times <- sum(unique.t0 <= predict.time)
  for (j in 1:n.times) {
    n <- sum(entry <= unique.t0[j] & times >= unique.t0[j])
    d <- sum((entry <= unique.t0[j]) & (times == unique.t0[j]) & 
               (status == 1))
    if (n > 0) 
      s0 <- s0 * (1 - d/n)
  }
  (s.pooled <- s0) # AB: this is the last survival probability (therefore cumulative by definition)
  ## AB: UNNECESSARY ----------------------------------------------------------------
  
  
  # AB: Pre-populating the artificial entries in the ROC-matrix
  roc.matrix <- matrix(NA, ncuts, 2)
  roc.matrix[ncuts, 1] <- 0
  roc.matrix[ncuts, 2] <- 1

  
  # AB: Largest portion of this if-then clause is dedicated to estimating S(t) using two kernels (0/1 vs exponential)
  if (method == "NNE") {
    if (is.null(lambda) & is.null(span)) {
      cat("method = NNE requires either lambda or span! \n")
      stop(0)
    }
    
    # AB: Simply obtain unique marker values and sort them similarly according to the indices of ordered unique failure times
    x.unique <- unique(x)
    x.unique <- x.unique[order(x.unique)]
    # AB: Initialise the eventual survival probability vector across unique ordered failure times
    S.t.x <- rep(0, length(x.unique))
    
    t.evaluate <- unique(times[status == 1])
    t.evaluate <- t.evaluate[order(t.evaluate)]
    t.evaluate <- t.evaluate[t.evaluate <= predict.time]
    
    # AB: UNIT TEST: 
    #all.equal(vEventTimes_Filtered, t.evaluate) #  IF TRUE, then parity achieved with our function that creates [vEventTimes_Filtered] 
    #ll.equal(vMarkers, x) #  IF TRUE, then parity achieved with our function that creates [vMarkers] 
    
    
    # AB: Starting point of the NN-estimator based method for estimating S(t) ------------------------------- 
    for (j in 1:length(x.unique)) { # AB: ----------------- Start of outer loop. 
      
      # AB: This appears to be the NN-estimator for S(t) using the 0/1-kernel, mainly in determining [wt] as the "weights"
      if (!is.null(span)) {
        if (window == "symmetric") {
          
          ddd <- (x - x.unique[j]) # AB This appears to correspond to [Diff] in our code
          n <- length(x)
          ddd <- ddd[order(ddd)] # sort again ascendantly
          # all.equal(vDiff, ddd) #  IF TRUE, then parity achieved with our function that creates [vDiff] 
          
          # AB: -- Establishing a symmetrical neighbourhood around each unique marker value
          # AB: Upper bound
          index0 <- sum(ddd < 0) + 1
          index1 <- index0 + trunc(n * span + 0.5) # Upper bound
          if (index1 > n) 
            index1 <- n
          lambda <- ddd[index1]
          wt <- as.integer(((x - x.unique[j]) <= lambda) & 
                             ((x - x.unique[j]) >= 0))
          # AB: Lower bound
          index0_b <- sum(ddd <= 0) # AB: I changed index0 to index0_b so that we can at least keep track of (and compare!) values during interactive execution within the loop
          index2 <- index0_b - trunc(n * span/2)
          if (index2 < 1) 
            index2 <- 1
          # AB: The following must surely be an error and should be index2?? Otherwise why create it? (I changed it)
          ### AB: We achieve parity in AUC between this function and our own when changing this to index2, as it should logically be anyways.
          lambda2 <- abs(ddd[index2]) # lower bound
          # AB: I also changed lambda to lambda2 here so that we can at least keep track of (and compare!) values during interactive execution within the loop
          set.index <- ((x - x.unique[j]) >= -lambda2) & 
            ((x - x.unique[j]) <= 0)
          wt[set.index] <- 1
          # AB: UNIT TEST: 
          #all.equal(weights, wt) # IF TRUE, then parity achieved with our function that creates [weights]
          
        }
        if (window == "asymmetric") {
          ddd <- (x - x.unique[j])
          n <- length(x)
          ddd <- ddd[order(ddd)]
          index0 <- sum(ddd < 0) + 1
          index <- index0 + trunc(n * span)
          if (index > n) 
            index <- n
          lambda <- ddd[index]
          wt <- as.integer(((x - x.unique[j]) <= lambda) & 
                             ((x - x.unique[j]) >= 0))
        }
      } else { # AB: This branch appears to be the exponential kernel. ignore for now 
        wt <- exp(-(x - x.unique[j])^2/lambda^2)
      }
      
      
      # AB: Calculate the necessary parts of the KM-estimator, i.e., incidence and at-risk populations
      # AB: NOTE: I have changed the iterator variable k to jj for parity with our tROC()-function
      # AB: UNIT TEST: 
      #all.equal(vStartTimes, entry) # IF TRUE, then parity achieved with our function that creates [weights]
      s0 <- 1
      for (jj in 1:length(t.evaluate)) {
        n <- sum(wt * (entry <= t.evaluate[jj]) & (times >= t.evaluate[jj]))
        d <- sum(wt * (entry <= t.evaluate[jj]) & (times == t.evaluate[jj]) * (status == 1))
        if (n > 0) 
          s0 <- s0 * (1 - d/n)
      } # AB: End of inner loop
      
      S.t.x[j] <- s0
      # AB: UNIT TEST: 
      #all.equal(S_0, s0) # IF TRUE, then parity achieved with our function that creates [weights]
      
    } # AB: ----------------- End of outer loop. 
    
    # AB: UNIT TEST: 
    #all.equal(S_t, S.t.x) # IF TRUE, then parity achieved with our function that creates [S_t]
    
    # AB: Simply allocate the estimated S(t)-values back to cases (or "markers"/rows) according to unique event times.
    # AB: Obtain an index mapping between the unique marker vector x' and the raw (non-distinct) marker vector x such that each element
    # in this mapping (corresponding to each value in x, i.e., the mapping has the same dimensions as x) denotes the index in x' that 
    # represents the unique value's position and/or rank in x'.
    vMatched2 <- match(x, x.unique)
    
    # AB: Allocate the estimated S(t)-values across subjects over prediction times t (same as unique event times), given
    # a specific marker value under which that particularly S(t) was estimated. This allocation is achieved using the 
    # previously-created map vector
    S.all.x <- S.t.x[vMatched2]
    
    
    n <- length(times)
    # AB Calculates the mean survival probability across scored cases, given an -\infty marker value
    S.marginal <- sum(S.all.x)/n
    # AB: UNIT TEST: 
    #all.equal(S.marginal, S_mean) # IF TRUE, then parity achieved with our function that creates [S_mean]
    
    # AB: Calculate TPR (true positive rate) and TNR (true negative rate) by iterating across each threshold (unique marker)
    for (c in 1:(ncuts - 1)) {
      p1 <- sum(x > cut.values[c])/n # AB: Complement distribution of 1-F_M(pc), i.e., the proportion  amongst all cases with markers greater than p_c 
      Sx <- sum(S.all.x[x > cut.values[c]])/n
      roc.matrix[c, 1] <- (p1 - Sx)/(1 - S.marginal)
      roc.matrix[c, 2] <- 1 - Sx/S.marginal
    }
    
    # AB: UNIT TEST: 
    #all.equal(roc.matrix[, 1], vTPR) # IF TRUE, then parity achieved with our function that creates [S_mean]
    #all.equal(roc.matrix[, 2], 1-vFPR) # IF TRUE, then parity achieved with our function that creates [S_mean]
  } # AB: ----------------- End of outer IF (method)
  
  
  # AB: Initialise and conversions and setting boundaries
  sensitivity = roc.matrix[, 1]
  specificity = roc.matrix[, 2]
  x <- 1 - c(0, specificity)
  y <- c(1, sensitivity)
  
  # AB: UNIT TESTS: 
  #all.equal(y, vTPR)
  #all.equal(x, vFPR)
  
  # Trapezoidal rule
  n <- length(x)
  dx <- x[-n] - x[-1]
  mid.y <- (y[-n] + y[-1])/2
  (area <- sum(dx * mid.y))
  ### AB: 57.68506% when error is fixed
  ###     56.15514% when error is left unfixed
  
  # AB: UNIT TESTS: 
  #all.equal(dx, vWidth)
  #all.equal(mid.y, vMidpoints)
  #all.equal(area, sArea)
  
  # AB: the returned output
  list(cut.values = c(-Inf, cut.values), TP = y, FP = x, predict.time = predict.time, 
       Survival = s.pooled, AUC = area)
  # AB: ----------- End of end point of the survivalROC()-function's logic
  
  
  
  # ----------------- 3. Run the function independently from above investigation for comparison
  
  # - Calculate AUC at median survival time for correctly-fitted Cox model | survivalROC
  survivalROC(Stime=dat2$End, status=dat2$Event_Ind, entry=dat2$Start, 
              method = "NNE", span=0.05, predict.time=predict.time,
              marker=round(predict(coxExample2, type="lp"),2))
  ### RESULTS: AUC: 67.06% up to time t (CD-approach)
  
  
  # --- Cleanup
  rm(Stime, status, entry, method, cut.values, span, window, lambda, predict.time, marker,
     dat2, coxExample2, x, ooo, ncuts, kmExample2, dat2, coxExample2)
  
  
  