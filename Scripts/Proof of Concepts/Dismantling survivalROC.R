# ========================= INVESTIGATING SURVIVALROC() ==========================
# As a poof of concept, we shall try to replicate the inner workings of the 
# survivalROC() function from Heagerty et al. in the package of the same name.
# --------------------------------------------------------------------------------
# PROJECT TITLE: Default Survival Modelling
# SCRIPT AUTHOR(S): Dr Arno Botha (AB)

# VERSION: 1.0 (November-2024)
# ================================================================================



# -------------- 1. Preliminaries
# The cgd dataset from the survival package contains survival data from a clinical
# trial on patients with chronic granulomatous disease (CGD), a rare immune deficiency.
force(data(cgd,package="survival"))
data(cgd) # Load data set
# Lightly prepare data into a generic format that can span our eventual credit dataset as well
dat <- as.data.table(cgd)[, .(id, tstart, tstop, status, sex, age, height, weight, inherit, enum, steroids, treat)] %>% 
  rename(ID=id, Start=tstart,End=tstop,Event_Ind=status)


# --- Fit Cox Regression Model correctly, where observations are clustered around a given ID without assuming independence
coxExample <- coxph(Surv(Start,End,Event_Ind) ~ weight + age + enum + steroids + treat,
                    data=dat, id=ID)
summary(coxExample)



# --- Assign explicitly the inner quantities found in survivalROC() function
Stime <- dat$End; status <- dat$Event_Ind; entry <- dat$Start
method <- "NNE"; cut.values <- NULL; span <- 0.05; window <- "symmetric"
lambda <- NULL; predict.time <- 334; marker <- round(predict(coxExample, type="lp"),2)
describe(marker); hist(marker, breaks="FD")


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
  (s.pooled <- s0)
  
  # AB: Prepopulating the artificial entries in the ROC-matrix
  roc.matrix <- matrix(NA, ncuts, 2)
  roc.matrix[ncuts, 1] <- 0
  roc.matrix[ncuts, 2] <- 1
  
  
  if (method == "NNE") {
    if (is.null(lambda) & is.null(span)) {
      cat("method = NNE requires either lambda or span! \n")
      stop(0)
    }
    x.unique <- unique(x)
    x.unique <- x.unique[order(x.unique)]
    S.t.x <- rep(0, length(x.unique))
    t.evaluate <- unique(times[status == 1])
    t.evaluate <- t.evaluate[order(t.evaluate)]
    t.evaluate <- t.evaluate[t.evaluate <= predict.time]
    for (j in 1:length(x.unique)) {
      if (!is.null(span)) {
        if (window == "symmetric") {
          ddd <- (x - x.unique[j])
          n <- length(x)
          ddd <- ddd[order(ddd)]
          index0 <- sum(ddd < 0) + 1
          index1 <- index0 + trunc(n * span + 0.5)
          if (index1 > n) 
            index1 <- n
          lambda <- ddd[index1]
          wt <- as.integer(((x - x.unique[j]) <= lambda) & 
                             ((x - x.unique[j]) >= 0))
          index0 <- sum(ddd <= 0)
          index2 <- index0 - trunc(n * span/2)
          if (index2 < 1) 
            index2 <- 1
          lambda <- abs(ddd[index1])
          set.index <- ((x - x.unique[j]) >= -lambda) & 
            ((x - x.unique[j]) <= 0)
          wt[set.index] <- 1
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
      }
      else {
        wt <- exp(-(x - x.unique[j])^2/lambda^2)
      }
      s0 <- 1
      for (k in 1:length(t.evaluate)) {
        n <- sum(wt * (entry <= t.evaluate[k]) & (times >= 
                                                    t.evaluate[k]))
        d <- sum(wt * (entry <= t.evaluate[k]) & (times == 
                                                    t.evaluate[k]) * (status == 1))
        if (n > 0) 
          s0 <- s0 * (1 - d/n)
      }
      S.t.x[j] <- s0
    }
    S.all.x <- S.t.x[match(x, x.unique)]
    n <- length(times)
    S.marginal <- sum(S.all.x)/n
    for (c in 1:(ncuts - 1)) {
      p1 <- sum(x > cut.values[c])/n
      Sx <- sum(S.all.x[x > cut.values[c]])/n
      roc.matrix[c, 1] <- (p1 - Sx)/(1 - S.marginal)
      roc.matrix[c, 2] <- 1 - Sx/S.marginal
    }
  }
  sensitivity = roc.matrix[, 1]
  specificity = roc.matrix[, 2]
  x <- 1 - c(0, specificity)
  y <- c(1, sensitivity)
  n <- length(x)
  dx <- x[-n] - x[-1]
  mid.y <- (y[-n] + y[-1])/2
  (area <- sum(dx * mid.y))
  
  
  
  
  # --- Unit test:
  
  # - Calculate AUC at median survival time for correctly-fitted Cox model | survivalROC
  survivalROC(Stime=dat$End, status=dat$Event_Ind, entry=dat$Start, 
              method = "NNE", span=0.05, predict.time=334,
              marker=round(predict(coxExample, type="lp"),2))
  ### RESULTS: Survival estimate at median survival time = 62% .. (should be 50%)
  #            This already shows the bias of neglecting the ID-variable in clustering observations 
  #            around each relevant subject. AUC: 67.06%
  
  
  # --- Cleanup
  rm(Stime, status, entry, method, cut.values, span, window, lambda, predict.time, marker)
  
  
  