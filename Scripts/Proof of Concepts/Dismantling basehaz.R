# ========================= INVESTIGATING BASEHAZ() ==========================
# As a poof of concept, we shall try to replicate the inner workings of the 
# basehaz() function from Therneau (2024) in the base package called 'survival'.
# --------------------------------------------------------------------------------
# PROJECT TITLE: Default Survival Modelling
# SCRIPT AUTHOR(S): Dr Arno Botha (AB)
# VERSION: 1.0 (December-2024)
# --------------------------------------------------------------------------------
# -- Script dependencies:
#   - 0.Setup.R
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


# ----------------- 2. Kaplan-Meier analysis

# --- Fit Kaplan-Meier (KM) nonparametric (and "empty-of-covariates") model
# Compute Kaplan-Meier survival estimates (product-limit) for main-event | Spell-level with right-censoring & left-truncation
# All competing events preclude the main event from happening and are therefore considered as censored
# ID is set as the spell key, with no stratification
kmExample2 <- survfit(Surv(time=Start, time2=End, event=Event_Ind==1,type="counting") ~ 1, 
                     id=ID, data=dat2)
summary(kmExample2)$table # overall summary statistics
### RESULTS: 76 events, with median survival probability at time 334 \in [280, 373] as a 95% Confidence Interval
(kmExample_survFitSummary <- surv_summary(kmExample2)) # Survival table
### RESULTS: Median survival time of 334 has standard error of 5.8%, which is relatively large



# --- Calculate various summary statistics
# - Chosen percentile of survival times, used later as a prediction time, e.g., 75%
# NOTE: Median would also work, though sample sizes typically dwindle at that point, which can complicate analyses
survPercentile <- 0.75
cat(paste0("Kaplan-Meier: Survival time at the xth percentile of the distribution of unique survival times: \n\t",
           percTimepoint <- min(kmExample2$time[kmExample2$surv <= survPercentile],na.rm=T)))
# - Calculate the Truncated Mean (or "Restricted Mean Survival Time" [RMST]) of survival times
# NOTE: The classical mean survival time is only defined when S(t) can reach 0
# NOTE2: Calculating RMST involves integrating the KM-curve up to the last (hence truncated) observed unique event time
#         I.e., sum the survival probabilities, as weighted by the time difference between two successive survival times
cat(paste0("Kaplan-Meier: Truncated/Restricted Mean Survival Time: \n\t", 
           sum(kmExample2$surv*diff(c(0,kmExample2$time)), na.rm=T)))
# - Retrieve the mean survival probability at the previously found percentile-basd prediction point
cat(paste0("Kaplan-Meier: Retrieved mean survival probability at a given time point (t=", percTimepoint, "): \n\t",
           percent(kmExample2$surv[max(which(kmExample2$time <= percTimepoint))], accuracy=0.01)))
# - Retrieve the mean survival probability at the last unique event time as prediction time, i.e., 
# the "pooled mean survival probability"
cat(paste0("Kaplan-Meier: Retrieved mean survival probability at last unique event point as time point (t=", max(kmExample2$time), "): \n\t",
           percent(kmExample2$surv[max(which(kmExample2$time <= max(kmExample2$time)))], accuracy=0.01)))



# --- Graphing survival and related quantities from fitted KM-model | S(t), h(t)

# -- Discrete baseline hazard function: h(t) | Empirical estimation method
# - create plotting data object
haz_dat <- data.table(Time=kmExample2$time, AtRisk_n=kmExample2$n.risk, 
                      Event_n = kmExample2$n.event, Censored_n=kmExample2$n.censor,
                      Actual_Hazard=kmExample2$n.event/kmExample2$n.risk, 
                      CumulHazard = kmExample2$cumhaz, #Nelson-Aalen estimator
                      Group="1",Surv_KM = kmExample2$surv) %>% 
  filter(Event_n > 0 | Censored_n >0) %>%
  # Discrete-time variants
  mutate(CumulHazard_Disc = cumsum(Actual_Hazard), Surv_KM_Disc = cumprod(1-Actual_Hazard)) %>% 
  mutate(Event_KM_Disc = 1-Surv_KM_Disc) %>% as.data.table()
haz_dat[, Surv_KM_Disc_prev:= shift(Surv_KM_Disc, n=1, type="lag"), by=list(Group)]
# Distributional analysis on hazard rate
describe(haz_dat$Actual_Hazard); hist(haz_dat$Actual_Hazard, breaks="FD")

# - create alternative versions for sanity checks
haz_dat[Time==Time[1], hazard2 := 1 - Surv_KM_Disc]
haz_dat[Time>Time[1], hazard2 := 1 - Surv_KM_Disc/Surv_KM_Disc_prev]

# - conduct sanity checks
all.equal(haz_dat$hazard, haz_dat$hazard2) # Should be TRUE
all.equal(haz_dat$Surv_KM, haz_dat$Surv_KM_Disc) # Should be TRUE
all.equal(haz_dat$CumulHazard, haz_dat$CumulHazard_Disc) # Should be TRUE
plot(kmExample2$time, haz_dat$CumulHazard - haz_dat$CumulHazard_Disc, type="b")
### RESULTS: The discrepancy is very small difference due to estimator method differences


# - Fit regression spline to summarise hazards
sKnots <- 10
vSpline <- spline_estimation(haz_dat$Time, haz_dat$Actual_Hazard, sKnots, 3) # Actual hazard spline
# Build data object for graphing purposes | must be same names as main graphing set
datSpline <- data.table(Time=haz_dat$Time,Actual_Hazard = vSpline)


# - Add as seperate records to graphing dataset for graphing purposes
datGraph <- rbind(haz_dat[,list(Time, HazardValue=Actual_Hazard, Type="A_Actual")],
                 data.table(Time = unique(haz_dat$Time), HazardValue=vSpline, Type="B_Spline"))

# - Graphing parameters
vCol <- brewer.pal(10, "Paired")[c(10,9)] # for h(t)
mainEventName <- "CGD"; chosenFont <- "Cambria"
vlabel <- c("A_Actual" = "Actual hazard", "B_Spline"= paste0("Cubic regression spline (knots: ", ))

# - Graph object for shorter time, informed by previous graphs
(gsurv1c_d <- ggplot(haz_dat[Time<=300,], aes(x=Time,y=Actual_Hazard)) + theme_minimal() +
    geom_point(colour=vCol2[1]) + 
    geom_line(data=datSpline, linetype="solid", colour=vCol2[2]) + 
    #geom_smooth(aes(colour=Group, fill=Group), se=T, method="loess", span=sSpan, alpha=0.25, linetype="dotted") +
    labs(y=bquote(plain(Estimated~hazard*" function ["*.(mainEventName)*"]"*~italic(h(t))*":  Kaplan-Meier (spell-level)")), 
         x=bquote(Discrete~time~italic(t)*" (months) in spell: Multi-spell")) + 
    theme(text=element_text(family=chosenFont),legend.position="bottom") + 
    scale_colour_manual(name="", values=vCol, labels=vlabel) + 
    scale_y_continuous(breaks=breaks_pretty(), label=percent) + 
    scale_x_continuous(breaks=breaks_pretty(n=8), label=comma))





# ----------------- 3. Cox PH model
# --- Fit Cox Regression Model correctly, where observations are clustered around a given ID without assuming independence
coxExample2 <- coxph(Surv(Start,End,Event_Ind) ~ weight + age + enum + steroids + treat, data=dat2, id=ID,
                     ties="efron") # Efron-method should be default, but we specify just in case for backwards compatibility
summary(coxExample2)


# --- Compare cumulative hazard functions: KM-based vs Cox
# Extract baseline hazard from Cox model, having assumed covariates are null (and therefore exp(0)=1)
vBaselineHaz1 <- basehaz(coxExample2) # calculates cumulative hazard
survfit(coxExample2)$cumhaz

# Plot Kaplan-Meier Hazard vs Baseline Hazard
plot(vBaselineHaz1$time, vBaselineHaz1$hazard, type = "l", col = "blue", lwd = 2,
     xlab = "Time", ylab = "Hazard")
title(main = "Cumulative Hazard comparison", sub= "Cox hazard vs Kaplan-Meier")
points(haz_dat$Time, haz_dat$CumulHazard_Disc, col = "red", pch = 20)  # Approximate cumulative hazard
legend("topright", legend = c("Baseline Cumulative Hazard (Cox)", "KM Hazard"), 
       col = c("blue", "red"), lwd = 2, pch = c(NA, 20))
### RESULTS: There is a clear difference between these two, as literature would suggest. h vs h_0, whereas
# Cox PH models the h, and would need an estimate of h_0. 


# --- Dismantle basehaz()

# - Obtain unique event times and predictions from a given cox model
vEventTimes <- sort(unique(dat2$End))
vPredictions <- predict(coxExample2, type="lp") # linear predictor type \lambda such that h(t, x) = h_0(t).exp(\lambda)
# Distributional analysis on predictions
describe(vPredictions); hist(vPredictions, breaks="FD")

# - Compute baseline hazard h_0(t) using Efron's method
vHazards <- sapply(vEventTimes, function(t, times, events, linPred) {
  # - Testing conditions
  # t = vEventTimes[18]; times=dat2$End; events=dat2$Event_Ind; linPred=vPredictions
  vAtRisk <- which(times > t) # Indices of at-risk subjects at t (Risk set)
  vEvents <- which(times == t & events==1) # indices of those who experienced the event together at t, i.e., tied events
  d_t <- length(vEvents) # number of tied events
  vLP_risk <- exp(linPred[vAtRisk]) # first term in denominator | risk set
  vLP_events <- exp(linPred[vEvents]) #second term in denominator | tied events set
  
  # Conditional execution based on whether ties were detected or not
  if (d_t == 1) {  # Handle no ties 
    hazard <- (d_t/ sum(vLP_risk))
    
  } else if (d_t == 0) { # Handle no events
    hazard <- 0
    
  } else { # handle ties
    # Efron-adjustment
    sum_Adjusted <- sum(vLP_risk)
    hazard <- sum(sapply(1:d_t, function(k) {
      # Adjust denominator progressively for each tied event
      currentAdj <- sum_Adjusted - ((k-1)/d_t) * sum(vLP_events)
      return( 1/ currentAdj)
    }))
    
    # Calculate and return incremental hazard at t using effron
    # d_t / sum(exp(linPred[at_risk])) # similar to Breslow
  }
  
  return(hazard)
  
}, times=dat2$End, events=dat2$Event_Ind, linPred=vPredictions)

# - Compute cumulative hazard from the incremental hazards produced by Efron's method earlier
vCumulHazard <- cumsum(vHazards)


# - Plot basehaz-results vs Bcustom implementation of Efron's method
plot(vBaselineHaz1$time, vBaselineHaz1$hazard, type = "l", col = "blue", lwd = 2,
     xlab = "Time", ylab = "Hazard")
title(main = "Cumulative Hazard comparison", sub= "Basehaz vs Own")
points(vEventTimes, vCumulHazard, col = "green", pch = 20)  # Approximate cumulative hazard
legend("topright", legend = c("Basehaz-based H_0 [Efron]", "Own Efron-estimated H_0"), 
       col = c("blue", "green"), lwd = 2, pch = c(NA, 20))
### RESULTS: The lines are different.

