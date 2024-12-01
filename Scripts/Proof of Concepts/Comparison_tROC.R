# ========================= INVESTIGATING SURVIVALROC() ==========================
# As a poof of concept, we shall compare various functions from various packages
# in conducting time-dependent ROC-analyses on the same fitted Cox regression model
# --------------------------------------------------------------------------------
# PROJECT TITLE: Default Survival Modelling
# SCRIPT AUTHOR(S): Dr Arno Botha (AB)
# VERSION: 1.0 (November-2024)
# --------------------------------------------------------------------------------
# -- Script dependencies:
#   - 0.Setup.R
# ================================================================================




# ----------------- 1. Extract and prepare an example dataset and models for unit tests
# The cgd dataset from the survival package contains survival data from a clinical
# trial on patients with chronic granulomatous disease (CGD), a rare immune deficiency.
force(data(cgd,package="survival"))
data(cgd) # Load data set
# Lightly prepare data into a generic format that can span our eventual credit dataset as well
dat <- as.data.table(cgd)[, .(id, tstart, tstop, status, sex, age, height, weight, inherit, enum, steroids, treat)] %>% 
  rename(ID=id, Start=tstart,End=tstop,Event_Ind=status)



# ----------------- 2. Conduct preliminary Kaplan-Meier analyses

# --- Fit Kaplan-Meier (KM) nonparametric (and "empty-of-covariates") model
# Compute Kaplan-Meier survival estimates (product-limit) for main-event | Spell-level with right-censoring & left-truncation
# All competing events preclude the main event from happening and are therefore considered as censored
# ID is set as the spell key, with no stratification
kmExample <- survfit(Surv(time=Start, time2=End, event=Event_Ind==1,type="counting") ~ 1, 
                     id=ID, data=dat)
summary(kmExample)$table # overall summary statistics
### RESULTS: 76 events, with median survival probability at time 334 \in [280, 373] as a 95% Confidence Interval
(kmExample_survFitSummary <- surv_summary(kmExample)) # Survival table
### RESULTS: Median survival time of 334 has standard error of 5.8%, which is relatively large

# -- Calculate various summary statistics
# - Chosen percentile of survival times, used later as a prediction time, e.g., 75%
# NOTE: Median would also work, though sample sizes typically dwindle at that point, which can complicate analyses
cat(paste0("Kaplan-Meier: Survival time at the xth percentile of the distribution of unique survival times: \n\t",
           percTimepoint <- min(kmExample$time[kmExample$surv <= 0.75],na.rm=T)))
# - Calculate the Truncated Mean (or "Restricted Mean Survival Time" [RMST]) of survival times
# NOTE: The classical mean survival time is only defined when S(t) can reach 0
# NOTE2: Calculating RMST involves integrating the KM-curve up to the last (hence truncated) observed unique event time
#         I.e., sum the survival probabilities, as weighted by the time difference between two successive survival times
cat(paste0("Kaplan-Meier: Truncated/Restricted Mean Survival Time: \n\t", 
           sum(kmExample$surv*diff(c(0,kmExample$time)), na.rm=T)))
# - Retrieve the mean survival probability at the previously found percentile-basd prediction point
cat(paste0("Kaplan-Meier: Retrieved mean survival probability at a given time point (t=", percTimepoint, "): \n\t",
           percent(kmExample$surv[max(which(kmExample$time <= percTimepoint))], accuracy=0.01)))
# - Retrieve the mean survival probability at the last unique event time as prediction time, i.e., 
# the "pooled mean survival probability"
cat(paste0("Kaplan-Meier: Retrieved mean survival probability at last unique event point as time point (t=", max(kmExample$time), "): \n\t",
           percent(kmExample$surv[max(which(kmExample$time <= max(kmExample$time)))], accuracy=0.01)))


# --- Graphing survival and related quantities from fitted KM-model | S(t), h(t)

# -- Graphing parameters
vCol <- brewer.pal(10, "Paired")[c(10)] # for S(t)
vCol2 <- brewer.pal(10, "Paired")[c(10,9)] # for h(t)
sSpan <- 0.1; # span for LOESS-smoother in h(t)
vlabel <- paste0("Loess-smoothed hazard [span: ", sSpan, "]") # for h(t)
mainEventName <- "CGD"
chosenFont <- "Cambria"

# -- Survival probability, S(t)=y
(gsurv1c_a <- ggsurvplot(kmExample, fun="pct", conf.int=T, legend="none", 
                         break.time.by=round(max(kmExample$time)/8), palette=vCol,
                         xlab = bquote(Discrete~time~italic(t)*" (months) in spell: Multi-spell"),
                         ylab = bquote(Survival~probability~"["*.(mainEventName)*"]"*~italic(S(t))*": spell-level (Kaplan-Meier)"), 
                         xlim=c(0, max(kmExample$time)+1), surv.median.line = "hv", censor=F, 
                         ggtheme = theme_bw(base_family=chosenFont), tables.theme = theme_cleantable(),
                         tables.height=0.10, tables.y.text=F, tables.y.text.col=T, risk.table = "abs_pct", risk.table.pos = "out",
                         cumevents=T, cumevents.title="Cumulative number of events", 
                         cumcensor=T, cumcensor.title="Cumulative number of censored observations (incl. competing risks)",
                         risk.table.title = "Number in (% of) sample at risk of main event", font.family=chosenFont, fontsize=2.5))


# -- Discrete baseline hazard function: h(t) | Empirical estimation method
# - create plotting data object
haz_dat <- data.table(Time=kmExample$time, AtRisk_n=kmExample$n.risk, 
                      Event_n = kmExample$n.event, Censored_n=kmExample$n.censor,
                      hazard=kmExample$n.event/kmExample$n.risk, 
                      CumulHazard = kmExample$cumhaz, #Nelson-Aalen estimator
                      Group="1",Surv_KM = kmExample$surv) %>% 
  filter(Event_n > 0 | Censored_n >0) %>%
  # Discrete-time variants
  mutate(CumulHazard_Disc = -cumsum(log(1-hazard)), Surv_KM_Disc = cumprod(1-hazard)) %>% 
  mutate(Event_KM_Disc = 1-Surv_KM_Disc) %>% as.data.table()
haz_dat[, Surv_KM_Disc_prev:= shift(Surv_KM_Disc, n=1, type="lag"), by=list(Group)]
# - create alternative versions for sanity checks
haz_dat[Time==Time[1], hazard2 := 1 - Surv_KM_Disc]
haz_dat[Time>Time[1], hazard2 := 1 - Surv_KM_Disc/Surv_KM_Disc_prev]
# - conduct sanity checks
all.equal(haz_dat$hazard, haz_dat$hazard2) # Should be TRUE
all.equal(haz_dat$Surv_KM, haz_dat$Surv_KM_Disc) # Should be TRUE
all.equal(haz_dat$CumulHazard, haz_dat$CumulHazard_Disc) # usually FALSE
plot(kmExample$time, haz_dat$CumulHazard - haz_dat$CumulHazard_Disc, type="b")
### RESULTS: The discrepancy is very small difference due to estimator method differences

# - Graph object for shorter time, informed by previous graphs
(gsurv1c_d <- ggplot(haz_dat[Time<=300,], aes(x=Time,y=hazard)) + theme_minimal() +
    geom_line(linetype="solid", colour=vCol2[1]) + geom_point(colour=vCol2[1]) + 
    geom_smooth(aes(colour=Group, fill=Group), se=T, method="loess", span=sSpan, alpha=0.25, linetype="dotted") +
    labs(y=bquote(plain(Estimated~hazard*" function ["*.(mainEventName)*"]"*~italic(h(t))*": spell-level (Kaplan-Meier)")), 
         x=bquote(Discrete~time~italic(t)*" (months) in spell: Multi-spell")) + 
    theme(text=element_text(family=chosenFont),legend.position="bottom") + 
    scale_colour_manual(name="", values=vCol2[2], labels=vlabel) + 
    scale_fill_manual(name="", values=vCol2[2], labels=vlabel) + 
    scale_y_continuous(breaks=breaks_pretty(), label=percent) + 
    scale_x_continuous(breaks=breaks_pretty(n=8), label=comma))
### RESULTS: The hazard appears to be near-constant over time, with some notable oscillation over some prediction periods.
# However, when viewed in tandem with S(t), itself almost a straight downward-sloping line, it makes sense for hazard
# to be near-flat. The oscillation also seems more pronounced towards later prediction periods than earlier ones.

# -- Save plots
dpi <- 150 # need to decrease size for risk tables' text
ggsave(print(gsurv1c_a,newpage=F), file=paste0(genFigPath,"Proof of Concepts/CoxExample1_SurvFig1c_a-", mainEventName,"_Surv-KaplanMeier-SpellLevel-MultiSpell-LatentComp-InclLeftTrunc_Correct.png"),
       width=1200/dpi, height=1000/dpi,dpi=dpi, bg="white")
dpi <- 180 # reset
ggsave(gsurv1c_d, file=paste0(genFigPath,"Proof of Concepts/CoxExample1_SurvFig1c_d-", mainEventName,"_Hazard-KaplanMeier-SpellLevel-MultiSpell-LatentComp-InclLeftTrunc_Correct.png"),
       width=1200/dpi, height=1000/dpi,dpi=dpi, bg="white")




# ----------------- 3. Fit Cox regression models
# NOTE: 'Best' model was fit interactively by changing input space until statistical significance is obtained.
# The goal is not to get the best model, but just a baseline upon which ROC-analyses

# --- Fit Cox Regression Model correctly, where observations are clustered around a given ID without assuming independence
coxExample <- coxph(Surv(Start,End,Event_Ind) ~ weight + age + enum + steroids + treat,
                    data=dat, id=ID)
summary(coxExample)

# --- Fit a Cox regression model incorrectly by ignoring the clustering and heterogeneous variance, assuming that
# all observations are subject-level (one row per subject) and therefore independent from one another (wrong)
coxExample_wrong <- coxph(Surv(Start,End,Event_Ind) ~ weight + age + enum + steroids + treat,
                          data=dat)
summary(coxExample_wrong)
### RESULTS: standard errors (as measured using the robust-variant of the log-rank test) are inflated
#   in this wrongly-fit model, as a result of assuming independence amongst observations within a cluster/spell per subject

  


# ----------------- 4. Compare ROC-analyses using the same fitted Cox model

# - Decide prediction point t at which we would like to conduct a time-dependent ROC-analysis
# assuming the cumulative/dynamic (CD) context, i.e., that TPR is calculated given that 
# raw event times T_{ijt} are less than or equal to the given prediction time t
predictTime <- 203 # 75% percentile survival time (203); median: 338


# ----- Package: survivalROC()
# Using survivalROC with NNE-estimator of S(t) under the CD-approach | NNE-kernel for S(t)

# - Calculate AUC up to given prediction time for correctly-fitted Cox model
# NOTE: Uses the superior Nearest Neighbour Estimator (NNE) method for S(t)
survivalROC::survivalROC(Stime=dat$End, status=dat$Event_Ind, entry=dat$Start, 
                         method = "NNE",  span=0.05, predict.time=predictTime,
                         marker=round(predict(coxExample, type="lp"),2))
### RESULTS: AUC: 56.16% up to t

# - Calculate AUC up to given prediction time for wrongly-fitted Cox model  | NNE-kernel for S(t)
# NOTE: Uses the superior Nearest Neighbour Estimator (NNE) method for S(t)
survivalROC::survivalROC(Stime=dat$End, status=dat$Event_Ind, entry=dat$Start, 
                         method = "NNE", span=0.05, predict.time=predictTime,
                         marker=round(predict(coxExample_wrong, type="lp"),2))
### RESULTS: AUC: 56.16% up to t. Therefore, prediction accuracy seems unaffected by model misfitting

# - Calculate AUC up to given prediction time for correctly-fitted Cox model | KM-kernel for S(t)
# NOTE: Uses the classical but flawed Kaplan-Meier (KM) estimator for S(t), purely for 
# parity with other packages
survivalROC::survivalROC(Stime=dat$End, status=dat$Event_Ind, entry=dat$Start, 
                         method = "KM", predict.time=predictTime,
                         marker=round(predict(coxExample, type="lp"),2))
### RESULTS: AUC: 59.53% up to t



# ----- Package: timeROC()
# Using timeROC with IPCW-estimator under the CD-approach

# - Calculate AUC up to given prediction time for correctly-fitted Cox model | KM-kernel of S(t)
# NOTE: This function calculates the "Inverse Probability of Censoring Weighting [IPCW]" version
# of the Cumulative/Dynamic (CD) time-dependence ROc-graph. Results may therefore differ inherently
# NOTE2: Only implements the the KM-estimator for S(t) instead of the superior
#   NNE-method from the survivalROC-package
# NOTE3: Since timeROC cannot accept an entry-argument (start time), the data has to be filtered
# such that the filtration contains only those at-risk subjects given a prediction time t, i.e., 
# the start times should at least be earlier than the t, while the end point still being beyond t.
# Doing so still enforces the CD-approach to constructing an ROC-curve and preserves parity (somewhat)
# However, doing so assumes start times of 0, which is still false. This implies a flawed approach.
datFiltered <- subset(dat, Start <= predictTime & End > predictTime)
timeROC::timeROC(T=datFiltered$End, delta=datFiltered$Event_Ind, cause=1, times=predictTime,
                 marker=round(predict(coxExample, newdata=datFiltered, type="risk"),2))
### RESULTS: AUC: 66.14% on (wrongfully) unfiltered data, though the difference is 
# expected given: 1) different estimation method;
# and 2) cannot accommodate non-zero starting times without filtering.
#   When filtered: AUC is NA, though the approach is misspecified and flawed.
### CONCLUSION: Discard and disavow the use of timeROC at present.



# ----- Package: risksetROC
# Using risksetROC from Heagerty2005 with various kernel-choices though under the ID-approach

# - Calculate AUC at given prediction time for correctly-fitted Cox model | Cox-kernel of S(t)
# NOTE1: This approach assumes the Incident/Dynamic (ID) ROC-curve is desired, which breaks parity with other methods
# I.e., TPR is estimated given that the raw end point = t exactly
risksetROC::risksetROC(Stime=dat$End, status=dat$Event_Ind, entry=dat$Start,
                       method="Cox", span=0.05, #only used when using the NNE-method with "LocalCox" or "Schoenfeld" (presumably different kernels)
                       predict.time=predictTime, marker=round(predict(coxExample, type="risk"),2))
### RESULTS: AUC: 63.06% up to t

# - Calculate AUC at given prediction time for correctly-fitted Cox model | risksetROC (LocalCox-kernel)
# NOTE1: This approach assumes the Incident/Dynamic (ID) ROC-curve is desired, which breaks parity with other methods
# I.e., TPR is estimated given that the raw end point = t exactly
risksetROC::risksetROC(Stime=dat$End, status=dat$Event_Ind, entry=dat$Start,
                       method="LocalCox", span=0.05, #only used when using the NNE-method with "LocalCox" or "Schoenfeld" (presumably different kernels)
                       predict.time=predictTime, marker=round(predict(coxExample, type="risk"),2))
### RESULTS: AUC: 94.16% at t

# - Calculate AUC at given prediction time for correctly-fitted Cox model | risksetROC (Cox-kernel)
# NOTE1: This approach assumes the Incident/Dynamic (ID) ROC-curve is desired, which breaks parity with other methods
# I.e., TPR is estimated given that the raw end point = t exactly
risksetROC::risksetROC(Stime=dat$End, status=dat$Event_Ind, entry=dat$Start,
                       method="Schoenfeld", span=0.05, #only used when using the NNE-method with "LocalCox" or "Schoenfeld" (presumably different kernels)
                       predict.time=predictTime, marker=round(predict(coxExample, type="risk"),2))
### RESULTS: AUC: 98.7% at t
### CONCLUSION: The AUC-results vary widely, which definitely hints at "functional misspecification", though it is
#               obvious how exactly, despite consulting the developer notes extensively.



# ----- Package: tROCkit() | custom "package"/function
# Using custom tROC()-function from script 0b(iii) under the CD-approach with an NN-estimator

# - Calculate AUC up to given prediction time for correctly-fitted Cox model, assuming independence | NNE-kernel for S(t)
# NOTE: Uses the superior Nearest Neighbour Estimator (NNE) method for S(t)
# NOTE2: Foregoes the ID-field for comparative purposes with survivalROC()
tROC1 <- tROC(datGiven=dat, cox=coxExample, month_End=predictTime, estMethod="NN-0/1", numDigits=2, 
              fld_Event="Event_Ind", eventVal=1, fld_StartTime="Start", fld_EndTime="End",
              graphName="coxExample1_cgd_independence", genFigPath=genFigPath)
tROC1$AUC; tROC1$ROC_graph
### RESULTS: AUC: 57.69% up to t, which is very close to that of survivalROC
# In fact, we discovered a small indexing error in survivalROC() which explain this small discrepancy exactly


# - Calculate AUC up to given prediction time for correctly-fitted Cox model, assuming independence | NNE-kernel for S(t)
# NOTE: Uses the superior Nearest Neighbour Estimator (NNE) method for S(t)
# NOTE2: Enforces the ID-field.
tROC2 <- tROC(datGiven=dat, cox=coxExample, month_End=predictTime, estMethod="NN-0/1", numDigits=2, 
              fld_ID="ID", fld_Event="Event_Ind", eventVal=1, fld_StartTime="Start", fld_EndTime="End",
              graphName="coxExample1_cgd_dependence", genFigPath=genFigPath)
tROC2$AUC; tROC2$ROC_graph
### RESULTS: AUC: 59.36% up to t


# - Multi-threaded calculation of the  AUC up to given prediction time for correctly-fitted Cox model, 
# assuming independence | NNE-kernel for S(t)
# NOTE: Uses the superior Nearest Neighbour Estimator (NNE) method for S(t)
# NOTE2: Enforces the ID-field.
tROC2.multi <- tROC.multi(datGiven=dat, cox=coxExample, month_End=predictTime, estMethod="NN-0/1", numDigits=2, 
              fld_ID="ID", fld_Event="Event_Ind", eventVal=1, fld_StartTime="Start", fld_EndTime="End",
              graphName="coxExample1_cgd_dependence", genFigPath=genFigPath, numThreads=6)
tROC2.multi$AUC; tROC2.multi$ROC_graph
### RESULTS: AUC: 59.36% up to t



### AB: Need to rewire the fields here, though I need an actual dataset to do that.
tROC(datGiven=datCredit_valid_TFD, cox=coxDelinq, month_End=12, sLambda=0.05, estMethod="NN-0/1", numDigits=2, 
     fld_ID="PerfSpell_Key", fld_Event="Default_Ind", eventVal=1, fld_StartTime="Start", fld_EndTime="Stop",
     graphName="DefaultSurvModel-Cox", genFigPath=paste0(genFigPath, "TFD/tdROC/"))



# - cleanup
suppressWarnings(rm(gsurv1c_a, gsurv1c_d, haz_dat, kmExample, kmExample_survFitSummary, coxExample, coxExample2))
  
  