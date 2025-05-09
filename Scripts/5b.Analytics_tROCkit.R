# ========================= TIME-DEPENDENT ROC-ANALYSIS ==========================
# Compare various functions from various packages in conducting time-dependent 
# ROC-analyses on the same fitted Cox regression model, having used the 
# prepared credit data
# --------------------------------------------------------------------------------
# PROJECT TITLE: Default Survival Modelling
# SCRIPT AUTHOR(S): Dr Arno Botha (AB)
# --------------------------------------------------------------------------------
# -- Script dependencies:
#   - 0.Setup.R
#   - 1.Data_Import.R
#   - 2a.Data_Prepare_Credit_Basic.R
#   - 2b.Data_Prepare_Credit_Advanced.R
#   - 2c.Data_Prepare_Credit_Advanced2.R
#   - 2d.Data_Enrich.R
#   - 2e.Data_Prepare_Macro.R
#   - 2f.Data_Fusion1.R
#   - 3a(i).Data_Transform.R
#   - 3c(i).Data_Fusion2_TFD.R
#   - 5a(i).InputSpace_TFD.R
#   - 5a(ii).InputSpace_AG.R
#   - 5a(iii).InputSpace_PWPST.R

# -- Inputs:
#   - datCredit_train_TFD | Prepared from script 3b
#   - datCredit_valid_TFD | Prepared from script 3b
#
# -- Outputs:
#   - <Analytics> | tROC-graphs
# ================================================================================





# ----------------- 1. Load data and conduct preliminary Kaplan-Meier (KM) analyses
# NOTE: Conducting KM-analyses helps identify interesting prediction times for the
# eventual ROC-analyses


# ------ Time to first Default (TFD) definition
# - Confirm prepared datasets are loaded into memory
if (!exists('datCredit_train_TFD')) unpack.ffdf(paste0(genPath,"creditdata_train_TFD"), tempPath);gc()
if (!exists('datCredit_valid_TFD')) unpack.ffdf(paste0(genPath,"creditdata_valid_TFD"), tempPath);gc()

# - Fit Kaplan-Meier (KM) nonparametric (and "empty-of-covariates") model
# Compute Kaplan-Meier survival estimates (product-limit) for main-event | Spell-level with right-censoring & left-truncation
# All competing events preclude the main event from happening and are therefore considered as censored
# ID is set as the spell key, with no stratification
km_TFD <- survfit(Surv(time=Start, time2=End, event=Default_Ind==1,type="counting") ~ 1, 
                     id=LoanID, data=datCredit_train_TFD)
summary(km_TFD)$table # overall summary statistics
### RESULTS: 62k events, with no median survival probability
(km_TFD_survFitSummary <- surv_summary(km_TFD)) # Survival table

# -- Calculate various summary statistics
# - Chosen percentile of survival times, used later as a prediction time, e.g., 90%
# NOTE: Median would also work, though sample sizes typically dwindle at that point, which can complicate analyses
cat(paste0("Kaplan-Meier: Survival time at the xth percentile of the distribution of unique survival times: \n\t",
           percTimepoint <- min(km_TFD$time[km_TFD$surv <= 0.90],na.rm=T)))
### RESULTS: Spell Period @ 90% mean survival prob: 36 months
# - Calculate the Truncated Mean (or "Restricted Mean Survival Time" [RMST]) of survival times
# NOTE: The classical mean survival time is only defined when S(t) can reach 0
# NOTE2: Calculating RMST involves integrating the KM-curve up to the last (hence truncated) observed unique event time
#         I.e., sum the survival probabilities, as weighted by the time difference between two successive survival times
cat(paste0("Kaplan-Meier: Truncated/Restricted Mean Survival Time: \n\t", 
           sum(km_TFD$surv*diff(c(0,km_TFD$time)), na.rm=T)))
### RESULTS: RMST: 156.9 months
# - Retrieve the mean survival probability at the previously found percentile-basd prediction point
cat(paste0("Kaplan-Meier: Retrieved mean survival probability at a given time point (t=", percTimepoint, "): \n\t",
           percent(km_TFD$surv[max(which(km_TFD$time <= percTimepoint))], accuracy=0.01)))
### RESULTS: Mean survival prob at prediction point: 89.92%
# - Retrieve the mean survival probability at the last unique event time as prediction time, i.e., 
# the "pooled mean survival probability"
cat(paste0("Kaplan-Meier: Retrieved mean survival probability at last unique event point as time point (t=", 
           max(km_TFD$time), "): \n\t", percent(km_TFD$surv[max(which(km_TFD$time <= max(km_TFD$time)))], accuracy=0.01)))
### RESULTS: 67.74% (Also last row in [km_TFD_survFitSummary])




# ------ Andersen-Gill (AG) Spell-time definition
# - Confirm prepared datasets are loaded into memory
if (!exists('datCredit_train_AG')) unpack.ffdf(paste0(genPath,"creditdata_train_AG"), tempPath);gc()
if (!exists('datCredit_valid_AG')) unpack.ffdf(paste0(genPath,"creditdata_valid_AG"), tempPath);gc()

# - Fit Kaplan-Meier (KM) nonparametric (and "empty-of-covariates") model
# Compute Kaplan-Meier survival estimates (product-limit) for main-event | Spell-level with right-censoring & left-truncation
# All competing events preclude the main event from happening and are therefore considered as censored
# ID is set as the spell key, with no stratification
km_AG <- survfit(Surv(time=Start, time2=End, event=Default_Ind==1,type="counting") ~ 1, 
                    id=PerfSpell_Key, data=datCredit_train_AG)
summary(km_AG)$table # overall summary statistics
### RESULTS: 62k events, with no median survival probability
(km_AG_survFitSummary <- surv_summary(km_AG)) # Survival table



# ------ Prentice-Williams-Peterson (PWP) Spell-time definition
# - Confirm prepared datasets are loaded into memory
if (!exists('datCredit_train_PWPST')) unpack.ffdf(paste0(genPath,"creditdata_train_PWPST"), tempPath);gc()
if (!exists('datCredit_valid_PWPST')) unpack.ffdf(paste0(genPath,"creditdata_valid_PWPST"), tempPath);gc()

# - Fit Kaplan-Meier (KM) nonparametric (and "empty-of-covariates") model
# Compute Kaplan-Meier survival estimates (product-limit) for main-event | Spell-level with right-censoring & left-truncation
# All competing events preclude the main event from happening and are therefore considered as censored
# ID is set as the spell key, with no stratification
km_PWPST <- survfit(Surv(time=Start, time2=End, event=Default_Ind==1,type="counting") ~ 1, 
                  id=PerfSpell_Key, data=datCredit_train_PWPST)
summary(km_PWPST)$table # overall summary statistics
### RESULTS: 62k events, with no median survival probability
(km_PWPST_survFitSummary <- surv_summary(km_PWPST)) # Survival table





# ----------------- 2. Fit a Cox regression model on the resampled prepared data

# ------ Time to first Default (TFD) definition
# - Initialize variables
vecVars_TFD <- c("g0_Delinq_SD_4", "Arrears", "g0_Delinq_Ave",      
                "slc_acct_arr_dir_3_Change_Ind", "slc_acct_roll_ever_24_imputed_mean",
                "slc_acct_pre_lim_perc_imputed_med", "M_DTI_Growth",
                "M_RealIncome_Growth","pmnt_method_grp",
                "InterestRate_Nom", "Principal")

# - Fit a Cox Proportional Hazards model with time-varying covariates, and clustered observations
# NOTE: Assume dependence (by specifying ID-field) amongst certain observations clustered around ID-values
cox_TFD <- coxph(as.formula(paste0("Surv(Start,End,Default_Ind) ~ ",
                                   paste(vecVars_TFD,collapse=" + "))), id=LoanID, ties="efron", datCredit_train_TFD)
summary(cox_TFD); AIC(cox_TFD); concordance(cox_TFD)



# ------ Time to first Default (TFD) definition
# - Initialize variables
vecVars_AG <- c("PerfSpell_Num_binned","g0_Delinq_SD_4", "Arrears", "g0_Delinq_Ave", "TimeInDelinqState_Lag_1",
                "slc_acct_arr_dir_3_Change_Ind", "slc_acct_roll_ever_24_imputed_mean","LN_TPE",
                "slc_acct_pre_lim_perc_imputed_med","pmnt_method_grp","M_Inflation_Growth",
                "InterestRate_Nom", "BalanceToPrincipal","AgeToTerm_Aggr_Mean","M_DTI_Growth_9")

# - Fit a Cox Proportional Hazards model with time-varying covariates, and clustered observations
# NOTE: Assume dependence (by specifying ID-field) amongst certain observations clustered around ID-values
cox_AG <- coxph(as.formula(paste0("Surv(Start,End,Default_Ind) ~ ",
                                   paste(vecVars_AG,collapse=" + "))), 
                id=PerfSpell_Key, ties="efron", datCredit_train_AG)
summary(cox_AG); AIC(cox_AG); concordance(cox_AG)



# ------ Prentice-Williams-Peterson (PWP) Spell-time definition
# - Initialize variables
vecVars_PWPST <- c( # Delinquency-themed
  "g0_Delinq_SD_4", "slc_acct_roll_ever_24_imputed_mean", "g0_Delinq_Ave", "Arrears", "PerfSpell_Num",
  # Portfolio-level variables
  "AgeToTerm_Aggr_Mean",
  # Loan-level variables
  "BalanceToPrincipal", "pmnt_method_grp", "InterestRate_Nom", "slc_acct_arr_dir_3_Change_Ind",
  # Macroeconomic variables
  "M_DTI_Growth_9", "M_Inflation_Growth_6", "M_Repo_Rate_6")

# # Fit a Cox Proportional Hazards model with time-varying covariates, and clustered observations
# # NOTE: Assume dependence (by specifying ID-field) amongst certain observations clustered around ID-values
cox_PWPST <- coxph(as.formula(paste0("Surv(Start,End,Default_Ind) ~ ", paste(vecVars_PWPST,collapse=" + "), 
                                     " + strata(PerfSpell_Num_binned)")),
                   id=PerfSpell_Key, datCredit_train_PWPST, ties="efron")
summary(cox_PWPST); AIC(cox_PWPST); concordance(cox_PWPST)





# ----------------- 3. Conduct & compare ROC-analyses at fixed prediction time
# NOTE: ROC-analyses are conducted using the same fitted Cox model at a specific prediction time t


# ------ Time to first Default (TFD) definition

# --- Package: survivalROC()
# Using survivalROC with NNE-estimator of S(t) under the CD-approach | NNE-kernel for S(t)

# - Decide prediction point t at which we would like to conduct a time-dependent ROC-analysis
# assuming the cumulative/dynamic (CD) context, i.e., that TPR is calculated given that 
# raw event times T_{ijt} are less than or equal to the given prediction time t
predictTime <- 12 # 90% percentile survival time (203); median: 338

# - Calculate AUC up to given prediction time for correctly-fitted Cox model
# NOTE: Uses the superior Nearest Neighbour Estimator (NNE) method for S(t)
ptm <- proc.time() #IGNORE: for computation time calculation;
survivalROC::survivalROC(Stime=datCredit_train_TFD$End, status=datCredit_train_TFD$Default_Ind, entry=datCredit_train_TFD$Start, 
                         method = "NNE",  span=0.05, predict.time=predictTime,
                         marker=round(predict(cox_TFD, type="lp"),2))
proc.time() - ptm
### RESULTS: AUC: 94% % up to t, achieved in 8781.23 sec (146.35 mins)



# --- Package: tROCkit() | custom "package"/function
# NOTE: Using custom tROC()-function from script 0b(iii) under the CD-approach with an NN-estimator and 0/1-kernel

# -- Multi-threaded calculation of the # AUC from given start up to given prediction time 3 in following the CD-approach
# NOTE: Uses the superior Nearest Neighbour Estimator (NNE) method for S(t) with a 0/1-kernelNNE-kernel for S(t)
# NOTE2: Assume dependence (by specifying ID-field) amongst certain observations clustered around ID-values
ptm <- proc.time() #IGNORE: for computation time calculation;
predictTime <- 3
objROC1_TFD <- tROC.multi(datGiven=datCredit_valid_TFD, cox=cox_TFD, month_End=predictTime, sLambda=0.05, estMethod="NN-0/1", numDigits=4, 
                          fld_ID="PerfSpell_Key", fld_Event="PerfSpell_Event", eventVal=1, fld_StartTime="Start", fld_EndTime="End",
                          graphName="DefaultSurvModel-Cox-TFD-ROC_Depedendence", genFigPath=paste0(genFigPath, "TFD/"), 
                          caseStudyName=paste0("TFD_", predictTime), numThreads=12, logPath=genPath)
objROC1_TFD$AUC; objROC1_TFD$ROC_graph
proc.time() - ptm
### RESULTS: AUC up to t: 96.38%, achieved in 746.16  secs ( 12.4 mins)


# -- Calculate AUC from given start up to given prediction time in following the CD-approach
# NOTE: Uses the superior Nearest Neighbour Estimator (NNE) method for S(t) with a 0/1-kernelNNE-kernel for S(t)
# NOTE2: Assume dependence (by specifying ID-field) amongst certain observations clustered around ID-values
# NOTE3: This iteration is obviously slow, but kept for quantifying the improvement in runtime speed.
# ptm <- proc.time() #IGNORE: for computation time calculation;
# objROC2_TFD <- tROC(datGiven=datCredit_valid_TFD, cox=cox_TFD, month_End=12, sLambda=0.05, estMethod="NN-0/1", numDigits=2, 
#                     fld_ID="PerfSpell_Key", fld_Event="PerfSpell_Event", eventVal=1, fld_StartTime="Start", fld_EndTime="End",
#                     graphName="DefaultSurvModel-Cox1_Depedendence", genFigPath=paste0(genFigPath, "TFD/"))
# objROC2_TFD$AUC; objROC2_TFD$ROC_graph
# proc.time() - ptm
### RESULTS: AUC up to t: 91.21%, achieved in 5633.98 secs (94 mins)


# -- Multi-threaded calculation of the # AUC from given start up to given prediction time 12 in following the CD-approach
# NOTE: Uses the superior Nearest Neighbour Estimator (NNE) method for S(t) with a 0/1-kernelNNE-kernel for S(t)
# NOTE2: Assume dependence (by specifying ID-field) amongst certain observations clustered around ID-values
ptm <- proc.time() #IGNORE: for computation time calculation;
predictTime <- 12
objROC2_TFD <- tROC.multi(datGiven=datCredit_valid_TFD, cox=cox_TFD, month_End=predictTime, sLambda=0.05, estMethod="NN-0/1", numDigits=4,
                      fld_ID="PerfSpell_Key", fld_Event="PerfSpell_Event", eventVal=1, fld_StartTime="Start", fld_EndTime="End",
                      graphName="DefaultSurvModel-Cox-TFD-ROC_Depedendence", genFigPath=paste0(genFigPath, "TFD/"),
                      caseStudyName=paste0("TFD_", predictTime), numThreads=12, logPath=genPath)
objROC2_TFD$AUC; objROC2_TFD$ROC_graph
proc.time() - ptm
### RESULTS: AUC up to t: 96.39%, achieved in 1490.08   secs (24.8 mins)


# -- Multi-threaded calculation of the # AUC from given start up to given prediction time 36 in following the CD-approach
# NOTE: Uses the superior Nearest Neighbour Estimator (NNE) method for S(t) with a 0/1-kernelNNE-kernel for S(t)
# NOTE2: Assume dependence (by specifying ID-field) amongst certain observations clustered around ID-values
ptm <- proc.time() #IGNORE: for computation time calculation;
predictTime <- 24
objROC3_TFD <- tROC.multi(datGiven=datCredit_valid_TFD, cox=cox_TFD, month_End=predictTime, sLambda=0.05, estMethod="NN-0/1", numDigits=4, 
                          fld_ID="PerfSpell_Key", fld_Event="PerfSpell_Event", eventVal=1, fld_StartTime="Start", fld_EndTime="End",
                          graphName="DefaultSurvModel-Cox-TFD-ROC_Depedendence", genFigPath=paste0(genFigPath, "TFD/"), 
                          caseStudyName=paste0("TFD_", predictTime), numThreads=12, logPath=genPath)
objROC3_TFD$AUC; objROC3_TFD$ROC_graph
proc.time() - ptm
### RESULTS: AUC up to t: 96.34%, achieved in 2362.68 secs (39.4 mins)


# -- Multi-threaded calculation of the # AUC from given start up to given prediction time 36 in following the CD-approach
# NOTE: Uses the superior Nearest Neighbour Estimator (NNE) method for S(t) with a 0/1-kernelNNE-kernel for S(t)
# NOTE2: Assume dependence (by specifying ID-field) amongst certain observations clustered around ID-values
ptm <- proc.time() #IGNORE: for computation time calculation;
predictTime <- 36
objROC4_TFD <- tROC.multi(datGiven=datCredit_valid_TFD, cox=cox_TFD, month_End=predictTime, sLambda=0.05, estMethod="NN-0/1", numDigits=4, 
                          fld_ID="PerfSpell_Key", fld_Event="PerfSpell_Event", eventVal=1, fld_StartTime="Start", fld_EndTime="End",
                          graphName="DefaultSurvModel-Cox-TFD-ROC_Depedendence", genFigPath=paste0(genFigPath, "TFD/"), 
                          caseStudyName=paste0("TFD_", predictTime), numThreads=12, logPath=genPath)
objROC4_TFD$AUC; objROC4_TFD$ROC_graph
proc.time() - ptm
### RESULTS: AUC up to t: 96.40%, achieved in 2362.68 secs (39.4 mins)


# -- Store experimental objects | Memory optimisation
pack.ffdf(paste0(genPath,"DefaultSurvModel-Cox-TFD-ROC_Depedendence_03"), objROC1_TFD);
pack.ffdf(paste0(genPath,"DefaultSurvModel-Cox-TFD-ROC_Depedendence_12"), objROC2_TFD);
pack.ffdf(paste0(genPath,"DefaultSurvModel-Cox-TFD-ROC_Depedendence_24"), objROC3_TFD);
pack.ffdf(paste0(genPath,"DefaultSurvModel-Cox-TFD-ROC_Depedendence_36"), objROC4_TFD);



# ------ Andersen-Gill (AG) Spell-time definition

# --- Package: tROCkit() | custom "package"/function
# NOTE: Using custom tROC()-function from script 0b(iii) under the CD-approach with an NN-estimator and 0/1-kernel


# -- Multi-threaded calculation of the # AUC from given start up to given prediction time 3 in following the CD-approach
# NOTE: Uses the superior Nearest Neighbour Estimator (NNE) method for S(t) with a 0/1-kernelNNE-kernel for S(t)
# NOTE2: Assume dependence (by specifying ID-field) amongst certain observations clustered around ID-values
ptm <- proc.time() #IGNORE: for computation time calculation;
predictTime <- 3
objROC1_AG <- tROC.multi(datGiven=datCredit_valid_AG, cox=cox_AG, month_End=predictTime, sLambda=0.05, estMethod="NN-0/1", numDigits=4, 
                            fld_ID="PerfSpell_Key", fld_Event="PerfSpell_Event", eventVal=1, fld_StartTime="Start", fld_EndTime="End",
                            graphName="DefaultSurvModel-Cox-AG-ROC_Depedendence", genFigPath=paste0(genFigPath, "AG/"), 
                            caseStudyName=paste0("AG_", predictTime), numThreads=12, logPath=genPath)
objROC1_AG$AUC; objROC1_AG$ROC_graph
proc.time() - ptm
### RESULTS: AUC up to t: 82.8%, achieved in 686.18  secs ( 11.4 mins)


# -- Multi-threaded calculation of the # AUC from given start up to given prediction time 12 in following the CD-approach
# NOTE: Uses the superior Nearest Neighbour Estimator (NNE) method for S(t) with a 0/1-kernelNNE-kernel for S(t)
# NOTE2: Assume dependence (by specifying ID-field) amongst certain observations clustered around ID-values
ptm <- proc.time() #IGNORE: for computation time calculation;
predictTime <- 12
objROC2_AG <- tROC.multi(datGiven=datCredit_valid_AG, cox=cox_AG, month_End=predictTime, sLambda=0.05, estMethod="NN-0/1", numDigits=4,
                            fld_ID="PerfSpell_Key", fld_Event="PerfSpell_Event", eventVal=1, fld_StartTime="Start", fld_EndTime="End",
                            graphName="DefaultSurvModel-Cox-AG-ROC_Depedendence", genFigPath=paste0(genFigPath, "AG/"),
                            caseStudyName=paste0("AG_", predictTime), numThreads=12, logPath=genPath)
objROC2_AG$AUC; objROC2_AG$ROC_graph
proc.time() - ptm
### RESULTS: AUC up to t: 95.70%, achieved in 1418.35  secs (23.6 mins)


# -- Multi-threaded calculation of the # AUC from given start up to given prediction time 36 in following the CD-approach
# NOTE: Uses the superior Nearest Neighbour Estimator (NNE) method for S(t) with a 0/1-kernelNNE-kernel for S(t)
# NOTE2: Assume dependence (by specifying ID-field) amongst certain observations clustered around ID-values
ptm <- proc.time() #IGNORE: for computation time calculation;
predictTime <- 24
objROC3_AG <- tROC.multi(datGiven=datCredit_valid_AG, cox=cox_AG, month_End=predictTime, sLambda=0.05, estMethod="NN-0/1", numDigits=4, 
                            fld_ID="PerfSpell_Key", fld_Event="PerfSpell_Event", eventVal=1, fld_StartTime="Start", fld_EndTime="End",
                            graphName="DefaultSurvModel-Cox-AG-ROC_Depedendence", genFigPath=paste0(genFigPath, "AG/"), 
                            caseStudyName=paste0("AG_", predictTime), numThreads=12, logPath=genPath)
objROC3_AG$AUC; objROC3_AG$ROC_graph
proc.time() - ptm
### RESULTS: AUC up to t: 95.91%, achieved in 2308.32 secs (38.5 mins)


# -- Multi-threaded calculation of the # AUC from given start up to given prediction time 36 in following the CD-approach
# NOTE: Uses the superior Nearest Neighbour Estimator (NNE) method for S(t) with a 0/1-kernelNNE-kernel for S(t)
# NOTE2: Assume dependence (by specifying ID-field) amongst certain observations clustered around ID-values
ptm <- proc.time() #IGNORE: for computation time calculation;
predictTime <- 36
objROC4_AG <- tROC.multi(datGiven=datCredit_valid_AG, cox=cox_AG, month_End=predictTime, sLambda=0.05, estMethod="NN-0/1", numDigits=4, 
                            fld_ID="PerfSpell_Key", fld_Event="PerfSpell_Event", eventVal=1, fld_StartTime="Start", fld_EndTime="End",
                            graphName="DefaultSurvModel-Cox-AG-ROC_Depedendence", genFigPath=paste0(genFigPath, "AG/"), 
                            caseStudyName=paste0("AG_", predictTime), numThreads=12, logPath=genPath)
objROC4_AG$AUC; objROC4_AG$ROC_graph
proc.time() - ptm
### RESULTS: AUC up to t: 95.98%, achieved in 3102.89 secs (51.7 mins)


# -- Store experimental objects | Memory optimisation
pack.ffdf(paste0(genPath,"DefaultSurvModel-Cox-AG-ROC_Depedendence_03"), objROC1_AG);
pack.ffdf(paste0(genPath,"DefaultSurvModel-Cox-AG-ROC_Depedendence_12"), objROC2_AG);
pack.ffdf(paste0(genPath,"DefaultSurvModel-Cox-AG-ROC_Depedendence_24"), objROC3_AG);
pack.ffdf(paste0(genPath,"DefaultSurvModel-Cox-AG-ROC_Depedendence_36"), objROC4_AG);




# ------ Prentice-Williams-Peterson (PWP) Spell-time definition

# --- Package: tROCkit() | custom "package"/function
# NOTE: Using custom tROC()-function from script 0b(iii) under the CD-approach with an NN-estimator and 0/1-kernel


# -- Multi-threaded calculation of the # AUC from given start up to given prediction time 3 in following the CD-approach
# NOTE: Uses the superior Nearest Neighbour Estimator (NNE) method for S(t) with a 0/1-kernelNNE-kernel for S(t)
# NOTE2: Assume dependence (by specifying ID-field) amongst certain observations clustered around ID-values
ptm <- proc.time() #IGNORE: for computation time calculation;
predictTime <- 3
objROC1_PWPST <- tROC.multi(datGiven=datCredit_valid_PWPST, cox=cox_PWPST, month_End=predictTime, sLambda=0.05, estMethod="NN-0/1", numDigits=4, 
                          fld_ID="PerfSpell_Key", fld_Event="PerfSpell_Event", eventVal=1, fld_StartTime="Start", fld_EndTime="End",
                          graphName="DefaultSurvModel-Cox-PWPST-ROC_Depedendence", genFigPath=paste0(genFigPath, "PWP ST/"), 
                          caseStudyName=paste0("PWPST_", predictTime), numThreads=12, logPath=genPath)
objROC1_PWPST$AUC; objROC1_PWPST$ROC_graph
proc.time() - ptm
### RESULTS: AUC up to t: 96.35%, achieved in 712.29 secs ( 11.9 mins)


# -- Multi-threaded calculation of the # AUC from given start up to given prediction time 12 in following the CD-approach
# NOTE: Uses the superior Nearest Neighbour Estimator (NNE) method for S(t) with a 0/1-kernelNNE-kernel for S(t)
# NOTE2: Assume dependence (by specifying ID-field) amongst certain observations clustered around ID-values
ptm <- proc.time() #IGNORE: for computation time calculation;
predictTime <- 12
objROC2_PWPST <- tROC.multi(datGiven=datCredit_valid_PWPST, cox=cox_PWPST, month_End=predictTime, sLambda=0.05, estMethod="NN-0/1", numDigits=4,
                          fld_ID="PerfSpell_Key", fld_Event="PerfSpell_Event", eventVal=1, fld_StartTime="Start", fld_EndTime="End",
                          graphName="DefaultSurvModel-Cox-PWPST-ROC_Depedendence", genFigPath=paste0(genFigPath, "PWP ST/"),
                          caseStudyName=paste0("PWPST_", predictTime), numThreads=12, logPath=genPath)
objROC2_PWPST$AUC; objROC2_PWPST$ROC_graph
proc.time() - ptm
### RESULTS: AUC up to t: 96.15%, achieved in 1370.21 secs (22.8 mins)


# -- Multi-threaded calculation of the # AUC from given start up to given prediction time 36 in following the CD-approach
# NOTE: Uses the superior Nearest Neighbour Estimator (NNE) method for S(t) with a 0/1-kernelNNE-kernel for S(t)
# NOTE2: Assume dependence (by specifying ID-field) amongst certain observations clustered around ID-values
ptm <- proc.time() #IGNORE: for computation time calculation;
predictTime <- 24
objROC3_PWPST <- tROC.multi(datGiven=datCredit_valid_PWPST, cox=cox_PWPST, month_End=predictTime, sLambda=0.05, estMethod="NN-0/1", numDigits=4, 
                          fld_ID="PerfSpell_Key", fld_Event="PerfSpell_Event", eventVal=1, fld_StartTime="Start", fld_EndTime="End",
                          graphName="DefaultSurvModel-Cox-PWPST-ROC_Depedendence", genFigPath=paste0(genFigPath, "PWP ST/"), 
                          caseStudyName=paste0("PWPST_", predictTime), numThreads=12, logPath=genPath)
objROC3_PWPST$AUC; objROC3_PWPST$ROC_graph
proc.time() - ptm
### RESULTS: AUC up to t: 96.18%, achieved in 2264.65 secs (37.7 mins)


# -- Multi-threaded calculation of the # AUC from given start up to given prediction time 36 in following the CD-approach
# NOTE: Uses the superior Nearest Neighbour Estimator (NNE) method for S(t) with a 0/1-kernelNNE-kernel for S(t)
# NOTE2: Assume dependence (by specifying ID-field) amongst certain observations clustered around ID-values
ptm <- proc.time() #IGNORE: for computation time calculation;
predictTime <- 36
objROC4_PWPST <- tROC.multi(datGiven=datCredit_valid_PWPST, cox=cox_PWPST, month_End=predictTime, sLambda=0.05, estMethod="NN-0/1", numDigits=4, 
                          fld_ID="PerfSpell_Key", fld_Event="PerfSpell_Event", eventVal=1, fld_StartTime="Start", fld_EndTime="End",
                          graphName="DefaultSurvModel-Cox-PWPST-ROC_Depedendence", genFigPath=paste0(genFigPath, "PWP ST/"), 
                          caseStudyName=paste0("PWPST_", predictTime), numThreads=12, logPath=genPath)
objROC4_PWPST$AUC; objROC4_PWPST$ROC_graph
proc.time() - ptm
### RESULTS: AUC up to t: 95.96%, achieved in 3166.55 secs (52.8 mins)


# -- Store experimental objects | Memory optimisation
pack.ffdf(paste0(genPath,"DefaultSurvModel-Cox-PWPST-ROC_Depedendence_03"), objROC1_PWPST);
pack.ffdf(paste0(genPath,"DefaultSurvModel-Cox-PWPST-ROC_Depedendence_12"), objROC2_PWPST);
pack.ffdf(paste0(genPath,"DefaultSurvModel-Cox-PWPST-ROC_Depedendence_24"), objROC3_PWPST);
pack.ffdf(paste0(genPath,"DefaultSurvModel-Cox-PWPST-ROC_Depedendence_36"), objROC4_PWPST);





# ----------------- 4. Create combined ROC-graph across multiple prediction times

# ------ Time to first Default (TFD) definition

# - Ensure required objects exist in memory
if (!exists('objROC1_TFD')) unpack.ffdf(paste0(genPath,"DefaultSurvModel-Cox-TFD-ROC_Depedendence_03"), tempPath);gc()
if (!exists('objROC2_TFD')) unpack.ffdf(paste0(genPath,"DefaultSurvModel-Cox-TFD-ROC_Depedendence_12"), tempPath);gc()
if (!exists('objROC3_TFD')) unpack.ffdf(paste0(genPath,"DefaultSurvModel-Cox-TFD-ROC_Depedendence_24"), tempPath);gc()
if (!exists('objROC4_TFD')) unpack.ffdf(paste0(genPath,"DefaultSurvModel-Cox-TFD-ROC_Depedendence_36"), tempPath);gc()

# - Set ROC-parameters and initialize data structures
vecPercTimepoint <- c(3,12,24,36)
vecTROC <- list(objROC1_TFD, objROC2_TFD, objROC3_TFD, objROC4_TFD)
vLabels <- vector("list", length=length(vecPercTimepoint))

# -- Create a combined data object for plotting purposes
for (i in 1:length(vecPercTimepoint)) {
  # i <-1 # testing condition
  
  # datGraph <- data.frame(x = vFPR[-(nThresh+1)], y=vTPR[-1])
  
  # - Create a data object for the current prediction time
  if (i == 1) {
    datGraph <- data.table(PredictTime=paste0(letters[i], "_", vecPercTimepoint[i]), Threshold=vecTROC[[i]]$Thresholds, 
                           x=vecTROC[[i]]$FPR, y=vecTROC[[i]]$TPR)
    
  } else {
    datGraph <- rbind(datGraph, 
                      data.table(PredictTime= paste0(letters[i], "_", vecPercTimepoint[i]), Threshold=vecTROC[[i]]$Thresholds, 
                                 x=vecTROC[[i]]$FPR, y=vecTROC[[i]]$TPR))
  }
  vLabels[[i]] <- bquote("Prediction time "*italic(t)==.(vecPercTimepoint[i])*"; AUC: "*.(percent(vecTROC[[i]]$AUC, accuracy=0.01)))
}


# -- Graph a combined ROC-graph across prediction times t
# - Aesthetic parameters
datGraph[, FacetLabel := "Time to First Default (TFD) model"]
vCol <- brewer.pal(8,"Set1")
vLabels_F <- setNames(vLabels, paste0(letters[1:length(vecPercTimepoint)],"_", vecPercTimepoint))
chosenFont <- "Cambria"

# - Create ROC-graph
(gg <- ggplot(datGraph, aes(x=x,y=y,group=PredictTime)) + theme_minimal() + 
    theme(text = element_text(family=chosenFont), legend.position="inside",
          strip.background=element_rect(fill="snow2", colour="snow2"),
          strip.text=element_text(size=8, colour="gray50"), strip.text.y.right=element_text(angle=90),
          legend.position.inside = c(0.55,0.4),
          legend.background = element_rect(fill="snow2", color="black",
                                           linetype="solid", linewidth=0.1)) +
    labs(x = bquote("False Positive Rate "*italic(F^"+")), y = 
           bquote("True Positive Rate "*italic(T^"+"))) + 
    # Add 45-degree line
    geom_segment(x = 0, y = 0, xend = 1, yend = 1, color = "grey", linewidth=0.2) +
    # Main line graph
    geom_step(aes(x=x, y=y, linetype=PredictTime, colour=PredictTime), linewidth=0.05) + 
    geom_point(aes(x=x, y=y, shape=PredictTime, colour=PredictTime), size=0.25) + 
    # Facets and scales
    facet_grid(FacetLabel ~ .) +  
    scale_color_manual(name=bquote("ROC"*(italic(t))), values=vCol, labels=vLabels) + 
    scale_linetype_discrete(name=bquote("ROC"*(italic(t))), labels=vLabels) + 
    scale_shape_discrete(name=bquote("ROC"*(italic(t))), labels=vLabels) + 
    scale_y_continuous(label=percent) + scale_x_continuous(label=percent))


# - Save graph
dpi <- 300
ggsave(gg, file=paste0(paste0(genFigPath, "TFD/DefaultSurvModel-Cox-TFD-CombinedROC_Depedendence.png")), 
       width=1200/dpi, height=1000/dpi, dpi=dpi, bg="white")




# ------ Andersen-Gill (AG) Total-time definition

# - Ensure required objects exist in memory
if (!exists('objROC1_AG')) unpack.ffdf(paste0(genPath,"DefaultSurvModel-Cox-AG-ROC_Depedendence_03"), tempPath);gc()
if (!exists('objROC2_AG')) unpack.ffdf(paste0(genPath,"DefaultSurvModel-Cox-AG-ROC_Depedendence_12"), tempPath);gc()
if (!exists('objROC3_AG')) unpack.ffdf(paste0(genPath,"DefaultSurvModel-Cox-AG-ROC_Depedendence_24"), tempPath);gc()
if (!exists('objROC4_AG')) unpack.ffdf(paste0(genPath,"DefaultSurvModel-Cox-AG-ROC_Depedendence_36"), tempPath);gc()

# - Set ROC-parameters and initialize data structures
vecPercTimepoint <- c(3,12,24,36)
vecTROC <- list(objROC1_AG, objROC2_AG, objROC3_AG, objROC4_AG)
vLabels <- vector("list", length=length(vecPercTimepoint))

# -- Create a combined data object for plotting purposes
for (i in 1:length(vecPercTimepoint)) {
  # i <-1 # testing condition
  
  # datGraph <- data.frame(x = vFPR[-(nThresh+1)], y=vTPR[-1])
  
  # - Create a data object for the current prediction time
  if (i == 1) {
    datGraph <- data.table(PredictTime=paste0(letters[i], "_", vecPercTimepoint[i]), Threshold=vecTROC[[i]]$Thresholds, 
                           x=vecTROC[[i]]$FPR, y=vecTROC[[i]]$TPR)
    
  } else {
    datGraph <- rbind(datGraph, 
                      data.table(PredictTime= paste0(letters[i], "_", vecPercTimepoint[i]), Threshold=vecTROC[[i]]$Thresholds, 
                                 x=vecTROC[[i]]$FPR, y=vecTROC[[i]]$TPR))
  }
  vLabels[[i]] <- bquote("Prediction time "*italic(t)==.(vecPercTimepoint[i])*"; AUC: "*.(percent(vecTROC[[i]]$AUC, accuracy=0.01)))
}


# -- Graph a combined ROC-graph across prediction times t
# - Aesthetic parameters
datGraph[, FacetLabel := "Andersen-Gill (AG) spell-time model"]
vCol <- brewer.pal(8,"Set1")
vLabels_F <- setNames(vLabels, paste0(letters[1:length(vecPercTimepoint)],"_", vecPercTimepoint))
chosenFont <- "Cambria"

# - Create ROC-graph
(gg <- ggplot(datGraph, aes(x=x,y=y,group=PredictTime)) + theme_minimal() + 
    theme(text = element_text(family=chosenFont), legend.position="inside", 
          strip.background=element_rect(fill="snow2", colour="snow2"),
          strip.text=element_text(size=8, colour="gray50"), strip.text.y.right=element_text(angle=90),
          legend.position.inside = c(0.55,0.45),
          legend.background = element_rect(fill="snow2", color="black",
                                           linetype="solid", linewidth=0.1)) +
    labs(x = bquote("False Positive Rate "*italic(F^"+")), y = 
           bquote("True Positive Rate "*italic(T^"+"))) + 
    # Add 45-degree line
    geom_segment(x = 0, y = 0, xend = 1, yend = 1, color = "grey", linewidth=0.2) +
    # Main line graph
    geom_step(aes(x=x, y=y, linetype=PredictTime, colour=PredictTime), linewidth=0.05) + 
    geom_point(aes(x=x, y=y, shape=PredictTime, colour=PredictTime), size=0.25) +
    # Facets and scales
    facet_grid(FacetLabel ~ .) +  
    scale_color_manual(name=bquote("ROC"*(italic(t))), values=vCol, labels=vLabels) + 
    scale_linetype_discrete(name=bquote("ROC"*(italic(t))), labels=vLabels) + 
    scale_shape_discrete(name=bquote("ROC"*(italic(t))), labels=vLabels) + 
    scale_y_continuous(label=percent) + scale_x_continuous(label=percent))


# - Save graph
dpi <- 300
ggsave(gg, file=paste0(paste0(genFigPath, "AG/DefaultSurvModel-Cox-AG-CombinedROC_Depedendence.png")), 
       width=1200/dpi, height=1000/dpi, dpi=dpi, bg="white")




# ------ Prentice-Williams-Peterson (PWP) Total-time definition

# - Ensure required objects exist in memory
if (!exists('objROC1_PWPST')) unpack.ffdf(paste0(genPath,"DefaultSurvModel-Cox-PWPST-ROC_Depedendence_03"), tempPath);gc()
if (!exists('objROC2_PWPST')) unpack.ffdf(paste0(genPath,"DefaultSurvModel-Cox-PWPST-ROC_Depedendence_12"), tempPath);gc()
if (!exists('objROC3_PWPST')) unpack.ffdf(paste0(genPath,"DefaultSurvModel-Cox-PWPST-ROC_Depedendence_24"), tempPath);gc()
if (!exists('objROC4_PWPST')) unpack.ffdf(paste0(genPath,"DefaultSurvModel-Cox-PWPST-ROC_Depedendence_36"), tempPath);gc()

# - Set ROC-parameters and initialize data structures
vecPercTimepoint <- c(3,12,24,36)
vecTROC <- list(objROC1_PWPST, objROC2_PWPST, objROC3_PWPST, objROC4_PWPST)
vLabels <- vector("list", length=length(vecPercTimepoint))

# -- Create a combined data object for plotting purposes
for (i in 1:length(vecPercTimepoint)) {
  # i <-1 # testing condition
  
  # datGraph <- data.frame(x = vFPR[-(nThresh+1)], y=vTPR[-1])
  
  # - Create a data object for the current prediction time
  if (i == 1) {
    datGraph <- data.table(PredictTime=paste0(letters[i], "_", vecPercTimepoint[i]), Threshold=vecTROC[[i]]$Thresholds, 
                           x=vecTROC[[i]]$FPR, y=vecTROC[[i]]$TPR)
    
  } else {
    datGraph <- rbind(datGraph, 
                      data.table(PredictTime= paste0(letters[i], "_", vecPercTimepoint[i]), Threshold=vecTROC[[i]]$Thresholds, 
                                 x=vecTROC[[i]]$FPR, y=vecTROC[[i]]$TPR))
  }
  vLabels[[i]] <- bquote("Prediction time "*italic(t)==.(vecPercTimepoint[i])*"; AUC: "*.(percent(vecTROC[[i]]$AUC, accuracy=0.01)))
}


# -- Graph a combined ROC-graph across prediction times t
# - Aesthetic parameters
datGraph[, FacetLabel := "Prentice-Williams-Peterson (PWP) spell-time model"]
vCol <- brewer.pal(8,"Set1")
vLabels_F <- setNames(vLabels, paste0(letters[1:length(vecPercTimepoint)],"_", vecPercTimepoint))
chosenFont <- "Cambria"

# - Create ROC-graph
(gg <- ggplot(datGraph, aes(x=x,y=y,group=PredictTime)) + theme_minimal() + 
    theme(text = element_text(family=chosenFont), legend.position="inside", 
          strip.background=element_rect(fill="snow2", colour="snow2"),
          strip.text=element_text(size=8, colour="gray50"), strip.text.y.right=element_text(angle=90),
          legend.position.inside = c(0.55,0.45),
          legend.background = element_rect(fill="snow2", color="black",
                                           linetype="solid", linewidth=0.1)) +
    labs(x = bquote("False Positive Rate "*italic(F^"+")), y = 
           bquote("True Positive Rate "*italic(T^"+"))) + 
    # Add 45-degree line
    geom_segment(x = 0, y = 0, xend = 1, yend = 1, color = "grey", linewidth=0.2) +
    # Main line graph
    geom_step(aes(x=x, y=y, linetype=PredictTime, colour=PredictTime), linewidth=0.05) + 
    geom_point(aes(x=x, y=y, shape=PredictTime, colour=PredictTime), size=0.25) +
    # Facets and scales
    facet_grid(FacetLabel ~ .) +  
    scale_color_manual(name=bquote("ROC"*(italic(t))), values=vCol, labels=vLabels) + 
    scale_linetype_discrete(name=bquote("ROC"*(italic(t))), labels=vLabels) + 
    scale_shape_discrete(name=bquote("ROC"*(italic(t))), labels=vLabels) + 
    scale_y_continuous(label=percent) + scale_x_continuous(label=percent))


# - Save graph
dpi <- 300
ggsave(gg, file=paste0(paste0(genFigPath, "PWP ST/DefaultSurvModel-Cox-PWPST-CombinedROC_Depedendence.png")), 
       width=1200/dpi, height=1000/dpi, dpi=dpi, bg="white")


# - cleanup
suppressWarnings( rm(gg, vLabels, vLabels_F, vecTROC, datGraph, dat, 
                     objROC1_TFD, objROC2_TFD, objROC3_TFD, objROC4_TFD, objROC1_AG, objROC2_AG, objROC3_AG, objROC4_AG,
                     objROC1_PWPST, objROC2_PWPST, objROC3_PWPST, objROC4_PWPST,
   km_TFD, km_TFD_survFitSummary, km_PWPST, km_PWPST_survFitSummary, cox_PWPST, cox_TFD, cox_AG,
   datCredit_train_TFD, datCredit_valid_TFD, datCredit_train_AG, datCredit_valid_AG,
   datCredit_train_PWPST, datCredit_valid_PWPST
   ) )
