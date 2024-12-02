# ========================= INVESTIGATING SURVIVALROC() ==========================
# As a poof of concept, we shall compare various functions from various packages
# in conducting time-dependent ROC-analyses on the same fitted Cox regression model,
# having used the prepared credit data
# --------------------------------------------------------------------------------
# PROJECT TITLE: Default Survival Modelling
# SCRIPT AUTHOR(S): Dr Arno Botha (AB)
# VERSION: 1.0 (November-2024)
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
#   - 3c.Data_Fusion2.R
#   - 5a(i).InputSpace_TFD.R

# -- Inputs:
#   - datCredit_train_TFD | Prepared from script 3b
#   - datCredit_valid_TFD | Prepared from script 3b

#
# -- Outputs:
#   - Input_Space
# ================================================================================





# ----------------- 1. Load data and conduct preliminary Kaplan-Meier (KM) analyses
# NOTE: Conducting KM-analyses helps identify interesting prediction times for the
# eventual ROC-analyses


# ------ Time to first Default (TFD) definition
# - Confirm prepared datasets are loaded into memory
if (!exists('datCredit_train_TFD')) unpack.ffdf(paste0(genPath,"creditdata_train_TFD"), tempPath);gc()
if (!exists('datCredit_valid_TFD')) unpack.ffdf(paste0(genPath,"creditdata_valid_TFD"), tempPath);gc()

# - Feature engineering
# NOTE: This will be moved to script 3c in due time
datCredit_train_TFD[, slc_acct_arr_dir_3_Change_Ind := ifelse(slc_acct_arr_dir_3 != "SAME", 1,0)]
datCredit_train_TFD[,TimeInDelinqState_Lag_1 := shift(TimeInDelinqState,fill=0),by=LoanID]
datCredit_valid_TFD[, slc_acct_arr_dir_3_Change_Ind := ifelse(slc_acct_arr_dir_3 != "SAME", 1,0)]
datCredit_valid_TFD[,TimeInDelinqState_Lag_1 := shift(TimeInDelinqState,fill=0),by=LoanID]

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



# ------ Prentice-Williams-Peterson (PWP) Total-time definition
### AB: Objects not yet created, but space reserved so long for the PWP-model
# - Confirm prepared datasets are loaded into memory
#if (!exists('datCredit_train_PWPST')) unpack.ffdf(paste0(genPath,"creditdata_train_PWPST"), tempPath);gc()
#if (!exists('datCredit_valid_PWPST')) unpack.ffdf(paste0(genPath,"creditdata_valid_PWPST"), tempPath);gc()





# ----------------- 2. Fit a Cox regression model on the resampled prepared data

# ------ Time to first Default (TFD) definition
# - Initialize variables | AB-variant
vecVars_TFD <- c("PerfSpell_g0_Delinq_Num", "Arrears", "g0_Delinq_Ave", "TimeInDelinqState_Lag_1",      
           "slc_acct_arr_dir_3_Change_Ind", "slc_acct_pre_lim_perc_imputed_med", 
           "InterestRate_Margin_Aggr_Med", "InterestRate_Margin_Aggr_Med_3", 
           "Principal", "LN_TPE", "M_DTI_Growth","M_Inflation_Growth",
           "M_Inflation_Growth_6", "M_RealIncome_Growth")

# - Fit a Cox Proportional Hazards model with time-varying covariates, and clustered observations
# NOTE: Assume dependence (by specifying ID-field) amongst certain observations clustered around ID-values
cox_TFD <- coxph(as.formula(paste0("Surv(Start,End,Default_Ind) ~ ",
                                   paste(vecVars_TFD,collapse=" + "))), id=LoanID, datCredit_train_TFD)
summary(cox_TFD)



# ------ Prentice-Williams-Peterson (PWP) Total-time definition
### AB: Objects not yet created, but space reserved so long for the PWP-model
# - Initialize variables | AB-variant
# vecVars_PWP <- c("PerfSpell_g0_Delinq_Num", "Arrears", "g0_Delinq_Ave", "TimeInDelinqState_Lag_1",      
#                  "slc_acct_arr_dir_3_Change_Ind", "slc_acct_pre_lim_perc_imputed_med", 
#                  "InterestRate_Margin_Aggr_Med", "InterestRate_Margin_Aggr_Med_3", 
#                  "Principal", "LN_TPE", "M_DTI_Growth","M_Inflation_Growth",
#                  "M_Inflation_Growth_6", "M_RealIncome_Growth")
# 
# # Fit a Cox Proportional Hazards model with time-varying covariates, and clustered observations
# # NOTE: Assume dependence (by specifying ID-field) amongst certain observations clustered around ID-values
# cox_PWP <- coxph(as.formula(paste0("Surv(Start,End,Default_Ind) ~ ",
#                                    paste(vecVars_PWP,collapse=" + "))), id=LoanID, datCredit_train_PWP)
# summary(cox_PWP)





# ----------------- 3. Conduct & compare ROC-analyses at fixed prediction time
# NOTE: ROC-analyses are conducted using the same fitted Cox model at a specific prediction time t

# - Decide prediction point t at which we would like to conduct a time-dependent ROC-analysis
# assuming the cumulative/dynamic (CD) context, i.e., that TPR is calculated given that 
# raw event times T_{ijt} are less than or equal to the given prediction time t
predictTime <- 12 # 90% percentile survival time (203); median: 338


# ------ Time to first Default (TFD) definition

# --- Package: survivalROC()
# Using survivalROC with NNE-estimator of S(t) under the CD-approach | NNE-kernel for S(t)

# - Calculate AUC up to given prediction time for correctly-fitted Cox model
# NOTE: Uses the superior Nearest Neighbour Estimator (NNE) method for S(t)
ptm <- proc.time() #IGNORE: for computation time calculation;
survivalROC::survivalROC(Stime=datCredit_train_TFD$End, status=datCredit_train_TFD$Default_Ind, entry=datCredit_train_TFD$Start, 
                         method = "NNE",  span=0.05, predict.time=predictTime,
                         marker=round(predict(cox_TFD, type="lp"),2))
proc.time() - ptm
### RESULTS: AUC: 94% % up to t, achieved in 8781.23 sec (146.35 mins)
### AB: When this is done running, note results, then restart PC for a fresh session, before testing tROC.multi()



# --- Package: tROCkit() | custom "package"/function
# NOTE: Using custom tROC()-function from script 0b(iii) under the CD-approach with an NN-estimator and 0/1-kernel

# -- Multi-threaded calculation of the # AUC from given start up to given prediction time 3 in following the CD-approach
# NOTE: Uses the superior Nearest Neighbour Estimator (NNE) method for S(t) with a 0/1-kernelNNE-kernel for S(t)
# NOTE2: Assume dependence (by specifying ID-field) amongst certain observations clustered around ID-values
ptm <- proc.time() #IGNORE: for computation time calculation;
predictTime <- 3
objROC1_TFD <- tROC.multi(datGiven=datCredit_valid_TFD, cox=cox_TFD, month_End=predictTime, sLambda=0.05, estMethod="NN-0/1", numDigits=2, 
                          fld_ID="PerfSpell_Key", fld_Event="PerfSpell_Event", eventVal=1, fld_StartTime="Start", fld_EndTime="End",
                          graphName="DefaultSurvModel-Cox1_Depedendence", genFigPath=paste0(genFigPath, "TFD/"), 
                          genObjPath=genObjPath, caseStudyName=paste0("TFD_", predictTime), numThreads=6)
objROC1_TFD$AUC; objROC1_TFD$ROC_graph
proc.time() - ptm
### RESULTS: AUC up to t: 91.21%, achieved in   secs ( mins)


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
objROC2_TFD <- tROC.multi(datGiven=datCredit_valid_TFD, cox=cox_TFD, month_End=predictTime, sLambda=0.05, estMethod="NN-0/1", numDigits=2,
                      fld_ID="PerfSpell_Key", fld_Event="PerfSpell_Event", eventVal=1, fld_StartTime="Start", fld_EndTime="End",
                      graphName="DefaultSurvModel-Cox1_Depedendence", genFigPath=paste0(genFigPath, "TFD/"),
                      genObjPath=genObjPath, caseStudyName=paste0("TFD_", predictTime), numThreads=6)
objROC2_TFD$AUC; objROC2_TFD$ROC_graph
proc.time() - ptm
### RESULTS: AUC up to t: %, achieved in   secs ( mins)


# -- Multi-threaded calculation of the # AUC from given start up to given prediction time 36 in following the CD-approach
# NOTE: Uses the superior Nearest Neighbour Estimator (NNE) method for S(t) with a 0/1-kernelNNE-kernel for S(t)
# NOTE2: Assume dependence (by specifying ID-field) amongst certain observations clustered around ID-values
ptm <- proc.time() #IGNORE: for computation time calculation;
predictTime <- 24
objROC3_TFD <- tROC.multi(datGiven=datCredit_valid_TFD, cox=cox_TFD, month_End=predictTime, sLambda=0.05, estMethod="NN-0/1", numDigits=2, 
                          fld_ID="PerfSpell_Key", fld_Event="PerfSpell_Event", eventVal=1, fld_StartTime="Start", fld_EndTime="End",
                          graphName="DefaultSurvModel-Cox1_Depedendence", genFigPath=paste0(genFigPath, "TFD/"), 
                          genObjPath=genObjPath, caseStudyName=paste0("TFD_", predictTime), numThreads=6)
objROC3_TFD$AUC; objROC3_TFD$ROC_graph
proc.time() - ptm
### RESULTS: AUC up to t: %, achieved in   secs ( mins)


# -- Multi-threaded calculation of the # AUC from given start up to given prediction time 36 in following the CD-approach
# NOTE: Uses the superior Nearest Neighbour Estimator (NNE) method for S(t) with a 0/1-kernelNNE-kernel for S(t)
# NOTE2: Assume dependence (by specifying ID-field) amongst certain observations clustered around ID-values
ptm <- proc.time() #IGNORE: for computation time calculation;
predictTime <- 36
objROC4_TFD <- tROC.multi(datGiven=datCredit_valid_TFD, cox=cox_TFD, month_End=predictTime, sLambda=0.05, estMethod="NN-0/1", numDigits=2, 
                          fld_ID="PerfSpell_Key", fld_Event="PerfSpell_Event", eventVal=1, fld_StartTime="Start", fld_EndTime="End",
                          graphName="DefaultSurvModel-Cox1_Depedendence", genFigPath=paste0(genFigPath, "TFD/"), 
                          genObjPath=genObjPath, caseStudyName=paste0("TFD_", predictTime), numThreads=6)
objROC4_TFD$AUC; objROC4_TFD$ROC_graph
proc.time() - ptm
### RESULTS: AUC up to t: %, achieved in   secs ( mins)







# ------ Prentice-Williams-Peterson (PWP) Total-time definition
### AB: Objects not yet created, but space reserved so long for the PWP-model

# --- Package: tROCkit() | custom "package"/function
# Using custom tROC()-function from script 0b(iii) under the CD-approach with an NN-estimator
# - Calculate AUC from given start up to given prediction time in following the CD-approach
# NOTE: Uses the superior Nearest Neighbour Estimator (NNE) method for S(t) with a 0/1-kernelNNE-kernel for S(t)
# NOTE2: Assume dependence (by specifying ID-field) amongst certain observations clustered around ID-values
# ptm <- proc.time() #IGNORE: for computation time calculation;
# objROC1_PWP <- tROC.multi(datGiven=datCredit_valid_PWP, cox=cox_PWP, month_End=12, sLambda=0.05, estMethod="NN-0/1", numDigits=2, 
#                           fld_ID="PerfSpell_Key", fld_Event="PerfSpell_Event", eventVal=1, fld_StartTime="Start", fld_EndTime="End",
#                           graphName="DefaultSurvModel-Cox1_Depedendence", genFigPath=paste0(genFigPath, "PWP/"), numThreads=5)
# objROC1_PWP$AUC; objROC1_PWP$ROC_graph
# proc.time() - ptm
### RESULTS: AUC up to t: %, achieved in  secs ( mins)
