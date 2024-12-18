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
#   - 5a(iv).InputSpace_PWPST.R

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
if (!exists('datCredit_train_PWPST')) unpack.ffdf(paste0(genPath,"creditdata_train_PWPST"), tempPath);gc()
if (!exists('datCredit_valid_PWPST')) unpack.ffdf(paste0(genPath,"creditdata_valid_PWPST"), tempPath);gc()

# - Feature engineering
# NOTE: This will be moved to script 3c in due time
datCredit_train_TFD[, slc_acct_arr_dir_3_Change_Ind := ifelse(slc_acct_arr_dir_3 != "SAME", 1,0)]
datCredit_valid_TFD[, slc_acct_arr_dir_3_Change_Ind := ifelse(slc_acct_arr_dir_3 != "SAME", 1,0)]
datCredit_train_PWPST[, slc_acct_arr_dir_3_Change_Ind := ifelse(slc_acct_arr_dir_3 != "SAME", 1,0)]
datCredit_valid_PWPST[, slc_acct_arr_dir_3_Change_Ind := ifelse(slc_acct_arr_dir_3 != "SAME", 1,0)]





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
# - Initialize variables | AB-variant
vecVars_PWPST <- c("g0_Delinq_SD_4", "Arrears", "g0_Delinq_Ave", "TimeInDelinqState_Lag_1",      
                 "slc_acct_arr_dir_3_Change_Ind", "slc_acct_pre_lim_perc_imputed_med", 
                 "InterestRate_Nom", "M_Repo_Rate", "M_Repo_Rate_9", 
                 "InstalmentToBalance_Aggr_Prop", "LN_TPE", "AgeToTerm",
                 "NewLoans_Aggr_Prop", "M_DTI_Growth","M_Inflation_Growth",
                 "M_Inflation_Growth_6", "M_RealIncome_Growth", "PerfSpell_Grp")

# # Fit a Cox Proportional Hazards model with time-varying covariates, and clustered observations
# # NOTE: Assume dependence (by specifying ID-field) amongst certain observations clustered around ID-values
cox_PWP <- coxph(as.formula(paste0("Surv(Start,End,Default_Ind) ~ ",
                                   paste(vecVars_PWPST,collapse=" + "))), id=LoanID, datCredit_train_PWPST)
summary(cox_PWP)





# ----------------- 3. Actual & Expected Term structure for a chosen cohort

# ------ Time to first Default (TFD) definition

# --- Preliminaries
# - Define cohort and filter data
sCohort <- "2007-11-30"
vSpellKeys <- datCredit_valid_TFD[Date == as.Date(sCohort) & Start == 0, PerfSpell_Key]
datCohort <- datCredit_valid_TFD[PerfSpell_Key %in% vSpellKeys]
check1 <- subset(datCohort, PerfSpell_Key %in% unique(datCohort[,PerfSpell_Key])[1])


# --- Fit Kaplan-Meier (KM) nonparametric model towards calculating "Actual Hazard"
# Compute Kaplan-Meier survival estimates (product-limit) for main-event | Spell-level with right-censoring & left-truncation
# All competing events preclude the main event from happening and are therefore considered as censored
# ID is set as the spell key, with no stratification
km_TFD <- survfit(Surv(time=Start, time2=End, event=Default_Ind==1,type="counting") ~ 1, 
                  id=LoanID, data=datCredit_train_TFD)
summary(km_TFD)$table # overall summary statistics
### RESULTS: 62k events, with no median survival probability
(km_TFD_survFitSummary <- surv_summary(km_TFD)) # Survival table

# - Create hazard table from KM-object towards creating "Actual Hazard"
datHaz <- data.table(Time = km_TFD$time, AtRisk_n = km_TFD$n.risk, Event_n = km_TFD$n.event, Censored_n=km_TFD$n.censor,
                     Actual_Hazard = km_TFD$n.event/km_TFD$n.risk, Surv_KM = km_TFD$surv) %>% 
  mutate(CHaz = cumsum(Actual_Hazard))  %>% 
  filter(Event_n > 0 | Censored_n >0)
describe(datHaz$Actual_Hazard); hist(datHaz$Actual_Hazard, breaks="FD")
# fuse with main dataset
datCohort <- merge(datCohort, datHaz,by.x="End", by.y="Time", all.x=T)

# - Obtain both cumulative and incremental baseline hazard functions H_0(t) and h_0(t) from Cox PH model
datCHaz <- basehaz(cox_TFD, centered=F) %>% as.data.table() # ensure Cumulative Hazard corresponds to unadjusted (baseline) subjects
# Approximate H_0'(t) [incremental hazard] over small periods by taking finite forward differences
datCHaz[, BaselineHaz :=  c(diff(hazard), NA)]

# - Match time points to baseline hazard at a specific times
datCohort <- merge(datCohort, datCHaz[,list(End=time, CHaz_CoxPH = hazard, BaselineHaz)],
                   by="End", all.x=T)

# - Score data using Cox PH model and hazard as conditional probability of defaulting in 1 month intervals
datCohort[, Risk_Score := predict(cox_TFD, newdata = datCohort, type = "risk", id=PerfSpell_Key)]
datCohort[, Predicted_SurvProb := exp(-CHaz), by=list(PerfSpell_Key)]
# Obtain incremental hazards by finite forward differencing
datCohort[, Predicted_Hazard := c(diff(CHaz),NA), by=list(PerfSpell_Key)]
datCohort[, Predicted_Hazard2 := BaselineHaz * Risk_Score] # Wrong definition somehow
describe(datCohort$Predicted_SurvProb); hist(datCohort$Predicted_SurvProb)
describe(datCohort$Predicted_Hazard); hist(datCohort$Predicted_Hazard)
describe(datCohort$Predicted_Hazard2); hist(datCohort$Predicted_Hazard2)


# - Fit spline regression to actual hazards
sKnots <- 10
datHaz[, Actual_Hazard_Spline := spline_estimation(Time, Actual_Hazard, sKnots, 2)] # Actual hazard spline








