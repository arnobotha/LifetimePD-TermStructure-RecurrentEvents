# ========================= INVESTIGATING SURVIVALROC() ==========================
# As a poof of concept, we shall compare various functions from various packages
# in conducting time-dependent ROC-analyses on the same fitted Cox regression model,
# having used the prepared credit data
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





# ----------------- 1. Load data

# - General parameters
sMaxSpellAge <- 300 # max for [PerfSpell_Age], as determined in earlier analyses (script 4a(i))

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
cox_TFD <- coxph(as.formula(paste0("Surv(Start,End,Default_Ind) ~ ", paste(vecVars_TFD,collapse=" + "))), 
                 ties="efron", id=LoanID, datCredit_train_TFD)
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
cox_PWP <- coxph(as.formula(paste0("Surv(Start,End,Default_Ind) ~ ",paste(vecVars_PWPST,collapse=" + "))), 
                 ties="efron", id=PerfSpell_Key, datCredit_train_PWPST)
summary(cox_PWP)





# ----------------- 3. Actual & Expected Term structure for a chosen cohort

# ------ Time to first Default (TFD) definition


# --- Preliminaries
# - Define cohort and filter data
sCohort <- "2007-11-30"
vSpellKeys <- datCredit_valid_TFD[Date == as.Date(sCohort) & Start == 0, PerfSpell_Key]
datCohort <- datCredit_valid_TFD[PerfSpell_Key %in% vSpellKeys]
check1 <- subset(datCohort, PerfSpell_Key %in% unique(datCohort[,PerfSpell_Key])[1],
                 select=c("LoanID","PerfSpell_Key", "PerfSpell_Num", "TimeInPerfSpell","PerfSpell_Age", "PerfSpellResol_Type_Hist"))
describe(datCohort[PerfSpell_Counter==1, PerfSpell_Age]); hist(datCohort[PerfSpell_Counter==1, PerfSpell_Age], breaks="FD")



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
                     Hazard_Actual = km_TFD$n.event/km_TFD$n.risk, Surv_KM = km_TFD$surv) %>% 
  mutate(CHaz = cumsum(Hazard_Actual),
         CHaz2 = -log(km_TFD$surv))  %>% # Very similar to CHaz, derived fundamentally from KM-estimate, kept for comparative purposes
  filter(Event_n > 0 | Censored_n >0)
describe(datHaz$Hazard_Actual); hist(datHaz$Hazard_Actual, breaks="FD")
plot(datHaz[Time <= sMaxSpellAge, Hazard_Actual], type="b") # restrict time just for quick plotting purpose
### RESULTS: Functional shape exhibits typical U-shape, as expected.


# --- Calculate fundamental quantities and enrich main data object accordingly


# - Fuse hazards with main dataset
#datCohort <- merge(datCohort, datHaz[, list(Time, Hazard_Actual, CHaz)], by.x="End", by.y="Time", all.x=T)

# - Obtain both cumulative and incremental baseline hazard functions H_0(t) and h_0(t) from Cox PH model
# NOTE: Uses Efron's method to obtain baseline hazard whereby risk sets are adjusted based on number of tied events
datSurv_TFD <- survfit(cox_TFD, newdata=datCredit_valid_TFD, id=LoanID)

datCHaz <- basehaz(cox_TFD, centered=F) %>% as.data.table() # ensure Cumulative Hazard corresponds to unadjusted (baseline) subjects

# - Obtain survival curve for new data


# Initialize interpolated baseline hazard
vBaselineHaz <- numeric(length(datCHaz$time))

# Loop through each time point to perform interpolation
for (i in seq_along(datCHaz$time)) {
  # i <- 889
  # Get list of indices where conditions hold true
  lower_index <- which(datCHaz$time <= datCHaz$time[i])
  upper_index <- which(datCHaz$time > datCHaz$time[i])
  
  # Apply min and max respectively on lower and upper indices towards getting bounds
  if (length(lower_index) > 1) lower_index <- max(lower_index) else lower_index <- NA
  if (length(upper_index) > 1)  upper_index <- min(upper_index) else upper_index <- NA
  
  if (is.na(lower_index)) {
    vBaselineHaz[i] <- datCHaz$hazard[1]
  } else if (is.na(upper_index)) {
    vBaselineHaz[i] <- datCHaz$hazard[length(datCHaz$hazard)]
  } else {
    t1 <- datCHaz$time[lower_index]
    t2 <- datCHaz$time[upper_index]
    h1 <- datCHaz$hazard[lower_index]
    h2 <- datCHaz$hazard[upper_index]
    
    # Interpolate linearly
    vBaselineHaz[i] <- h1 + ( (h2 - h1) * (datCHaz$time[i] - t1) / (t2 - t1) )
  }
}

datCHaz[, BaselineHaz := vBaselineHaz]








### AB [2024-12-18]: The code below is still very messy. This is just proof/disproof of concepts in trying to untangle this problem


# --- Score data using Cox PH model 
# We want to focus on the discrete-time hazard as the conditional probability of defaulting at t in 1 month intervals,
# having survived up to time t
datCohort[, RiskScore := predict(cox_TFD, newdata = datCohort, type = "risk", id=PerfSpell_Key)]

# - Obtain unique event times for those cases that experienced the event
#vEventTimes <- sort(unique(datCohort$End[datCohort$Default_Ind==1]))
vEventTimes <- sort(unique(datCohort$End))

# - Initialize data object for baseline hazard calculation
datEventTimes_Haz <- data.table(time = vEventTimes, n_event = numeric(length(vEventTimes)),
                                total_risk = numeric(length(vEventTimes))
)

# - Calculate number of events and total risk scores for each ordered failure time
for (i in seq_along(vEventTimes)) {
  # current event time t
  t <- vEventTimes[i]
  # Number of events
  datEventTimes_Haz$n_event[i] <- sum(datCohort$End == t & datCohort$Default_Ind == 1)
  # Get indices of those at risk at current event time t
  vAtRisk <- which(datCohort$Start <= t & datCohort$End > t)
  
  # Get Total risk scores for subjects at risk at t
  datEventTimes_Haz$total_risk[i] <- sum(sim_data$RiskScore[vAtRisk])
}


# - Match time points to baseline hazard at a specific times
datCohort <- merge(datCohort, datCHaz[,list(End=time, CHaz_CoxPH = hazard, BaselineHaz)],
                   by="End", all.x=T)



datCohort[, ]
# Obtain survival probability  as a fundamental quantity
datCohort[, Predicted_SurvProb := exp(-CHaz), by=list(PerfSpell_Key)]

# - Obtain incremental baseline hazards h_0 from cumulative baseline hazard H_0
# Note: h_0(t) is the derivative of H_0(t) over t. In obtaining this numerically,
# we can use finite (forward) differencing as an approximation as H_0(t_{i+1}) - H_0(t_i) / (t_{i+1} - t_i)
datCohort[, Hazard_Predicted := c(diff(CHaz_CoxPH),NA), by=list(PerfSpell_Key)]
#datCohort[, Hazard_Predicted2 := BaselineHaz * Risk_Score] # Wrong definition somehow
# Conduct distributional analyses
describe(datCohort$Predicted_SurvProb); hist(datCohort$Predicted_SurvProb)
describe(datCohort$Hazard_Predicted); hist(datCohort$Hazard_Predicted)
describe(datCohort$Hazard_Actual); hist(datCohort$Hazard_Actual)
plot(datCohort[PerfSpell_Age<=300, Hazard_Actual], type="b")
#describe(datCohort$Hazard_Predicted2); hist(datCohort$Hazard_Predicted2)
### RESULTS: Different distributions of hazard functions

# --- Data filtration and aggregation

# Filter for outlying spell ages, as determined in earlier analyses ()
datCohort_sub <- datCohort %>% subset(PerfSpell_Age <= sMaxSpellAge)
percent(datCohort_sub[, .N] / datCohort[,.N], accuracy=0.01)
### RESULTS: Still 100% of data

# - Get unique event times, 
vEventTimes <- unique(datCohort_sub[, End])


# --- Fitting regression splines as way to summarise point-estimates at each event time t
# - Fit cubic regression splines to actual hazard for graphing purposes
sKnots_Act <- 10
splHazard_Act <- lm(Hazard_Actual ~ bs(Hazard_Actual, knots=sKnots_Act, degree=3), data=datCohort)
summary(splHazard_Act)
datCohort[!is.na(Hazard_Actual),, Hazard_Actual_Spline := predict(splHazard_Act, newx=Hazard_Actual)]

# - Fit cubic regression splines to expected hazard for graphing purposes
sKnots_Exp <- 10
splHazard_Exp <- lm(Hazard_Predicted ~ bs(Hazard_Predicted, knots=sKnots_Exp, degree=3), data=datCohort)
summary(splHazard_Exp)
datCohort[!is.na(Hazard_Predicted), Hazard_Predicted_Spline := predict(splHazard_Exp, newx=Hazard_Predicted)]


# --- Data aggregation and graphing preparation
# - Select and pivot dataset to long format for graphing purposes
datAggr <- subset(datCohort, select=c ("PerfSpell_Key", "End", "PerfSpellResol_Type_Hist", "PerfSpell_Age",
                                "CHaz", "CHaz_CoxPH", "Hazard_Actual", "BaselineHaz", "Hazard_Predicted",
                                "Hazard_Actual_Spline", "Hazard_Predicted_Spline")) %>% 
  setorder(PerfSpell_Key, End) %>%
  pivot_longer(cols=!(c("PerfSpell_Key", "End", "PerfSpellResol_Type_Hist", "PerfSpell_Age")), 
               names_to = "Type", values_to = "Value") %>% as.data.table()
check2 <- subset(datAggr, PerfSpell_Key %in% unique(datCohort[,PerfSpell_Key])[1])

# - Create different groupings for graphing purposes
datAggr[, ColourVar1 := paste0(Type, "_Line")]
datAggr[, ColourVar2 := paste0(Type, "_Point")]


# --- Creating the main graph

# - Graphing parameters
vSelected <- c("Hazard_Actual", "Hazard_Actual_Spline", "Hazard_Predicted", "Hazard_Predicted_Spline")
vCol_Point <- brewer.pal(8, "Pastel1")[c(1,2)]
vCol_Line <- brewer.pal(8, "Set1")[c(1,2)]
mainEventName <- "Default"
chosenFont <- "Cambria"
vLabel <- c("Actual Probability", "Actual Probability Cubic Spline")
vLabel2 <- c("A_Actual_Line"="Spline: Actual", "B_Predicted_Line"="Spline: Scored",
             "A_Actual_Point"="Case: Actual", "B_Predicted_Point"="Case: Scored")
vValues <- c("A_Actual_Line"=vCol_Line[1], "B_Predicted_Line"=vCol_Line[2],
             "A_Actual_Point"=vCol_Point[1], "B_Predicted_Point"=vCol_Point[2])

ggplot(subset(datAggr, Type %in% vSelected & PerfSpell_Age <= 300), aes(y=Value, x=End, group=Type)) + theme_minimal() +
  geom_point(aes(y=Value, color=ColourVar2))



