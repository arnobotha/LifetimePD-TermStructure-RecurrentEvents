# ============================================ INPUT SPACE: TFD ========================================
# Divide data into thematic groups and perform data analysis on them to compile an input space for the TFD model
# ------------------------------------------------------------------------------------------------------
# PROJECT TITLE: Default survival modelling
# SCRIPT AUTHOR(S): Bernard Scheepers (BS)
# ------------------------------------------------------------------------------------------------------
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

# -- Inputs:
#   - datCredit_train_TFD | Prepared from script 3c
#   - datCredit_valid_TFD | Prepared from script 3c

#
# -- Outputs:
#   - Input_Space
# ------------------------------------------------------------------------------------------------------

# ------ 0. Preliminaries
# - Confirm prepared datasets are loaded into memory
if (!exists('datCredit_train_TFD')) unpack.ffdf(paste0(genPath,"creditdata_train_TFD"), tempPath);gc()
if (!exists('datCredit_valid_TFD')) unpack.ffdf(paste0(genPath,"creditdata_valid_TFD"), tempPath);gc()





#============================================================================================
# ------ 1. Delinquency measures
# - Vector storing the remaining thematic variable names and type
varlist <- data.table(vars=c("g0_Delinq","g0_Delinq_fac","PerfSpell_g0_Delinq_Num",
                             "Arrears" ,"TimeInDelinqState","g0_Delinq_Any_Aggr_Prop",
                             "g0_Delinq_Ave","slc_acct_arr_dir_3",
                             "slc_acct_roll_ever_24_imputed_mean"),
                      vartypes=c("acc", "cat", "acc", "acc", "acc", "dte", "dte",
                                 "cat","num"))

#=========================================================================================



# ------ 1.1 Which time window length is the best in calculating Delinquency volatility?

# Initialize variables to be tested
vars <- c("g0_Delinq_SD_4", "g0_Delinq_SD_5", "g0_Delinq_SD_6", "g0_Delinq_SD_9", "g0_Delinq_SD_12")

# Goodness of fit test
csTable(datCredit_train_TFD, vars, TimeDef="TFD") # [g0_Delinq_SD_4] has the lowest KS-statistic

aicTable(datCredit_train_TFD, vars, TimeDef="TFD", genPath=genObjPath) # [g0_Delinq_SD_5] has the lowest AIC

# Accuracy test
concTable(datCredit_train_TFD, datCredit_valid_TFD, vars, TimeDef="TFD", genPath=genObjPath) # [g0_Delinq_SD_4] has the highest concordance

### RESULTS: As the SD period increase, there is a slight decrease in concordance.

### Conclusion: Larger window are less influenced by large changes therefore significant changes
###             are less pronounced. Include [g0_Delinq_SD_4] in the varlist.

varlist <- vecChange(varlist,Remove="PerfSpell_g0_Delinq_SD",Add=data.table(vars=c("g0_Delinq_SD_4"), vartypes=c("acc")))



# ------ 1.2 Which variables are highly correlated?

# Correlation analysis
corrAnalysis(datCredit_train_TFD, varlist[vartypes!="cat"]$vars, corrThresh = 0.6, method = 'spearman') # Obtain correlation groups

### RESULTS:  1) [g0_Delinq_Any_Aggr_Prop] and [g0_Delinq_Ave] with a correlation of 1
###           2) [g0_Delinq] and [Arrears]

### CONCLUSION: A single variable from each group must be retained while the rest are removed.


# ------ 1.2.1 Which variable should be kept from group 1) [g0_Delinq_Any_Aggr_Prop] and [g0_Delinq_Ave]

vars <- c("g0_Delinq_Any_Aggr_Prop", "g0_Delinq_Ave")

# Goodness of fit
csTable(datCredit_train_TFD, vars, TimeDef="TFD") # No noticeable difference

AICTable(datCredit_train_TFD, vars, TimeDef="TFD", genPath=genObjPath) # [g0_Delinq_Ave] has the lowest AIC

### RESULTS: [g0_Delinq_Ave] seems to have the better goodness of fit.

# Accuracy
concTable(datCredit_train_TFD, datCredit_valid_TFD, vars, TimeDef="TFD", genPath=NA) # [g0_Delinq_Ave] has the highest concordance

### CONCLUSION: [g0-Delinq_Ave] seems to outperform [g0_Delinq_Any_Aggr_Prop] and therefore is kept in the model
###             and [g0_Delinq_Any_Aggr_Prop_Lag_3] is removed along with [g0_Delinq_Any_Aggr_Prop] due to the high correlation

varlist <- vecChange(varlist,Remove=c("g0_Delinq_Any_Aggr_Prop"))



# ------ 1.3 Which version of [g0_Delinq] should be kept in the model?

# [g0_Delinq]
cox <- coxph(Surv(Start,End,Default_Ind) ~ g0_Delinq, id=LoanID, datCredit_train_TFD)
summary(cox);rm(cox)
### RESULTS: Beta is unstable with high variability.

# [g0_Delinq_fac]
cox <- coxph(Surv(Start,End,Default_Ind) ~ g0_Delinq_fac, id=LoanID, datCredit_train_TFD)
summary(cox);rm(cox)
### RESULTS: coef is unstable with high variability.

### INVESTIGATE: Is quasi-complete seperation present?
datCredit_train_TFD[g0_Delinq==3 & Default_Ind==0, .N]
### RESULTS: 0
datCredit_train_TFD[g0_Delinq!=3 & Default_Ind==1, .N]
### RESULTS: 0

### CONCLUSION: Quasi-complete separation is present for [g0_Delinq]=3

varlist <- vecChange(varlist,Remove=c("g0_Delinq","g0_Delinq_fac"))

# Arrears
cox <- coxph(Surv(Start,End,Default_Ind) ~ Arrears, id=LoanID, datCredit_train_TFD)
summary(cox); rm(cox)
### RESULTS: Obtained a stable model

### CONCLUSION: Unable to add [Delinq_0] to the model, since the various forms'
###             coef are unstable, however, [Arrears] fits a seemingly stable
###             model, therefore it can serve as a proxy for [g0_Delinq.]



# ------ 1.4 How does [PerfSpell_g0_Delinq_Num] compare to [g0_Delinq_SD_4]?

# Initialize variables
vars <- c("PerfSpell_g0_Delinq_Num", "g0_Delinq_SD_4")

# Goodness of fit
csTable(datCredit_train_TFD, vars, TimeDef="TFD") # [PerfSpell_g0_Delinq_Num] has the lowest KS-statistic

aicTable(datCredit_train_TFD, vars, TimeDef="TFD", genPath=genObjPath) # [g0_Delinq_SD_4] has the lowest AIC

# Accuracy
concTable(datCredit_train_TFD, datCredit_valid_TFD, vars, TimeDef="TFD", genPath=genObjPath) # [g0_Delinq_SD_4] has the highest concordance

### Conclusion: Remove [PerfSpell_g0_Delinq_Num]

varlist <- vecChange(varlist,Remove="PerfSpell_g0_Delinq_Num")



# ------ 1.5 What is the performance of current thematic variables in univariate models?

# Build thematic model based on remaining delinquency variables.
vars <- c("PerfSpell_g0_Delinq_Num","TimeInDelinqState","g0_Delinq_Ave",
          "slc_acct_arr_dir_3", "Arrears")

# Goodness of fit test
csTable(datCredit_train_TFD, vars, TimeDef="TFD")

### NOTE: [TimeInDelinqState] ran out of iterations and did not converge
### NOTE: [slc_acct_arr_dir_3] exp overflow due to covariates



# ------ 1.5.1 Why does [TimeInDelinqState] not converge?

cox <- coxph(Surv(Start,End,Default_Ind) ~ TimeInDelinqState, id=LoanID, datCredit_train_TFD)
### RESULTS:  Beta tends to Inf. After some inspection on the data it relates to all
###           Defaulting events starting in a new delinquency state,
###           i.e. [TimeInDelinqState] = 1. Therefore quasi-complete separation
###           seems to be present.

### INVESTIGATE: WHETHER QUASI-COMPLETE SEPERATION IS PRESENT
datCredit_train_TFD[TimeInDelinqState!=1 & Default_Ind==1, .N]
### RESULTS: 0
datCredit_train_TFD[TimeInDelinqState==1 & Default_Ind==0, .N]
### RESULTS: 169272

### CONCLUSION: Quasi-complete separation seems to be present for [g0_Delinq]=3
###             and should therefore be removed.

# Test a lagged version of TimeInDelinqState
cox <- coxph(Surv(Start,End,Default_Ind) ~ TimeInDelinqState_Lag_1, id=LoanID,
             datCredit_train_TFD)
summary(cox); rm(cox)

### CONCLUSION: Replace old variable with new variable, since it fits a stable model.

varlist <- vecChange(varlist,Remove="TimeInDelinqState",
                     Add=data.table(vars="TimeInDelinqState_Lag_1",vartypes="acc"))



# ------ 1.5.2 Why does [slc_acct_arr_dir_3] exp overflow?

describe(datCredit_train_TFD[,slc_acct_arr_dir_3])
# Value            CURING MISSING_DATA      ROLLING         SAME
# Proportion        0.014        0.125        0.022        0.838

describe(datCredit_train_TFD[Default_Ind==1,slc_acct_arr_dir_3])
# Value            CURING MISSING_DATA      ROLLING         SAME
# Proportion        0.004        0.157        0.788        0.051

### RESULTS:  Although the ROLLING level is not dominant in datCredit_train_TFD,
###           we can clearly see that it is highly predictive of a default event
###           occurring, given its high proportion if Default_Ind == 1.

# Use an indicator variable for a change in account
cox <- coxph(Surv(Start,End,Default_Ind) ~ slc_acct_arr_dir_3_Change_Ind,
             datCredit_valid_TFD, id=LoanID)
summary(cox); rm(cox)

### CONCLUSION: Replace old variable with new variable, since resulting model is stable.

varlist <- vecChange(varlist,Remove=c("slc_acct_arr_dir_3") ,
                     Add=data.table(vars=c("slc_acct_arr_dir_3_Change_Ind"),
                                    vartypes=c("bin")))



# ------ 1.5.3 What is the performance of current thematic variables in univariate models?

vars <- c("g0_Delinq_SD_4","g0_Delinq_Ave", "Arrears",
          "TimeInDelinqState_Lag_1","slc_acct_arr_dir_3_Change_Ind",
          "slc_acct_roll_ever_24_imputed_mean")

# Goodness of fit test
csTable(datCredit_train_TFD, vars, TimeDef="TFD")
#                         Variable      D
# 2                      g0_Delinq_Ave 0.6540
# 3                            Arrears 0.6539
# 1            PerfSpell_g0_Delinq_Num 0.6520
# 5      slc_acct_arr_dir_3_Change_Ind 0.6514
# 4            TimeInDelinqState_Lag_1 0.6445
# 6 slc_acct_roll_ever_24_imputed_mean 0.6336

#HERE
aicTable(datCredit_train_TFD, vars, TimeDef="TFD", genPath=genObjPath)
#                         Variable      AIC
# 1: slc_acct_roll_ever_24_imputed_mean 144336.7
# 2:            TimeInDelinqState_Lag_1 156436.3
# 3:      slc_acct_arr_dir_3_Change_Ind 160294.1
# 4:            PerfSpell_g0_Delinq_Num 182572.5
# 5:                            Arrears 183412.2
# 6:                      g0_Delinq_Ave 187989.6

# Accuracy test
concTable(datCredit_train_TFD, datCredit_valid_TFD, vars, TimeDef="TFD", genPath=genObjPath)
#                             Variable Concordance           SD LR_Statistic
# 1:                            Arrears   0.9606159 0.0016931039         6450
# 2:            PerfSpell_g0_Delinq_Num   0.9468068 0.0006694657         7290
# 3: slc_acct_roll_ever_24_imputed_mean   0.9205182 0.0017207379        45526
# 4:            TimeInDelinqState_Lag_1   0.9020457 0.0029892384        33426
# 5:      slc_acct_arr_dir_3_Change_Ind   0.8770983 0.0013364884        29569
# 6:                      g0_Delinq_Ave   0.5837088 0.0042886285         1873

### CONCLUSION: Leave all variables in the model.



# ------ 1.5 What is the performance of current thematic cox ph model?

# Build cox model based on all thematic variables
coxDelinq <- coxph(Surv(Start,End,Default_Ind) ~ g0_Delinq_Ave +
                           g0_Delinq_SD_4 + slc_acct_arr_dir_3_Change_Ind +
                           TimeInDelinqState_Lag_1 + Arrears + slc_acct_roll_ever_24_imputed_mean
                   , id=LoanID, data=datCredit_train_TFD)

summary(coxDelinq)

### RESULTS: All variables are significant (p-value < 0.001)

# Test goodness of fit
GoF_CoxSnell_KS(coxDelinq,datCredit_train_TFD,GraphInd = FALSE) # 0.6202

AIC(coxDelinq) # 94583.67

# Test accuracy
concordance(coxDelinq, newdata=datCredit_valid_TFD) # Concordance= 0.9963 se= 0.0002313

### CONCLUSION: Keep all variables in the model

# House keeping
rm(coxDelinq)

#============================================================================================

modelVar <- c("g0_Delinq_SD_4", "Arrears", "g0_Delinq_Ave",
              "TimeInDelinqState_Lag_1", "slc_acct_arr_dir_3_Change_Ind",
              "slc_acct_roll_ever_24_imputed_mean")

#============================================================================================





#============================================================================================
# ------ 2. Engineered measures
varlist <- data.table(vars=c("slc_acct_pre_lim_perc_imputed_med",
                             "slc_acct_prepaid_perc_dir_12_imputed_med",
                             "pmnt_method_grp") ,
                      vartypes=c("prc", "dec", "cat"))
#=========================================================================================



# ------ 2.1 Which variables should be removed due to high correlation?

# - Correlation analysis
corGroups <- corrAnalysis(datCredit_train_TFD, varlist[vartypes!="cat"]$vars, corrThresh = 0.6, method = 'spearman') # Obtain correlation groups

### RESULTS:  1) [slc_acct_pre_lim_perc_imputed_med] and [slc_acct_prepaid_perc_dir_12_imputed_med] with a correlation of 0.82

# Initialize variables to be tested
vars <- c("slc_acct_pre_lim_perc_imputed_med", "slc_acct_prepaid_perc_dir_12_imputed_med")

# Goodness of fit test
csTable(datCredit_train_TFD, vars, TimeDef="TFD") # No noticeable difference found

aicTable(datCredit_train_TFD, vars, TimeDef="TFD", genPath=genObjPath) # [slc_acct_pre_lim_perc_imputed_med] has the lowest AIC

# Accuracy test
concTable(datCredit_train_TFD, datCredit_valid_TFD, vars, TimeDef="TFD", genPath=genObjPath) # [slc_acct_pre_lim_perc_imputed_med] has the highest concordance

varlist <- vecChange(varlist,Remove="slc_acct_prepaid_perc_dir_12_imputed_med")



# ------ 2.2 What is the predictive power of the variables left in varlist?

# Initialize variables to be tested
vars <- c("pmnt_method_grp", "slc_acct_pre_lim_perc_imputed_med")

csTable(datCredit_train_TFD, vars, TimeDef="TFD") # [slc_acct_pre_lim_perc_imputed_med ] has the lowest KS-statistic

aicTable(datCredit_train_TFD, vars, TimeDef="TFD", genPath=genObjPath) # [pmnt_method_grp] has the lowest AIC

# Accuracy test
concTable(datCredit_train_TFD, datCredit_valid_TFD, vars, TimeDef="TFD", genPath=genObjPath) # [pmnt_method_grp] has the highest concordance

### Conclusion: All values should be kept in the model



# ------ 2.3 How predictive is a single model based on the thematic variables?

# Initialize variables
vars <- c("pmnt_method_grp", "slc_acct_pre_lim_perc_imputed_med")

# Build Cox model based on each variable
coxEngineered <- coxph(Surv(Start, End, Default_Ind) ~ pmnt_method_grp +
                               slc_acct_pre_lim_perc_imputed_med, id=LoanID,
                             data=datCredit_train_TFD)
summary(coxEngineered)


# Goodness of fit
GoF_CoxSnell_KS(coxEngineered,datCredit_train_TFD, GraphInd=FALSE) # 0.6463

AIC(coxEngineered) # 177400.4

# Accuracy
concordance(coxEngineered,newdata=datCredit_valid_TFD) # Concordance= 0.8111 se= 0.002782

# House keeping
rm(coxEngineered); gc()

#============================================================================================

modelVar <- c(modelVar,"slc_acct_pre_lim_perc_imputed_med","pmnt_method_grp")

#============================================================================================




#============================================================================================
# ------ 3. Interest Rate
varlist <- data.table(vars=c("InterestRate_Nom", "InterestRate_Margin"),
                      vartypes=c("prc","prc"))
#=========================================================================================



# ------ 3.1 How correlated are the two variables?

# - Correlation analysis
corrAnalysis(datCredit_train_TFD, varlist[vartypes!="cat"]$vars, corrThresh = 0.6, method = 'spearman') # Obtain correlation groups

### RESULTS: Correlation of 0.41, therefore no significant correlations.

# Goodness of fit test
csTable(datCredit_train_TFD, varlist$vars, TimeDef="TFD") # Insignificant difference.

aicTable(datCredit_train_TFD, varlist$vars, TimeDef="TFD", genPath=genObjPath) # [InterestRate_Nom] has the lowest AIC

# Accuracy test
concTable(datCredit_train_TFD, datCredit_valid_TFD,  varlist$vars, TimeDef="TFD", genPath=genObjPath) # [InterestRate_Nom] has the highest concordance

### CONCLUSION: Keep both variables in the model.



# ------ 3.2 How does [InterestRate_Margin] compare with [InterestRate_Margin_imputed_bin]?

# Initialize variables
vars <- c("InterestRate_Margin", "InterestRate_Margin_imputed_bin")

# Goodness of fit test
csTable(datCredit_train_TFD, vars, TimeDef="TFD") # No noticeable differences.

aicTable(datCredit_train_TFD, vars, TimeDef="TFD", genPath=genObjPath) # [InterestRate_Margin_imputed_bin] has the lowest AIC

# Accuracy test
concTable(datCredit_train_TFD, datCredit_valid_TFD, vars, TimeDef="TFD", genPath=genObjPath) # [InterestRate_Margin_imputed_bin] has the highest concordance

### CONCLUSIONS: Replace [InterestRate_Margin] with [InterestRate_Margin_imputed_bin]

varlist <- vecChange(varlist, Remove="InterestRate_Margin", Add=data.table(vars="InterestRate_Margin_imputed_bin", vartypes="cat"))



# ------ 3.3 How does [InterestRate_Margin_imputed_bin], InterestRate_Margin_Aggr_Med and [M_Repo_Rate] compare with one another?

vars <- c("InterestRate_Margin_imputed_bin", "M_Repo_Rate", "InterestRate_Margin_Aggr_Med")

# Goodness of fit test
csTable(datCredit_train_TFD, vars, TimeDef="TFD") # No noticable differences

aicTable(datCredit_train_TFD, vars, TimeDef="TFD", genPath=genObjPath) # [M_Repo_Rate] has the lowest AIC

# Accuracy test
concTable(datCredit_train_TFD, datCredit_valid_TFD, vars, TimeDef="TFD", genPath=genObjPath) # [InterestRate_Margin_Aggr_Med] has the highest concordance

### CONCLUSIONS: Keep [InterestRate_Margin_Aggr_Med] in the model.

varlist <- vecChange(varlist,Remove="InterestRate_Margin_imputed_bin", Add=data.table(vars="InterestRate_Margin_Aggr_Med", vartypes="prc"))



# ------ 3.4 Should we add lagging variables to the model?

vars <- c("InterestRate_Margin_Aggr_Med","InterestRate_Margin_Aggr_Med_1",
          "InterestRate_Margin_Aggr_Med_3","InterestRate_Margin_Aggr_Med_2",
          "InterestRate_Margin_Aggr_Med_9")

# Goodness of fit test
csTable(datCredit_train_TFD, vars, TimeDef="TFD") # No noticeable difference

aicTable(datCredit_train_TFD, vars, TimeDef="TFD", genPath=genObjPath) # [InterestRate_Margin_Aggr_Med] has the lowest AIC

# Accuracy test
concTable(datCredit_train_TFD, datCredit_valid_TFD, vars, TimeDef="TFD", genPath=genObjPath) # [InterestRate_Margin_Aggr_Med_3] has the highest concordance

### CONCLUSION: Add [InterestRate_Margin_Aggr_Med_3] to the model since it seems
###             to have the second best goodness of fit and concordance.

varlist <- vecChange(varlist,Add=data.table(vars=c("InterestRate_Margin_Aggr_Med_3"),
                                            vartypes=c("prc")))



# ------ 3.5 What is the performance of current thematic variables in univariate models?

vars <- c("InterestRate_Nom", "InterestRate_Margin_Aggr_Med", "InterestRate_Margin_Aggr_Med_3")

# Goodness of fit test
csTable(datCredit_train_TFD, vars, TimeDef="TFD") # Not a noticeable difference
#                     Variable      D
# 2   InterestRate_Margin_Aggr_Med 0.6540
# 3 InterestRate_Margin_Aggr_Med_3 0.6540
# 1               InterestRate_Nom 0.6539

aicTable(datCredit_train_TFD, vars, TimeDef="TFD", genPath=genObjPath) # [InterestRate_Nom] has the lowest AIC
#                     Variables      AIC
# 1:               InterestRate_Nom 188115.5
# 2:   InterestRate_Margin_Aggr_Med 188595.0
# 3: InterestRate_Margin_Aggr_Med_3 188678.2

# Accuracy test
concTable(datCredit_train_TFD, datCredit_valid_TFD, vars, TimeDef="TFD", genPath=genObjPath) # [InterestRate_Margin_Aggr_Med_3] has the highest concordance
#             Variables Concordance          SD LR_Statistic
# 1: InterestRate_Margin_Aggr_Med_3   0.5871809 0.004424536         1184
# 2:   InterestRate_Margin_Aggr_Med   0.5859953 0.004457969         1268
# 3:               InterestRate_Nom   0.5723056 0.004548748         1747

### CONCLUSION: Keep all variables in the model.



# ------ 3.6 What is the performance of current thematic cox ph model?

# Build cox model based on all thematic variables
coxInterest <- coxph(Surv(Start,End,Default_Ind) ~ InterestRate_Nom +
                           InterestRate_Margin_Aggr_Med + InterestRate_Margin_Aggr_Med_3,
                         id=LoanID, data=datCredit_train_TFD)
summary(coxInterest)

### RESULTS: All variables are significant.

# Test goodness of fit
GoF_CoxSnell_KS(coxInterest,datCredit_train_TFD,GraphInd = FALSE) # 0.6539

AIC(coxInterest) # 187743.1

# Accuracy
concordance(coxInterest, newdata=datCredit_valid_TFD)# Concordance= 0.5916 se= 0.004439
### RESULTS: Improved accuracy, but not particularly good.

# House keeping
rm(coxInterest)

#===========================================================================================

modelVar <- c(modelVar,"InterestRate_Nom", "InterestRate_Margin_Aggr_Med",
              "InterestRate_Margin_Aggr_Med_3")

#============================================================================================





#============================================================================================
# ------ 4. General
varlist <- data.table(vars=c("Balance","Instalment",
                             "Principal","LN_TPE","Term"),
                      vartypes=c("int", "fin", "fin","cat", "int"))
#=========================================================================================



# ------ 4.1 Which variables are highly correlated in the group?

# - Correlation analysis
corrAnalysis(datCredit_train_TFD, varlist[vartypes!="cat"]$vars, corrThresh = 0.6, method = 'spearman') # Obtain correlation groups

### RESULTS: 1) [Balance], [Installment] and [Principal] are highly correlated.



# ------ 4.2 Which variable in the correlated group should be kept in the model?

# Initialize variable
vars <- c("Balance", "Instalment", "Principal")

# Goodness of fit test
csTable(datCredit_train_TFD, vars, TimeDef="TFD") # No noticeable differences.

aicTable(datCredit_train_TFD, vars, TimeDef="TFD", genPath=genObjPath) # [Principal] has the lowest AIC

# Accuracy test
concTable(datCredit_train_TFD, datCredit_valid_TFD, vars, TimeDef="TFD", genPath=genObjPath) # [Principal] has the highest concordance

### CONCOLUSION:  [Principal] outperformed other variables, therefore it should be kept in the model.


# ------ 4.2.1 Can [Balance] and [Instalment] be replaced with [InstalmentToBalance_Aggr_Prop]?

# Initialize variable
vars <- c("Balance", "Instalment", "InstalmentToBalance_Aggr_Prop")

# Goodness of fit test of variables
aicTable(datCredit_train_TFD, vars, TimeDef="TFD", genPath=genObjPath) # [InstalmentToBalance_Aggr_Prop ] has the lowest AIC

# Accuracy test
concTable(datCredit_train_TFD, datCredit_valid_TFD, vars, TimeDef="TFD", genPath=genObjPath) # [Instalment] has the highest concordance

### CONCLUSION: [InstalmentToBalance_Aggr_Prop] cannot replace [Instalment] and [Balance].


# ------ 4.2.2 Can [Balance] be replaced with [ArrearsToBalance_Aggr_Prop]?

# Initialize variable
vars <- c("Balance", "ArrearsToBalance_Aggr_Prop")

# Goodness of fit test of variables
aicTable(datCredit_train_TFD, vars, TimeDef="TFD", genPath=genObjPath) # [ArrearsToBalance_Aggr_Prop] has the lowest AIC

# Accuracy test
concTable(datCredit_train_TFD, datCredit_valid_TFD, vars, TimeDef="TFD", genPath=genObjPath) # [Balance] has the highest concordance

### CONCOLUSION:  [ArrearsToBalance_Aggr_Prop] cannot reaplace [Balance]

### CONCLUSION: Remove [Instalment] and [Balance].

varlist <- vecChange(varlist,Remove=c("Instalment","Balance"))



# ------ 4.3 How does [Term] compare to [AgeToTerm_Aggr_Mean]?

# Initialize variables
vars <- c("Term", "AgeToTerm_Aggr_Mean")

# Goodness of fit
aicTable(datCredit_train_TFD, vars, TimeDef="TFD", genPath=genObjPath) # [AgeToTerm_Aggr_Mean] has the lowest AIC

# Accuracy test
concTable(datCredit_train_TFD, datCredit_valid_TFD, vars, TimeDef="TFD", genPath=genObjPath) # [AgeToTerm_Aggr_Mean] has the highest concordance

### CONCLUSION: Remove [Term].

varlist <- vecChange(varlist,Remove=c("Term"))



# ------ 4.4 What is the predictive power of the variables in univariate models?

vars <- c("Principal", "LN_TPE")

# Goodness of fit test
csTable(datCredit_train_TFD, vars, TimeDef="TFD") # No noticeable difference

aicTable(datCredit_train_TFD, vars, TimeDef="TFD", genPath=genObjPath) # [Principal] has the lowest AIC

# Accuracy test
concTable(datCredit_train_TFD, datCredit_valid_TFD, vars, TimeDef="TFD", genPath=genObjPath) # [Principal] has the highest concordance

### RESULTS:  [Principal] appears to have significant more concordance than the other variables.
###           [LN_TPE] could possibly be removed.



# ------ 4.5 What is the predictive power of the variables in a single model?

vars <- c("Principal", "LN_TPE")

# Goodness of fit
coxGen <- coxph(as.formula(paste0("Surv(Start,End,Default_Ind) ~ ",
                                  paste(vars,collapse=" + "))), id=LoanID,
                      datCredit_train_TFD)

summary(coxGen)# All variables are significant (p-value < 0.001)

GoF_CoxSnell_KS(coxGen,datCredit_train_TFD,GraphInd = FALSE) # 0.6534

AIC(coxGen) # 189418.9

# Accuracy
concordance(coxGen, newdata=datCredit_valid_TFD) # Concordance= 0.6159 se= 0.004391

# House keeping
rm(coxGen)

### CONCLUSION: Leave all variables in the model.


#============================================================================================

modelVar <- c(modelVar, "Principal", "LN_TPE")

#============================================================================================





#============================================================================================
# ------ 5. Portfolio Level
varlist <- data.table(vars=c("CuringEvents_Aggr_Prop","NewLoans_Aggr_Prop"),
                      vartypes=c("dec", "dec"))
#============================================================================================



# ------ 5.1 Are the values correlated?

# - Correlation analysis
corrAnalysis(datCredit_train_TFD, varlist[vartypes!="cat"]$vars, corrThresh = 0.6, method = 'spearman') # Obtain correlation groups

### RESULTS: The variables are not correlated significantly.



# ------ 5.2 What is the performance of the variables in univariate models?

# Initialize variable
vars <- c("CuringEvents_Aggr_Prop","NewLoans_Aggr_Prop")

# Goodness of fit test
csTable(datCredit_train_TFD, vars, TimeDef="TFD") # No noticeable difference

aicTable(datCredit_train_TFD, vars, TimeDef="TFD", genPath=genObjPath) # [NewLoans_Aggr_Prop] has the lowest AIC

# Accuracy test
concTable(datCredit_train_TFD, datCredit_valid_TFD, vars, TimeDef="TFD", genPath=genObjPath) # [NewLoans_Aggr_Prop] has the highest concordance

### CONCOLUSION:  Keep both variables in the model



# ------ 5.3 What is the performance of the variables in a single model?

# Goodness of fit
coxPort <- coxph(as.formula(paste0("Surv(Start,End,Default_Ind) ~ ",
                                        paste(vars,collapse=" + "))), id=LoanID,
                      datCredit_train_TFD)

summary(coxPort)
### RESULS: [CuringEvents_Aggr_Prop] has a high se(coef) ((19)realative to coef (-64))

### CONCLUSION: Remove [CuringEvents_Aggr_Prop]

varlist <- vecChange(varlist, Remove="CuringEvents_Aggr_Prop")



#============================================================================================

modelVar <- c(modelVar,"NewLoans_Aggr_Prop")

#============================================================================================





#============================================================================================
# ------ 6. Macro Economic
varlist <- data.table(vars=c("M_DTI_Growth","M_Emp_Growth","M_Inflation_Growth",
                             "M_RealGDP_Growth","M_RealIncome_Growth"),
                      vartypes=c("prc", "prc", "prc", "prc", "prc"))
#============================================================================================



# ------ 6.1 Which economic variables are highly correlated with one another?

# - Correlation analysis
corrAnalysis(datCredit_train_TFD, varlist$vars[varlist$vartypes != 'cat'],
             corrThresh = 0.6, method = 'spearman') # Obtain correlation groups

### RESULTS: [M_Emp_Growth], [M_RealGDP_Growth] and [M_RealIncome_Growth] are highly correlated



# ------ 6.2 Which correlated variable should remain in the model?

# Initialize variables
vars <- c("M_Emp_Growth", "M_RealGDP_Growth", "M_RealIncome_Growth")

# Goodness of fit test
csTable(datCredit_train_TFD, vars, TimeDef="TFD") # No noticeable difference

aicTable(datCredit_train_TFD, vars, TimeDef="TFD", genPath=genObjPath) # [M_RealGDP_Growth] has the lowest AIC

# Accuracy test
concTable(datCredit_train_TFD, datCredit_valid_TFD, vars, TimeDef="TFD", genPath=genObjPath) # [M_RealGDP_Growth] has the highest concordance

### DISCUSS CHANGE FROM [M_RealIncome_Growth] TO [M_RealGDP_Growth]   

### CONCLUSION: Keep [M_RealGDP_Growth] in the model.

varlist = vecChange(varlist, Remove=c("M_Emp_Growth", "M_RealIncome_Growth"))



# ------ 6.3 What are the predictive powers of the current variables in the models?

# Initialize variables
vars <- c("M_DTI_Growth", "M_Inflation_Growth", "M_RealGDP_Growth")

# Goodness of fit test
csTable(datCredit_train_TFD, vars, TimeDef="TFD") # No noticeable difference.

aicTable(datCredit_train_TFD, vars, TimeDef="TFD", genPath=genObjPath) # [] has the lowest AIC
#               Variables Iteration_1 Iteration_2 Iteration_3 Iteration_4 Iteration_5 Average
# 1:  M_Inflation_Growth      0.6479      0.6488      0.6452      0.6466      0.6503 0.64776
# 2: M_RealIncome_Growth      0.6448      0.6474      0.6470      0.6473      0.6471 0.64672
# 3:        M_DTI_Growth      0.6454      0.6441      0.6466      0.6484      0.6471 0.64632
# 4:               Range      0.0031      0.0060      0.0037      0.0018      0.0054 0.00400

### RESULTS:  After 5 iterations, there is too much variability in B-statistic, therefore no conclusions can be made.

# Accuracy test
concTable(datCredit_train_TFD, datCredit_valid_TFD, vars, TimeDef="TFD", genPath=genObjPath) # [] has the highest concordance
#               Variable Concordance          SD LR_Statistic
# 1:        M_DTI_Growth   0.5350355 0.003955393          117
# 2: M_RealIncome_Growth   0.5341913 0.003951199           68
# 3:  M_Inflation_Growth   0.5202318 0.004253188           55

### RESULTS: [M_DTI_Growth] has the highest concordance, albeit not that high.
### CONCLUSION: The low predictive power may be enhanced with lags.



# ------ 6.4 Which lags for current variables be included in the models?


# ------ 6.4.1 Which lags for [M_DTI_Growth] be included in the models?
# Initialize variables
vars <- c("M_DTI_Growth_1","M_DTI_Growth_12","M_DTI_Growth_3", "M_DTI_Growth_6",
          "M_DTI_Growth_9","M_DTI_Growth")

# Goodness of fit test
csTable(datCredit_train_TFD, vars, TimeDef="TFD") # No noticeable difference

aicTable(datCredit_train_TFD, vars, TimeDef="TFD", genPath=genObjPath) # [M_DTI_Growth] has the lowest AIC

# Accuracy test
concTable(datCredit_train_TFD, datCredit_valid_TFD, vars, TimeDef="TFD", genPath=genObjPath)

### RESULTS: Only [M_DTI_Growth_12], [M_DTI_Growth_9] and [M_DTI_Growth_6] appears to have the best fit.

### CONCLUSION: Include the following lags in the model: [M_DTI_Growth_6] and [M_DTI_Growth_9]

varlist <- vecChange(varlist,Add=data.table(vars=c("M_DTI_Growth_6", "M_DTI_Growth_9"),
                                            vartypes=c("prc", "prc")))


# ------ 6.4.2 Should lags for [M_Inflation_Growth] be included in the models?
# Initialize variables
vars <- c("M_Inflation_Growth_1","M_Inflation_Growth_12","M_Inflation_Growth_3",
          "M_Inflation_Growth_6","M_Inflation_Growth_9","M_Inflation_Growth")

# Goodness of fit test
csTable(datCredit_train_TFD, vars, TimeDef="TFD") # No noticeable difference

aicTable(datCredit_train_TFD, vars, TimeDef="TFD", genPath=genObjPath)

### RESULTS: [M_Inflation_Growth_6], [M_Inflation_Growth_1] and [M_Inflation_Growth_3] appear to have the best fits.

# Accuracy test
concTable(datCredit_train_TFD, datCredit_valid_TFD, vars, TimeDef="TFD", genPath=genObjPath)

### RESULTS: [M_Inflation_Growth_12], [M_Inflation_Growth_6] and [M_Inflation_Growth_9] appear to be the most accurate.

### CONCLUSION: Add [M_Inflation_Growth_6] to the model since it seems
###             to have the first best goodness of fit and second best accuracy

varlist <- vecChange(varlist,Add=data.table(vars=c("M_Inflation_Growth_6"),
                                            vartypes=c("prc")))


# ------ 6.4.3 Should lags for [M_RealIncome_Growth] be included in the models?
vars <- c("M_RealIncome_Growth_1","M_RealIncome_Growth_12","M_RealIncome_Growth_3",
          "M_RealIncome_Growth_6","M_RealIncome_Growth_9","M_RealIncome_Growth")

# Goodness of fit test
csTable(datCredit_train_TFD, vars, TimeDef="TFD") # No noticeable difference

aicTable(datCredit_train_TFD, vars, TimeDef="TFD", genPath=genObjPath)

### RESULTS: [M_RealIncome_Growth_12], [M_RealIncome_Growth_9] and [M_RealIncome_Growth_6] appear to have the best fits.

# Accuracy test
concTable(datCredit_train_TFD, datCredit_valid_TFD, vars, TimeDef="TFD", genPath=genObjPath)

### RESULTS: [M_RealIncome_Growth_12], [M_RealIncome_Growth_9] and [M_RealIncome_Growth_6] appear to be the most accurate.

### CONCLUSION: Do not include the lags for [M_RealIncome_Growth].



# ------ 6.5 What is the predictive performance of the current univariate thematic models?

# Initialize macro-economic thematic variables
vars <- c("M_DTI_Growth", "M_Inflation_Growth","M_DTI_Growth_6", "M_DTI_Growth_9",
          "M_Inflation_Growth_6", "M_Inflation_Growth_12", "M_RealIncome_Growth")

# Goodness of fit test
csTable(datCredit_train_TFD, vars, TimeDef="TFD") # No noticeable differences

aicTable(datCredit_train_TFD, vars, TimeDef="TFD", genPath=genObjPath)

### RESULTS: [M_DTI_Growth], [M_DTI_Growth_9] and [M_DTI_Growth_6] appear to be the most accurate.

# Accuracy test
concTable(datCredit_train_TFD, datCredit_valid_TFD, vars, TimeDef="TFD", genPath=genObjPath)

### RESULTS: [M_Inflation_Growth_6], [M_DTI_Growth_9] and [M_DTI_Growth_6] appear to be the most accurate.

### RESULTS: All variables should be kept in the model. Although [M_RealIncome_Growth]
###           has a low concordance, it is kept in the model as a proxy for [M_Emp_Growth].



# ------ 6.6 What is the predictive performance of the current thematic model?
# Initialize macro-economic thematic variables
vars <- c("M_DTI_Growth", "M_Inflation_Growth","M_DTI_Growth_6", "M_DTI_Growth_9",
          "M_Inflation_Growth_6", "M_RealIncome_Growth", "M_Inflation_Growth_12")

# Build model based on macro economic variables
coxMacro <- coxph(as.formula(paste0("Surv(Start,End,Default_Ind) ~ ",
                                    paste(vars,collapse=" + "))), id=LoanID,
                  datCredit_train_TFD)
summary(coxMacro)

### RESULTS: [M_DTI_Growth_6] and [M_DTI_Growth_9] have opposite coefficients, which
###           seems counter intuitive, therefore they should be removed. Similar case 
###           can be concluded for [M_Inflation_Growth_12].

vars <- c("M_DTI_Growth", "M_Inflation_Growth", "M_Inflation_Growth_6", "M_RealIncome_Growth")

# Build model based on macro economic variables
coxMacro <- coxph(as.formula(paste0("Surv(Start,End,Default_Ind) ~ ",
                                    paste(vars,collapse=" + "))), id=LoanID,
                  datCredit_train_TFD)
summary(coxMacro)

### RESULTS:  All values have significant p-value.

# Test goodness of fit
GoF_CoxSnell_KS(coxMacro,datCredit_train_TFD,GraphInd = FALSE) # 0.6523

AIC(coxMacro) # 187854.8

# Test accuracy
concordance(coxMacro, newdata=datCredit_valid_TFD) # Concordance= 0.5876 se= 0.004408



#============================================================================================

modelVar <- c(modelVar,"M_DTI_Growth", "M_Inflation_Growth", "M_Inflation_Growth_6",
              "M_RealIncome_Growth")

#============================================================================================





# ------ 7. Check semi-final model for statistical significance and parsimony

# - Initialize variables
vars2 <- c("g0_Delinq_SD_4", "Arrears", "g0_Delinq_Ave", "TimeInDelinqState_Lag_1",      
           "slc_acct_arr_dir_3_Change_Ind", "slc_acct_roll_ever_24_imputed_mean",
           "slc_acct_pre_lim_perc_imputed_med","InterestRate_Margin_Aggr_Med",
           "InterestRate_Margin_Aggr_Med_3","M_DTI_Growth","M_Inflation_Growth",
           "M_Inflation_Growth_6", "M_RealIncome_Growth","pmnt_method_grp",
           "InterestRate_Nom", "Principal","LN_TPE", "NewLoans_Aggr_Prop")

# - Build model based on variables
cox_TFD <- coxph(as.formula(paste0("Surv(Start,End,Default_Ind) ~ ", paste(vars2,collapse=" + "))),
                 id=LoanID, datCredit_train_TFD, ties="efron")
# NOTE: Default option for handling ties is recently "Efron's method", but we specify this for backwards compatability
summary(cox_TFD)
### RESULTS: [TimeInDelinqState_Lag_1], [InterestRate_Margin_Aggr_Med], [InterestRate_Margin_Aggr_Med_3],
###           [M_Inflation_Growth], [M_Inflation_Growth_6], [LN_TPEWHL], [NewLoans_Aggr_Prop] 
###           have insignificant p-value (> 0.05).

vars2 <- vars2[!(vars2 %in% c("TimeInDelinqState_Lag_1", "InterestRate_Margin_Aggr_Med",
                            "InterestRate_Margin_Aggr_Med_3","M_Inflation_Growth",
                            "M_Inflation_Growth_6", "LN_TPE", "NewLoans_Aggr_Prop"))]
# - Build model based on variables
cox_TFD <- coxph(as.formula(paste0("Surv(Start,End,Default_Ind) ~ ", paste(vars2,collapse=" + "))),
                 id=LoanID, datCredit_train_TFD, ties="efron")
summary(cox_TFD)
### RESULTS: All variables are significant

#============================================================================================





# ------ 8. Final model

# - Confirm prepared datasets are loaded into memory
if (!exists('datCredit_train_TFD')) unpack.ffdf(paste0(genPath,"creditdata_train_TFD"), tempPath);gc()
if (!exists('datCredit_valid_TFD')) unpack.ffdf(paste0(genPath,"creditdata_valid_TFD"), tempPath);gc()

# - Initialize variables
vars2 <- c("g0_Delinq_SD_4", "Arrears", "g0_Delinq_Ave",      
           "slc_acct_arr_dir_3_Change_Ind", "slc_acct_roll_ever_24_imputed_mean",
           "slc_acct_pre_lim_perc_imputed_med", "M_DTI_Growth",
           "M_RealIncome_Growth","pmnt_method_grp",
           "InterestRate_Nom", "Principal")

# - Build model based on variables
cox_TFD <- coxph(as.formula(paste0("Surv(Start,End,Default_Ind) ~ ", paste(vars2,collapse=" + "))),
                 id=LoanID, datCredit_train_TFD, ties="efron")
summary(cox_TFD); AIC(cox_TFD); concordance(cox_TFD)
### RESULTS: AIC: 92825.18 Harrell's c: 0.9971

c <- coefficients(cox_TFD)
(c <- data.table(Variable=names(c),Coefficient=c))
#                           Variable      Coefficient
# 1:                     g0_Delinq_SD_4  6.7992534010147
# 2:                            Arrears  0.0000094732796
# 3:                      g0_Delinq_Ave -4.7616827115497
# 4:      slc_acct_arr_dir_3_Change_Ind  0.7100710046257
# 5: slc_acct_roll_ever_24_imputed_mean  0.8471021545785
# 6:  slc_acct_pre_lim_perc_imputed_med -2.7530111775470
# 7:                       M_DTI_Growth  1.8646494414421
# 8:                M_RealIncome_Growth  4.1549018668866
# 9:        pmnt_method_grpMISSING_DATA  1.4697040486667
# 10:     pmnt_method_grpSalary/Suspense  0.6942347599052
# 11:           pmnt_method_grpStatement  0.3341055117906
# 12:                   InterestRate_Nom  5.3517761028485
# 13:                          Principal -0.0000001306697


# -Test Goodness of fit using bootstrapped B-statistics (1-KS statistic) over single-factor models
csTable_TFD <- csTable(datCredit_train_TFD,vars2, TimeDef="TFD", seedVal=1, numIt=10,
                       fldSpellID="PerfSpell_Key", fldLstRowInd="PerfSpell_Exit_Ind", fldEventInd="Default_Ind", genPath=genPath)
### RESULTS: Top 3 single-factor models: M_RealIncome_Growth + g0_Delinq_Ave  + slc_acct_pre_lim_perc_imputed_med  

aicTable_TFD <- aicTable(datCredit_train_TFD, variables=vars2,fldSpellID="PerfSpell_Key",
                         TimeDef="TFD", genPath=genPath)
### RESULTS: Top 3 single-factor models: g0_Delinq_SD_4 + slc_acct_roll_ever_24_imputed_mean + slc_acct_arr_dir_3_Change_Ind  
# Where the first three results have AIC values significantly different from the rest.

# Test accuracy using Harrell's c-statistic over single-factor models
concTable_TFD <- concTable(datCredit_train_TFD, datCredit_valid_TFD, variables=vars2, 
                           fldSpellID="PerfSpell_Key", TimeDef="TFD", genPath=genPath)
### RESULTS: Top x single-factor models (>80%):
# Arrears + PerfSpell_g0_Delinq_Num + TimeInDelinqState_Lag_1 + slc_acct_arr_dir_3_Change_Ind  

# - Combine results into a single object
Table_TFD <- concTable_TFD[,1:2] %>% left_join(aicTable_TFD, by ="Variable") %>% left_join(data.table(csTable_TFD$Results), by="Variable")

# - Test Goodnes-of-fit using Cox-Snell, having measured distance between residual distribution and unit exponential using KS-statistic
GoF_CoxSnell_KS(cox_TFD,datCredit_train_TFD, GraphInd=TRUE, legPos=c(0.6,0.4), panelTitle="Time to First Default (TFD) model",
                fileName = paste0(genFigPath, "TFD/KS_Test_CoxSnellResiduals_Exp_TFD", ".png"), dpi=280) # 0.6167
### RESULTS: Goodness of fit for the model seems to be a bit low.


concordance(cox_TFD) # Concordance= 0.9971 se= 0.0001901

# Save objects
pack.ffdf(paste0(genObjPath,"TFD_Univariate_Models"), Table_TFD)
pack.ffdf(paste0(genPath,"TFD_Cox_Model"), cox_TFD)

# - Cleanup
rm(datCredit_train_TFD, datCredit_valid_TFD, cox_TFD, c); gc()