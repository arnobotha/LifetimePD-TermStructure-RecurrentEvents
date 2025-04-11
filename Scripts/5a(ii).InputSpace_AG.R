# ============================================ INPUT SPACE: AG =========================================
# Divide data into thematic groups and perform data analysis on them to compile an input space for the AG model
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
#   - 3c(ii).Data_Fusion2_AG.R

# -- Inputs:
#   - datCredit_train_AG | Prepared from script 3c
#   - datCredit_valid_AG | Prepared from script 3c

#
# -- Outputs:
#   - Input_Space
# ------------------------------------------------------------------------------------------------------

# ------ 1. Preliminaries
# - Confirm prepared datasets are loaded into memory
if (!exists('datCredit_train_AG')) unpack.ffdf(paste0(genPath,"creditdata_train_AG"), tempPath);gc()
if (!exists('datCredit_valid_AG')) unpack.ffdf(paste0(genPath,"creditdata_valid_AG"), tempPath);gc()


#============================================================================================
# ------ 1. Delinquency measures
varlist <- data.table(vars=c("PerfSpell_g0_Delinq_Num","Arrears" ,"TimeInDelinqState_Lag_1",
                             "g0_Delinq_Any_Aggr_Prop","g0_Delinq_Ave"),
                      vartypes=c("acc", "acc", "acc", "dte", "dte"))

#==========================================================================================



# ------ 1.1 Which time window length is the best in calculating Delinquency volatility?

# Initialize variables to be tested
vars <- c("g0_Delinq_SD_4", "g0_Delinq_SD_5", "g0_Delinq_SD_6", "g0_Delinq_SD_9", "g0_Delinq_SD_12")

# Goodness of fit test
csTable(datCredit_train_AG, vars, TimeDef="AG") # [g0_Delinq_SD_12] has the lowest KS-statistic

aicTable(datCredit_train_AG, vars, TimeDef="AG", genPath=genObjPath) # [g0_Delinq_SD_5] has the lowest AIC

# Accuracy test
concTable(datCredit_train_AG, datCredit_valid_AG, vars, TimeDef="AG", genPath=genObjPath) # [g0_Delinq_SD_4] has the highest concordance

### Conclusion: Larger window are less influenced by large changes therefore significant changes
###             are less pronounced. Include [g0_Delinq_SD_4] in the varlist.

varlist <- vecChange(varlist,Remove="PerfSpell_g0_Delinq_SD",Add=data.table(vars=c("g0_Delinq_SD_4"), vartypes=c("acc")))

# ------ 1.2 Which variables are highly correlated?

# Correlation analysis
corrAnalysis(datCredit_train_AG, varlist[vartypes!="cat"]$vars, corrThresh = 0.6, method = 'spearman') # Obtain correlation groups

### RESULTS:  [g0_Delinq_Any_Aggr_Prop] and [g0_Delinq_Ave] with a correlation of 100%

### CONCLUSION: A single variable from each group must be retained while the rest are removed.



# ------ 1.2.1 Which variable should be kept from group 1) [g0_Delinq_Any_Aggr_Prop] and [g0_Delinq_Ave]

vars <- c("g0_Delinq_Any_Aggr_Prop", "g0_Delinq_Ave")

# Goodness of fit
csTable(datCredit_train_AG, vars, TimeDef="AG") # No noticeable difference

aicTable(datCredit_train_AG, vars, TimeDef="AG", genPath=genObjPath) # [g0_Delinq_Ave] has the lowest AIC

### RESULTS: [g0_Delinq_Ave] seems to have the better goodness of fit.

# Accuracy
concTable(datCredit_train_AG, datCredit_valid_AG, vars, TimeDef="AG", genPath=NA) # [g0_Delinq_Ave] has the highest concordance

### CONCLUSION: [g0-Delinq_Ave] seems to outperform [g0_Delinq_Any_Aggr_Prop] and therefore is kept in the model.

varlist <- vecChange(varlist,Remove=c("g0_Delinq_Any_Aggr_Prop"))

# ------ 1.3 How does [PerfSpell_g0_Delinq_Num] compare to [g0_Delinq_SD_4]?

# Initialize variables
vars <- c("PerfSpell_g0_Delinq_Num", "g0_Delinq_SD_4")

# Goodness of fit
csTable(datCredit_train_AG, vars, TimeDef="AG") # [PerfSpell_g0_Delinq_Num] has the lowest KS-statistic

aicTable(datCredit_train_AG, vars, TimeDef="AG", genPath=genObjPath) # [g0_Delinq_SD_4] has the lowest AIC

# Accuracy
concTable(datCredit_train_AG, datCredit_valid_AG, vars, TimeDef="AG", genPath=genObjPath) # [g0_Delinq_SD_4] has the highest concordance

### CONCLUSION: Keep [g0_Delinq_SD_4] in the model and remove [PerfSpell_g0_Delinq_Num]

varlist <- vecChange(varlist,Remove="PerfSpell_g0_Delinq_Num")



# ------ 1.4 What is the performance of current thematic variables in univariate models?

# Build thematic model based on remaining delinquency variables.
vars <- c("PerfSpell_g0_Delinq_Num","TimeInDelinqState","g0_Delinq_Ave",
          "slc_acct_arr_dir_3", "Arrears")

# Goodness of fit test
csTable(datCredit_train_AG, vars, TimeDef="AG")

### NOTE: [TimeInDelinqState] ran out of iterations and did not converge
### NOTE: [slc_acct_arr_dir_3] exp overflow due to covariates

# ------ 1.4.1 Why does [TimeInDelinqState] not converge?

cox <- coxph(Surv(Start,End,Default_Ind) ~ TimeInDelinqState, id=LoanID, datCredit_train_AG)
### RESULTS:  Beta tends to Inf. After some inspection on the data it relates to all
###           Defaulting events starting in a new delinquency state,
###           i.e. [TimeInDelinqState] = 1. Therefore quasi-complete separation
###           seems to be present.

### INVESTIGATE: WHETHER QUASI-COMPLETE SEPERATION IS PRESENT
datCredit_train_AG[TimeInDelinqState!=1 & Default_Ind==1, .N]
### RESULTS: 0
datCredit_train_AG[TimeInDelinqState==1 & Default_Ind==0, .N]
### RESULTS: 169272

### CONCLUSION: Quasi-complete separation seems to be present for [g0_Delinq]=3
###             and should therefore be removed.

# Test a lagged version of TimeInDelinqState
cox <- coxph(Surv(Start,End,Default_Ind) ~ TimeInDelinqState_Lag_1, id=LoanID,
             datCredit_train_AG)
summary(cox); rm(cox)

### CONCLUSION: Replace old variable with new variable, since it fits a stable model.

varlist <- vecChange(varlist,Remove="TimeInDelinqState",
                     Add=data.table(vars="TimeInDelinqState_Lag_1",vartypes="acc"))


# ------ 1.4.2 Why does [slc_acct_arr_dir_3] exp overflow?

describe(datCredit_train_AG[,slc_acct_arr_dir_3])
# Value            CURING MISSING_DATA      ROLLING         SAME
# Proportion        0.014        0.125        0.022        0.838

describe(datCredit_train_AG[Default_Ind==1,slc_acct_arr_dir_3])
# Value            CURING MISSING_DATA      ROLLING         SAME
# Proportion        0.004        0.157        0.788        0.051

### RESULTS:  Although the ROLLING level is not dominant in datCredit_train_AG,
###           we can clearly see that it is highly predictive of a default event
###           occurring, given its high proportion if Default_Ind == 1.

# Use an indicator variable for a change in account
cox <- coxph(Surv(Start,End,Default_Ind) ~ slc_acct_arr_dir_3_Change_Ind,
             datCredit_valid_AG, id=LoanID)
summary(cox); rm(cox)

### CONCLUSION: Replace old variable with new variable, since resulting model is stable.

varlist <- vecChange(varlist,Remove=c("slc_acct_arr_dir_3") ,
                     Add=data.table(vars=c("slc_acct_arr_dir_3_Change_Ind"),
                                    vartypes=c("bin")))



# ------ 1.4.3 What is the performance of current thematic variables in univariate models?

vars <- c("PerfSpell_g0_Delinq_Num","g0_Delinq_Ave", "Arrears",
          "TimeInDelinqState_Lag_1","slc_acct_arr_dir_3_Change_Ind")

# Goodness of fit test
csTable(datCredit_train_AG, vars, TimeDef="AG")

#                        Variable      D
# 2                 g0_Delinq_Ave 0.6540
# 3                       Arrears 0.6539
# 1       PerfSpell_g0_Delinq_Num 0.6520
# 5 slc_acct_arr_dir_3_Change_Ind 0.6514
# 4       TimeInDelinqState_Lag_1 0.6445

aicTable(datCredit_train_AG, vars, TimeDef="AG", genPath=genObjPath)

#                   Variables      AIC
# 1:       TimeInDelinqState_Lag_1 156436.3
# 2: slc_acct_arr_dir_3_Change_Ind 160294.1
# 3:       PerfSpell_g0_Delinq_Num 182572.5
# 4:                       Arrears 183412.2
# 5:                 g0_Delinq_Ave 187989.6

# Accuracy test
concTable(datCredit_train_AG, datCredit_valid_AG, vars, TimeDef="AG", genPath=genObjPath) # [] has the highest concordance

#                       Variables Concordance           SD LR_Statistic
# 1:                       Arrears   0.9606159 0.0016931039         6450
# 2:       PerfSpell_g0_Delinq_Num   0.9468068 0.0006694657         7290
# 3:       TimeInDelinqState_Lag_1   0.9020457 0.0029892384        33426
# 4: slc_acct_arr_dir_3_Change_Ind   0.8770983 0.0013364884        29569
# 5:                 g0_Delinq_Ave   0.5837088 0.0042886285         1873

### CONCLUSION: Leave all variables in the model.



# ------ 1.5 What is the performance of current thematic cox ph model?

# Build cox model based on all thematic variables
coxDelinq <- coxph(Surv(Start,End,Default_Ind) ~ g0_Delinq_Ave +
                     PerfSpell_g0_Delinq_Num + slc_acct_arr_dir_3_Change_Ind +
                     TimeInDelinqState_Lag_1 + Arrears, id=LoanID,
                   data=datCredit_train_AG)

summary(coxDelinq)

### RESULTS: All variables are significant (p-value < 0.001)

# Test goodness of fit
GoF_CoxSnell_KS(coxDelinq,datCredit_train_AG,GraphInd = FALSE) # 0.6451

# Test accuracy
concordance(coxDelinq, newdata=datCredit_valid_AG) # Concordance= 0.9638 se= 0.001422

### CONCLUSION: Keep all variables in the model

# House keeping
rm(coxDelinq)

#============================================================================================

modelVar <- c("PerfSpell_g0_Delinq_Num", "Arrears", "g0_Delinq_Ave",
              "TimeInDelinqState_Lag_1", "slc_acct_arr_dir_3_Change_Ind")

#============================================================================================
# ------ 2. Engineered measures
varlist <- data.table(vars=c("slc_acct_pre_lim_perc_imputed_med",
                             "slc_acct_prepaid_perc_dir_12_imputed_med",
                             "pmnt_method_grp") ,
                      vartypes=c("prc", "dec", "cat"))

### HW: Remove pmnt_method_grp

#=========================================================================================

# ------ 2.1 Which variables should be removed due to high correlation?

# - Correlation analysis
corGroups <- corrAnalysis(datCredit_train_AG, varlist[vartypes!="cat"]$vars, corrThresh = 0.6, method = 'spearman') # Obtain correlation groups

### RESULTS:  1) [slc_acct_pre_lim_perc_imputed_med] and [slc_acct_prepaid_perc_dir_12_imputed_med] with a correlation of 0.82

# Initialize variables to be tested
vars <- c("slc_acct_pre_lim_perc_imputed_med", "slc_acct_prepaid_perc_dir_12_imputed_med")

# Goodness of fit test
csTable(datCredit_train_AG, vars, TimeDef="AG") # No noticeable difference found

aicTable(datCredit_train_AG, vars, TimeDef="AG", genPath=genObjPath) # [slc_acct_pre_lim_perc_imputed_med] has the lowest AIC

# Accuracy test
concTable(datCredit_train_AG, datCredit_valid_AG, vars, TimeDef="AG", genPath=genObjPath) # [slc_acct_pre_lim_perc_imputed_med] has the highest concordance

varlist <- vecChange(varlist,Remove="slc_acct_prepaid_perc_dir_12_imputed_med")

# ------ 2.2 What is the predictive power of the variables left in varlist?

# Initialize variables to be tested
vars <- c("pmnt_method_grp", "slc_acct_pre_lim_perc_imputed_med")

csTable(datCredit_train_AG, vars, TimeDef="AG") # [slc_acct_pre_lim_perc_imputed_med ] has the lowest KS-statistic

aicTable(datCredit_train_AG, vars, TimeDef="AG", genPath=genObjPath) # [pmnt_method_grp] has the lowest AIC

# Accuracy test
concTable(datCredit_train_AG, datCredit_valid_AG, vars, TimeDef="AG", genPath=genObjPath) # [pmnt_method_grp] has the highest concordance

### Conclusion: All values should be kept in the model



# ------ 2.3 How predictive is a single model based on the thematic variables?

# Initialize variables
vars <- c("pmnt_method_grp", "slc_acct_pre_lim_perc_imputed_med")

# Build Cox model based on each variable
coxEngineered <- coxph(Surv(Start, End, Default_Ind) ~ pmnt_method_grp +
                         slc_acct_pre_lim_perc_imputed_med, id=LoanID,
                       data=datCredit_train_AG)
summary(coxEngineered)

### RESULTS: [slc_acct_pre_lim_perc_imputed_med]  has a high coef of 59

# Goodness of fit
GoF_CoxSnell_KS(coxEngineered,datCredit_train_AG, GraphInd=FALSE) # 0.6446

# Accuracy
concordance(coxEngineered,newdata=datCredit_valid_AG) # Concordance= 0.8111 se= 0.002782

# House keeping
rm(coxEngineered); gc()

#============================================================================================

modelVar <- c(modelVar,"slc_acct_pre_lim_perc_imputed_med","pmnt_method_grp")

#============================================================================================



#============================================================================================
# ------ 7. Final Model

# - Confirm prepared datasets are loaded into memory
if (!exists('datCredit_train_AG')) unpack.ffdf(paste0(genPath,"creditdata_train_AG"), tempPath);gc()
if (!exists('datCredit_valid_AG')) unpack.ffdf(paste0(genPath,"creditdata_valid_AG"), tempPath);gc()

# - Initialize variables | AB-variant
vars2 <- c("PerfSpell_g0_Delinq_Num", "Arrears", "g0_Delinq_Ave", "TimeInDelinqState_Lag_1",      
           "slc_acct_arr_dir_3_Change_Ind", "slc_acct_pre_lim_perc_imputed_med", 
           "InterestRate_Margin_Aggr_Med", "InterestRate_Margin_Aggr_Med_3", 
           "M_DTI_Growth","M_Inflation_Growth", "M_Inflation_Growth_6", "M_RealIncome_Growth")

# - Build model based on variables
cox_AG <- coxph(as.formula(paste0("Surv(Start,End,Default_Ind) ~ ", paste(vars2,collapse=" + "))),
                 id=LoanID, datCredit_train_AG, ties="efron")
# NOTE: Default option for handling ties is recently "Efron's method", but we specify this for backwards compatability
summary(cox_AG)
### RESULTS: All variables are statistically significant

c <- coefficients(cox_AG)
(c <- data.table(Variable=names(c),Coefficient=c))
#                           Variable      Coefficient
#1:           PerfSpell_g0_Delinq_Num    0.01088066805
#2:                           Arrears    0.00002245914
#3:                     g0_Delinq_Ave   -6.29712945291
#4:           TimeInDelinqState_Lag_1   -0.16289352028
#5:     slc_acct_arr_dir_3_Change_Ind    3.16314704845
#6: slc_acct_pre_lim_perc_imputed_med   -8.77841646561
#7:      InterestRate_Margin_Aggr_Med  250.40823968411
#8:    InterestRate_Margin_Aggr_Med_3 -333.82470866537
#9:                      M_DTI_Growth   -5.45866171416
#10:                M_Inflation_Growth   10.26239576449
#11:              M_Inflation_Growth_6    6.86623111352
#12:               M_RealIncome_Growth  -17.83127157512

# -T Test Goodness of fit using bootstrapped B-statistics (1-KS statistic) over single-factor models
csTable_AG <- csTable(datCredit_train_AG,vars2, TimeDef="AG", seedVal=1, numIt=10,
                       fldSpellID="PerfSpell_Key", fldLstRowInd="PerfSpell_Exit_Ind", fldEventInd="Default_Ind", genPath=genPath)
### RESULTS: Top 3 single-factor models: M_RealIncome_Growth  + M_Inflation_Growth  + Arrears  
# Results are, however, very close one another such that meaningful analysis is not possible

# Test accuracy using Harrell's c-statistic over single-factor models
concTable_AG <- concTable(datCredit_train_AG, datCredit_valid_AG, variables=vars2, 
                           fldSpellID="PerfSpell_Key", TimeDef="AG", genPath=genPath)
### RESULTS: Top x single-factor models (>80%):
# Arrears + PerfSpell_g0_Delinq_Num + TimeInDelinqState_Lag_1 + slc_acct_arr_dir_3_Change_Ind  

# - Combine results into a single object
Table_AG <- left_join(csTable_AG$Results, concTable_AG, by="Variable")

# - Test Goodnes-of-fit using Cox-Snell, having measured distance between residual distribution and unit exponential using KS-statistic
GoF_CoxSnell_KS(cox_AG,datCredit_train_AG, GraphInd=TRUE, legPos=c(0.6,0.4), panelTitle="Andersen-Gill (AG) model",
                fileName = paste0(genFigPath, "AG/KS_Test_CoxSnellResiduals_Exp_AG", ".png"), dpi=280) # 0.6372
### RESULTS: Goodness of fit for the model seems to be a bit low.

# Save objects
pack.ffdf(paste0(genObjPath,"AG_Univariate_Models"), Table_AG)
pack.ffdf(paste0(genPath,"AG_Cox_Model"), cox_AG)

# - Cleanup
rm(datCredit_train_AG, datCredit_valid_AG, cox_AG, c); gc()
