# ============================================ INPUT SPACE =========================================
# Divide data into thematic groups and perform data analysis on them to compile an input space for the TTFD model
# ------------------------------------------------------------------------------------------------------
# PROJECT TITLE: Default survival modelling
# SCRIPT AUTHOR(S): Bernard Scheepers
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
#   - datCredit_train_TFD | Prepared from script 3b
#   - datCredit_valid_TFD | Prepared from script 3b

#
# -- Outputs:
#   - Input_Space
# ------------------------------------------------------------------------------------------------------

# ------ 1. Preliminaries
# - Confirm prepared datasets are loaded into memory
if (!exists('datCredit_train_TFD')) unpack.ffdf(paste0(genPath,"creditdata_train_TFD"), tempPath);gc()
if (!exists('datCredit_valid_TFD')) unpack.ffdf(paste0(genPath,"creditdata_valid_TFD"), tempPath);gc()

# BS: will work these changes into script c once master's is done.

# BS: Detected a problem with [PerfSpell_Exit_Ind] where some loans never leave state (although they do) (Problem lies with PerfSpell_Max_Date)
datCredit_train_TFD[PerfSpell_Exit_Ind==1,.N]
# 62132
datCredit_train_TFD[!duplicated(PerfSpell_Key),.N]
# 62723

# ### HOMEWORK: Create [Removed]
datCredit_train_TFD$Removed <- with(datCredit_train_TFD, ave(seq_along(PerfSpell_Key), PerfSpell_Key, FUN = function(x) x == max(x)))
datCredit_valid_TFD$Removed <- with(datCredit_valid_TFD, ave(seq_along(PerfSpell_Key), PerfSpell_Key, FUN = function(x) x == max(x)))

datCredit_train_TFD[Removed==1,.N]
# 62723

### HOMEWORK: Create [slc_acct_arr_dir_3_Change_Ind]
datCredit_train_TFD[, slc_acct_arr_dir_3_Change_Ind := ifelse(slc_acct_arr_dir_3 != "SAME", 1,0)]
datCredit_valid_TFD[, slc_acct_arr_dir_3_Change_Ind := ifelse(slc_acct_arr_dir_3 != "SAME", 1,0)]

# ### AB: This [Removed]-variable already exists as [PerfSpell_Exit_Ind], created in script 3a(i).. Remove removed from all subsequent scripts. ;-)
# datCredit_train_TFD[, slc_acct_arr_dir_3_Change_Ind := ifelse(slc_acct_arr_dir_3 != "SAME", 1,0)]
# ### AB: At this point, we really don't want to perform data preparation and feature engineering within the modelling scripts.
# # By doing so, we are breaking our own "mould"/pattern. Move these into the appropriate section within script 3c.
# 
# datCredit_valid_TFD[,Removed := ifelse(Date==PerfSpell_Max_Date,T,F)]
# datCredit_valid_TFD[, slc_acct_arr_dir_3_Change_Ind := ifelse(slc_acct_arr_dir_3 != "SAME", 1,0)]
# datCredit_valid_TFD[,TimeInDelinqState_Lag_1 := shift(TimeInDelinqState,fill=0),by=LoanID]
# ### AB: Same comment here as before
# 
# datCredit_train_TFD[,slc_acct_roll_ever_24_imputed_med_f := factor(slc_acct_roll_ever_24_imputed_med)]
# ### AB: This variable is created but seemingly never used? Rather remote if true.
# 
# ### AB [2024-12-01]: I moved the function definitions to the most appropriate script, i.e., 0b(i).
# # Please go and craft decent comments for them, i.e., inputs and outputs, as one would do with any function

#============================================================================================
# ------ 1. Delinquency measures
varlist <- data.table(vars=c("g0_Delinq","g0_Delinq_fac","PerfSpell_g0_Delinq_Num",
                             "Arrears" ,"TimeInDelinqState","g0_Delinq_Any_Aggr_Prop",
                             "g0_Delinq_Ave","slc_acct_arr_dir_3",
                             "slc_acct_roll_ever_24_imputed_mean"),
                      vartypes=c("acc", "cat", "acc", "acc", "acc", "dte", "dte",
                                 "cat", "acc"))

#=========================================================================================



# ------ 1.1 Which time window length is the best in calculating Delinquency volatility?

# Initialize variables to be tested
vars <- c("g0_Delinq_SD_4", "g0_Delinq_SD_5", "g0_Delinq_SD_6", "g0_Delinq_SD_9", "g0_Delinq_SD_12")

# Goodness of fit test
csTable(datCredit_train_TFD,vars)
#         Variable     B
# Variable B_Statistic
# 1  g0_Delinq_SD_4      0.6186
# 5 g0_Delinq_SD_12      0.6157
# 2  g0_Delinq_SD_5      0.6101
# 4  g0_Delinq_SD_9      0.5955
# 3  g0_Delinq_SD_6      0.5901

### RESULTS:  [g0_Delinq_SD_4] fits the data the best, slightly better than [g0_Delinq_SD_5]

# Accuracy test
concTable(datCredit_train_TFD, datCredit_valid_TFD, vars)
#           Variable Concordance           SD LR_Statistic
# 1:  g0_Delinq_SD_4   0.9928422 0.0004451396        84702
# 2:  g0_Delinq_SD_5   0.9908799 0.0007682680        91976
# 3:  g0_Delinq_SD_6   0.9753519 0.0015888086        91431
# 4:  g0_Delinq_SD_9   0.9531019 0.0021947146        87688
# 5: g0_Delinq_SD_12   0.9238008 0.0026908258        81014

### RESULTS: As the SD period increase, there is a slight decrease in concordance.
### NOTE: Concordance is extremely high with low variability

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
csTable(datCredit_train_TFD,vars,seedVal = NA)
#                   Variables Iteration_1 Iteration_2 Iteration_3 Iteration_4 Iteration_5 Average
# 1:           g0_Delinq_Ave      0.6482      0.6434      0.6482      0.6493       0.646 0.64702
# 2: g0_Delinq_Any_Aggr_Prop      0.6436      0.6504      0.6475      0.6476       0.644 0.64662
# 3:                   Range      0.0046      0.0070      0.0007      0.0017       0.002 0.00320


### RESULTS: [g0-Delinq_Ave] seems to have the better goodness of fit over 5 iterations.

# Accuracy
concTable(datCredit_train_TFD, datCredit_valid_TFD,vars)
#                   Variable Concordance          SD LR_Statistic
# 1:           g0_Delinq_Ave   0.5862387 0.004229343         1872
# 2: g0_Delinq_Any_Aggr_Prop   0.5849507 0.004230491         1840

### RESULTS: [g0-Delinq_Ave] has a slightly better concordance.

### CONCLUSION: [g0-Delinq_Ave] seems to outperform [g0_Delinq_Any_Aggr_Prop] and therefore is kept in the model
###             and [g0_Delinq_Any_Aggr_Prop_Lag_3] is removed along with [g0_Delinq_Any_Aggr_Prop] due to the high correlation
###             I expect similar results.

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

### INVESTIGATE: Should an indicator version of [g0_Delinq] be included in the model.

# An indicator for when the value is greater than 0
datCredit_train_TFD[,g0_Delinq_Ind := ifelse(g0_Delinq > 0, 1, 0)]
cox <- coxph(Surv(Start,End,Default_Ind) ~ g0_Delinq_Ind, id=LoanID, datCredit_train_TFD)
summary(cox); rm(cox)
### RESULTS: coef is unstable with high variability.
datCredit_train_TFD[,g0_Delinq_Ind := NULL]

### INVESTIGATE: Should an lagged version of [g0_Delinq] be included in the model.

# A lag version for [g0_Delinq]
datCredit_train_TFD[,g0_Delinq_Lag_1 := shift(g0_Delinq,fill=0),by=LoanID]
cox <- coxph(Surv(Start,End,Default_Ind) ~ g0_Delinq_Lag_1, id=LoanID, datCredit_train_TFD)
### RESULTS: exp(coef) tends to Inf
datCredit_train_TFD[,g0_Delinq_Lag_1 := NULL]

# Arrears
cox <- coxph(Surv(Start,End,Default_Ind) ~ Arrears, id=LoanID, datCredit_train_TFD)
summary(cox); rm(cox)
### RESULTS: Obtained a stable model

### CONCLUSION: Unable to add [Delinq_0] to the model, since the various forms'
###             coef is unstable, however, [Arrears] compiles a seemingly stable
###             model, therefore it can serve as a proxy for g0_Delinq.

# ------ 1.4 How does [PerfSpell_g0_Delinq_Num] compare to [g0_Delinq_SD_4]?

# Initialize variables
vars <- c("PerfSpell_g0_Delinq_Num", "g0_Delinq_SD_4")

# Goodness of fit
csTable(datCredit_train_TFD,vars,seedVal = NA)
#                   Variables Iteration_1 Iteration_2 Iteration_3 Iteration_4 Iteration_5 Average
# 1: PerfSpell_g0_Delinq_Num      0.6448      0.6472      0.6486      0.6470      0.6431 0.64614
# 2:          g0_Delinq_SD_4      0.6175      0.6180      0.6201      0.6176      0.6192 0.61848
# 3:                   Range      0.0273      0.0292      0.0285      0.0294      0.0239 0.02766

### RESULTS: [PerfSpell_g0_Delinq_Num] seems to have a much better goodness of fit.

# Accuracy
concTable(datCredit_train_TFD, datCredit_valid_TFD,vars)
#                   Variable Concordance           SD LR_Statistic
# 1:          g0_Delinq_SD_4   0.9928422 0.0004451396        84702
# 2: PerfSpell_g0_Delinq_Num   0.9439437 0.0006710513         6307

### RESULTS: [g0_Delinq_SD_4] seems to have better concordance.

### CONCLUSION: Keep [PerfSpell_g0_Delinq_Num] in the model and remove [g0_Delinq_SD_4] since it 
###             such a better fit.

varlist <- vecChange(varlist,Remove="g0_Delinq_SD_4")



# ------ 1.4 What is the performance of current thematic variables in univariate models?

# Build thematic model based on remaining delinquency variables.
vars <- c("PerfSpell_g0_Delinq_Num","TimeInDelinqState","g0_Delinq_Ave",
          "slc_acct_arr_dir_3", "Arrears")

# Test Goodness of fit
csTable(datCredit_train_TFD,vars)
#                   Variable B_Statistic
# 3           g0_Delinq_Ave      0.6489
# 1 PerfSpell_g0_Delinq_Num      0.6475
# 5                 Arrears      0.6455
# 2       TimeInDelinqState          NA
# 4      slc_acct_arr_dir_3          NA

### RESULTS: Fits are close to on another.
### NOTE: [TimeInDelinqState] ran out of iterations and did not converge
### NOTE: [slc_acct_arr_dir_3] exp overflow due to covariates

# ------ 1.4.1 Why does [TimeInDelinqState] not converge?

cox <- coxph(Surv(Start,End,Default_Ind) ~ TimeInDelinqState, id=LoanID, datCredit_train_TFD)
### RESULTS:  Beta tends to Inf. After some inspection on the data it relates to all
###           Defaulting events starting in a new delinquency state,
###           i.e. [TimeInDelinqState] = 1. Therefore quasi-complete separation
###           seems to be present.

### INVESTIGATE: WHETHER QUASI-COMPLETE SEPERATION IS PRESENT
datCredit_train_TFD[TimeInDelinqState!=1 & Default_Ind==1, .N]
### RESULTS: 0
datCredit_train_TFD[TimeInDelinqState==1 & Default_Ind==0, .N]
### RESULTS: 167718

### CONCLUSION: Quasi-complete separation seems to be present for [g0_Delinq]=3
###             and should therefore be removed.

# Test a lagged version of TimeInDelinqState
cox <- coxph(Surv(Start,End,Default_Ind) ~ TimeInDelinqState_Lag_1, id=LoanID,
             datCredit_train_TFD)
summary(cox); rm(cox)

### CONCLUSION: Replace old variable with new variable

### RESULTS: Variable is significant with a high concordance of 0.942, thus replace
###           [TimeInDelinqState] with [TimeInDelinqState_Lag_1]

varlist <- vecChange(varlist,Remove="TimeInDelinqState",
                     Add=data.table(vars="TimeInDelinqState_Lag_1",vartypes="acc"))


# ------ 1.4.2 Why does [slc_acct_arr_dir_3] exp overflow?

describe(datCredit_train_TFD[,slc_acct_arr_dir_3])
# Value            CURING MISSING_DATA      ROLLING         SAME
# Proportion        0.015        0.126        0.022        0.837

describe(datCredit_train_TFD[Default_Ind==1,slc_acct_arr_dir_3])
# Value            CURING MISSING_DATA      ROLLING         SAME
# Proportion        0.007        0.160        0.783        0.049

### RESULTS:  Although the ROLLING level is not dominant in datCredit_train_TFD,
###           we can clearly see that it is highly predictive of a default event
###           occurring, given its high proportion if Default_Ind == 1.

# Create an indicator function
datCredit_train_TFD[, slc_acct_arr_dir_3_ROLLING_Ind := 
                      ifelse(slc_acct_arr_dir_3 == "ROLLING", 1,0)]
cox <- coxph(Surv(Start,End,Default_Ind) ~ slc_acct_arr_dir_3_ROLLING_Ind,
             datCredit_train_TFD, id=LoanID)
### RESULTS: exp overflows

# Make the indicator categorical
cox <- coxph(Surv(Start,End,Default_Ind) ~ factor(slc_acct_arr_dir_3_ROLLING_Ind),
             datCredit_train_TFD, id=LoanID)
### RESULTS: exp overflows
datCredit_train_TFD[,slc_acct_arr_dir_3_ROLLING_Ind := NULL]

# Use an indicator variable for a change in account
cox <- coxph(Surv(Start,End,Default_Ind) ~ slc_acct_arr_dir_3_Change_Ind,
             datCredit_valid_TFD, id=LoanID)
summary(cox); rm(cox)

### CONCLUSION: Replace old variable with new variable

varlist <- vecChange(varlist,Remove=c("slc_acct_arr_dir_3") ,
                     Add=data.table(vars=c("slc_acct_arr_dir_3_Change_Ind"),
                                    vartypes=c("bin")))



# ------ 1.4.3 What is the performance of current thematic variables in univariate models?

vars <- c("PerfSpell_g0_Delinq_Num","g0_Delinq_Ave", "Arrears",
          "TimeInDelinqState_Lag_1","slc_acct_arr_dir_3_Change_Ind")

# Goodness of fit
csTable(datCredit_train_TFD,vars)
#                         Variable B_Statistic
# 5 slc_acct_arr_dir_3_Change_Ind      0.6501
# 2                 g0_Delinq_Ave      0.6489
# 1       PerfSpell_g0_Delinq_Num      0.6475
# 3                       Arrears      0.6455
# 4       TimeInDelinqState_Lag_1      0.6385

### RESULTS: All variables seems to be a good fit.

# Accuracy
concTable(datCredit_train_TFD, datCredit_valid_TFD,vars)
#                         Variable Concordance           SD LR_Statistic
# 1:                       Arrears   0.9568891 0.0017938453        13502
# 2:       PerfSpell_g0_Delinq_Num   0.9439437 0.0006710513         6307
# 3:       TimeInDelinqState_Lag_1   0.8933628 0.0031888664        31583
# 4: slc_acct_arr_dir_3_Change_Ind   0.8722179 0.0014858895        28431
# 5:                 g0_Delinq_Ave   0.5862387 0.0042293435         1872

### RESULTS: [g0_Delinq_Ave] has significant less concordance than the other variables.

### CONCLUSION: Leave all variables in the model (including [g0_Delinq_Ave], since it seems to have the best fit for the data).



# ------ 1.5 What is the performance of current thematic cox ph model?

# Build cox model based on all thematic variables
coxDelinq <- coxph(Surv(Start,End,Default_Ind) ~ g0_Delinq_Ave +
                           PerfSpell_g0_Delinq_Num + slc_acct_arr_dir_3_Change_Ind +
                           TimeInDelinqState_Lag_1 + Arrears, id=LoanID,
                         data=datCredit_train_TFD)

summary(coxDelinq)

### RESULTS: All variables are significant (p-value < 0.001)

# Test goodness of fit
GoF_CoxSnell_KS(coxDelinq,datCredit_train_TFD,GraphInd = FALSE) # 0.6394

### RESULTS: Goodness of fit is a bit low

# Test accuracy
concordance(coxDelinq, newdata=datCredit_valid_TFD) # Concordance= 0.9638 se= 0.00161

### RESUTLS: extremely good prediction accuracy.

### CONCLUSION: Keep all variables in the model

# House keeping
rm(coxDelinq)

# (0,3) (4,12) (13,24) (0,12) (0,36)
#timedROC(datCredit_valid_TFD, coxDelinq_valid, month_Start=0, month_End=36,
         # fld_ID="LoanID", fld_Event="Default_Ind",fld_StartTime="Start",
         # fld_EndTime="End", numDigits=0, Graph=FALSE)
# AUC: 0.9760028
#HW: combine ggplot objects into graph facets

### RESULTS: Graph is extremely accurate with a AUC of 97.6%

# # ------ 1.6 Compare results with that of an algorithm
# 
# # My final cox model
# myCox <- coxph(Surv(Start, End, Default_Ind) ~ g0_Delinq_Ave)
# 
# # Cox model with all the other variables.
# Fullcox <- coxph(Surv(Start, End, Default_Ind) ~
#                    PerfSpell_g0_Delinq_Num + g0_Delinq_Any_Aggr_Prop +
#                    g0_Delinq_Ave + g0_Delinq_SD_4 +  g0_Delinq_SD_5 + 
#                    g0_Delinq_Any_Aggr_Prop + g0_Delinq_Any_Aggr_Prop_Lag_1 + 
#                    g0_Delinq_Any_Aggr_Prop_Lag_2 + g0_Delinq_Any_Aggr_Prop_Lag_3 + 
#                    g0_Delinq_Any_Aggr_Prop_Lag_5, datCredit_train_TFD)
# 
# compCox <- stepAIC(Fullcox, direction = "backward", trace = TRUE)
# #                                 Df    AIC
# # <none>                              98216
# # - g0_Delinq_Any_Aggr_Prop_Lag_5  1  98218
# # - g0_Delinq_Any_Aggr_Prop_Lag_1  1  98226
# # - g0_Delinq_SD_4                 1  98235
# # - PerfSpell_g0_Delinq_Num        1  98601
# # - g0_Delinq_SD_5                 1 105649
#============================================================================================

modelVar <- c("PerfSpell_g0_Delinq_Num", "Arrears", "g0_Delinq_Ave",
              "TimeInDelinqState_Lag_1", "slc_acct_arr_dir_3_Change_Ind")

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

# Compare goodness of fit of different variables
csTable(datCredit_train_TFD,vars)
#                                   Variable B_Statistic
# 1        slc_acct_pre_lim_perc_imputed_med      0.6492
# 2 slc_acct_prepaid_perc_dir_12_imputed_med          NA

### RESULTS:  [slc_acct_prepaid_perc_dir_12_imputed_med] ran out of iterations and did not converge

### INVESTIGATE: Why did slc_acct_prepaid_perc_dir_12_imputed_med] ran out of iterations and not converge?

hist(datCredit_train_TFD[,slc_acct_prepaid_perc_dir_12_imputed_med])
### RESULTS: A significant portion of the values are 0 with possible extreme values.

# Proportion of values being 0
datCredit_train_TFD[slc_acct_prepaid_perc_dir_12_imputed_med==0,.N]/datCredit_train_TFD[,.N]
# 0.734379

hist(datCredit_train_TFD[Default_Ind==1, slc_acct_prepaid_perc_dir_12_imputed_med])
### Similar distribution with most of the values being 0, but all values are below 3.5

# Proportion of values being 0 | Default_Ind == 1
datCredit_train_TFD[slc_acct_prepaid_perc_dir_12_imputed_med==0 & Default_Ind == 1,.N]/datCredit_train_TFD[Default_Ind == 1,.N]
# 0.9980394

# Proportion of values defaulted | slc_acct_prepaid_perc_dir_12_imputed_med == 0
datCredit_train_TFD[slc_acct_prepaid_perc_dir_12_imputed_med==0 & Default_Ind == 1,.N]/datCredit_train_TFD[slc_acct_prepaid_perc_dir_12_imputed_med==0,.N]
# 0.002991395

# RESULTS: (Default_Ind == 1) => (slc_acct_prepaid_perc_dir_12_imputed_med == 0)

### CONCLUSION: Quasi-complete separation seems to be present for [slc_acct_prepaid_perc_dir_12_imputed_med],
###             therefore remove it from the varlist

varlist <- vecChange(varlist,Remove="slc_acct_prepaid_perc_dir_12_imputed_med")

# ------ 2.2 What is the predictive power of the variables left in varlist?

# Initialize variables to be tested
vars <- c("pmnt_method_grp", "slc_acct_pre_lim_perc_imputed_med")

csTable(datCredit_train_TFD,vars)
#                             Variable B_Statistic
# 2 slc_acct_pre_lim_perc_imputed_med      0.6489
# 1                   pmnt_method_grp      0.6461

### RESULTS:  [slc_acct_pre_lim_perc_imputed_med] fits the better and there is an improvement for [pmnt_method_grp]

concTable(datCredit_train_TFD, datCredit_valid_TFD,vars)
#                             Variable Concordance           SD LR_Statistic
# 1:                   pmnt_method_grp   0.7643746 0.0036472739         8864
# 2: slc_acct_pre_lim_perc_imputed_med   0.6496281 0.0005780128         5558

### RESULTS: [pmnt_method_grp]  has a much better concordance than [slc_acct_pre_lim_perc_imputed_med],
###           but both concordances are reasonably high.

### Conclusion: All values should be kept in the model



# ------ 2.3 How predictive is a single model based on the thematic variables?

# Initialize variables
vars <- c("pmnt_method_grp", "slc_acct_pre_lim_perc_imputed_med")

# Build Cox model based on each variable
coxEngineered <- coxph(Surv(Start, End, Default_Ind) ~ pmnt_method_grp +
                               slc_acct_pre_lim_perc_imputed_med, id=LoanID,
                             data=datCredit_train_TFD)
summary(coxEngineered)

### RESULTS: [slc_acct_pre_lim_perc_imputed_med]  has an insignificant coef of -154

# Goodness of fit
GoF_CoxSnell_KS(coxEngineered,datCredit_train_TFD, GraphInd=FALSE) # 0.6421

### RESULTS: slight decrease in goodness of fit, although it could be due to random sampling

# Accuracy
concordance(coxEngineered,newdata=datCredit_valid_TFD) # Concordance= 0.8112 se= 0.002737

### CONCLUSION: Due to increase in Concordance keep all current variables in the model.

# House keeping
rm(coxEngineered); gc()

#============================================================================================

modelVar <- c(modelVar,"slc_acct_pre_lim_perc_imputed_med","pmnt_method_grp")

#============================================================================================
# ------ 3. Interest Rate
varlist <- data.table(vars=c("InterestRate_Nom", "InterestRate_Margin"),
                      vartypes=c("prc","prc"))

#=========================================================================================

# ------ 3.1 How correlated are the two variables?

# - Correlation analysis
corrAnalysis(datCredit_train_TFD, varlist[vartypes!="cat"]$vars, corrThresh = 0.6, method = 'spearman') # Obtain correlation groups

### RESULTS: Correlation of 0.41, therefore no significant correlations.

# Obtain the B-statistics for goodness of fit tests.
csTable(datCredit_train_TFD, varlist$vars, seedVal=NA)
#             Variables Iteration_1 Iteration_2 Iteration_3 Iteration_4 Iteration_5 Average
# 1:    InterestRate_Nom      0.6448      0.6499      0.6499      0.6503      0.6483 0.64864
# 2: InterestRate_Margin      0.6461      0.6481      0.6464      0.6474      0.6520 0.64800
# 3:               Range      0.0013      0.0018      0.0035      0.0029      0.0037 0.00264

### RESULTS:  [InterestRate_Nom] is a better fit than [InterestRate_Margin]

# Obtain the concordance for different variables.
concTable(datCredit_train_TFD, datCredit_valid_TFD, varlist$vars)
#                 Variable Concordance LR_Statistic
# 1:    InterestRate_Nom   0.5497724          165
# 2: InterestRate_Margin   0.5452677          174

### RESULTS:  [InterestRate_Nom] has a slightly higher concordance than [InterestRate_Margin].

### CONCLUSION: Variables are quite similar in predictive power, however are dissimilar
###              in correlation, therefore both can be kept in the model.



# ------ 3.2 How does [InterestRate_Margin] compare with [InterestRate_Margin_imputed_bin]?

# Initialize variables
vars <- c("InterestRate_Margin", "InterestRate_Margin_imputed_bin")

# Obtain the B-statistics for goodness of fit tests.
csTable(datCredit_train_TFD, vars, seedVal=NA)
#                           Variables Iteration_1 Iteration_2 Iteration_3 Iteration_4 Iteration_5 Average
# 1: InterestRate_Margin_imputed_bin      0.6483      0.6479      0.6487      0.6495      0.6455 0.64798
# 2:             InterestRate_Margin      0.6480      0.6461      0.6473      0.6497      0.6453 0.64728
# 3:                           Range      0.0003      0.0018      0.0014      0.0002      0.0002 0.00078
### RESULTS: [InterestRate_Margin_imputed_bin] has a better fit than [InterestRate_Margin]

# Obtain the concordance for different variables.
concTable(datCredit_train_TFD, datCredit_valid_TFD, vars)
#                           Variable Concordance          SD LR_Statistic
# 1: InterestRate_Margin_imputed_bin   0.5480786 0.003948766          454
# 2:             InterestRate_Margin   0.5452677 0.004064693          279

### RESULTS:  [InterestRate_Margin_imputed_bin] has the higher concordance.

### CONCLUSIONS: Replace [InterestRate_Margin] with [InterestRate_Margin_imputed_bin]

varlist <- vecChange(varlist, Remove="InterestRate_Margin", Add=data.table(vars="InterestRate_Margin_imputed_bin", vartypes="cat"))



# ------ 3.3 How does [InterestRate_Margin], InterestRate_Margin_Aggr_Med and [M_Repo_Rate] compare with one another?

vars <- c("InterestRate_Margin_imputed_bin", "M_Repo_Rate", "InterestRate_Margin_Aggr_Med")

# Evaluate their correlation
corrAnalysis(datCredit_train_TFD,vars)
### RESULTS: no significant correlations were detected between the variables.

# Obtain the B-statistics for goodness of fit tests.
csTable(datCredit_train_TFD, vars, seedVal=NA)
#                           Variables Iteration_1 Iteration_2 Iteration_3 Iteration_4 Iteration_5 Average
# 1:    InterestRate_Margin_Aggr_Med      0.6506      0.6491      0.6523      0.6459      0.6472 0.64902
# 2: InterestRate_Margin_imputed_bin      0.6488      0.6453      0.6488      0.6489      0.6450 0.64736
# 3:                     M_Repo_Rate      0.6449      0.6489      0.6494      0.6487      0.6438 0.64714
# 4:                           Range      0.0057      0.0038      0.0035      0.0030      0.0034 0.00388
### RESULTS:  [InterestRate_Margin_Aggr_Med] seems to be have the best fit.

# Obtain the concordance for different variables.
concTable(datCredit_train_TFD, datCredit_valid_TFD, vars)
#                           Variable Concordance          SD LR_Statistic
# 1:    InterestRate_Margin_Aggr_Med   0.5506784 0.004195174         1198
# 2: InterestRate_Margin_imputed_bin   0.5480786 0.003948766          454
# 3:                     M_Repo_Rate   0.5180882 0.004234430         1078

### RESULTS:  [InterestRate_Margin_Aggr_Med] has the highest concordance.

### CONCLUSIONS: Keep [InterestRate_Margin_Aggr_Med] in the model.

varlist <- vecChange(varlist,Remove="InterestRate_Margin_imputed_bin", Add=data.table(vars="InterestRate_Margin_Aggr_Med", vartypes="prc"))


# ------ 3.2 Should we add lagging variables to the model?

vars <- c("InterestRate_Margin_Aggr_Med","InterestRate_Margin_Aggr_Med_1",
          "InterestRate_Margin_Aggr_Med_12","InterestRate_Margin_Aggr_Med_2",
          "InterestRate_Margin_Aggr_Med_3","InterestRate_Margin_Aggr_Med_4",
          "InterestRate_Margin_Aggr_Med_5","InterestRate_Margin_Aggr_Med_6",
          "InterestRate_Margin_Aggr_Med_9")

# Obtain the B-statistics for goodness of fit tests.
csTable(datCredit_train_TFD, vars, seedVal=NA)
#                           Variables Iteration_1 Iteration_2 Iteration_3 Iteration_4 Iteration_5 Average
# 1: InterestRate_Margin_Aggr_Med_12      0.6449      0.6520      0.6510      0.6472      0.6511 0.64924
# 2:  InterestRate_Margin_Aggr_Med_3      0.6500      0.6484      0.6492      0.6491      0.6492 0.64918
# 3:    InterestRate_Margin_Aggr_Med      0.6479      0.6474      0.6485      0.6501      0.6497 0.64872
# 4:  InterestRate_Margin_Aggr_Med_9      0.6504      0.6476      0.6471      0.6482      0.6495 0.64856
# 5:  InterestRate_Margin_Aggr_Med_1      0.6477      0.6483      0.6469      0.6506      0.6480 0.64830
# 6:  InterestRate_Margin_Aggr_Med_5      0.6499      0.6455      0.6463      0.6459      0.6507 0.64766
# 7:  InterestRate_Margin_Aggr_Med_2      0.6461      0.6460      0.6453      0.6482      0.6490 0.64692
# 8:  InterestRate_Margin_Aggr_Med_6      0.6464      0.6457      0.6462      0.6475      0.6466 0.64648
# 9:  InterestRate_Margin_Aggr_Med_4      0.6481      0.6468      0.6434      0.6478      0.6448 0.64618
# 10:                           Range      0.0055      0.0065      0.0076      0.0047      0.0063 0.00612

### RESULTS:  [InterestRate_Margin_Aggr_Med_12] and [InterestRate_Margin_Aggr_Med_3] seems to be better fits than [InterestRate_Margin_Aggr_Med]

# Obtain the concordance for different variables.
concTable(datCredit_train_TFD, datCredit_valid_TFD, vars)
#                           Variable Concordance          SD LR_Statistic
# 1:  InterestRate_Margin_Aggr_Med_4   0.5611257 0.004199905          104
# 2:  InterestRate_Margin_Aggr_Med_3   0.5603674 0.004199684          105
# 3:  InterestRate_Margin_Aggr_Med_2   0.5575290 0.004200634          102
# 4:  InterestRate_Margin_Aggr_Med_5   0.5566768 0.004176251           99
# 5:  InterestRate_Margin_Aggr_Med_6   0.5566370 0.004185777           98
# 6:  InterestRate_Margin_Aggr_Med_1   0.5555107 0.004193798           99
# 7:  InterestRate_Margin_Aggr_Med_9   0.5546448 0.004170861           85
# 8:    InterestRate_Margin_Aggr_Med   0.5506784 0.004195174           96
# 9: InterestRate_Margin_Aggr_Med_12   0.5481947 0.004155018           68

### RESULTS:  [InterestRate_Margin_Aggr_Med_4] seems to have the best concordance, 
###           with [InterestRate_Margin_Aggr_Med_3] falling close behind.

### CONCLUSION: Add [InterestRate_Margin_Aggr_Med_3] to the model since it seems
###             to have the second best goodness of fit and concordance.

varlist <- vecChange(varlist,Add=data.table(vars=c("InterestRate_Margin_Aggr_Med_3"),
                                            vartypes=c("prc")))



# ------ 3.3 What is the performance of current thematic variables in univariate models?

vars <- c("InterestRate_Nom", "InterestRate_Margin_Aggr_Med", "InterestRate_Margin_Aggr_Med_3")

# Goodness of fit
csTable(datCredit_train_TFD,vars)
#                           Variable B_Statistic
# 1               InterestRate_Nom       0.6479
# 2   InterestRate_Margin_Aggr_Med       0.6477
# 3 InterestRate_Margin_Aggr_Med_3       0.6449

### RESULTS: Variables seem to be a good fit.

# Accuracy
concTable(datCredit_train_TFD, datCredit_valid_TFD,vars)
#                         Variable Concordance          SD LR_Statistic
# 1: InterestRate_Margin_Aggr_Med_3   0.5603674 0.004199684          105
# 2:   InterestRate_Margin_Aggr_Med   0.5506784 0.004195174           96
# 3:               InterestRate_Nom   0.5497724 0.004486873          165

### RESULTS: Concordances are quite low, this may lead to not including the variables into the model

### CONCLUSION: Take note of low concordances, but evaluate the concordance of the full model.



# ------ 3.4 What is the performance of current thematic cox ph model?

# Build cox model based on all thematic variables
coxInterest <- coxph(Surv(Start,End,Default_Ind) ~ InterestRate_Nom +
                           InterestRate_Margin_Aggr_Med + InterestRate_Margin_Aggr_Med_3,
                         id=LoanID, data=datCredit_train_TFD)
summary(coxInterest)

### RESULTS: All variables are significant.

# Kolmogorov-Smirnof of coxDelinq
GoF_CoxSnell_KS(coxInterest,datCredit_train_TFD,GraphInd = FALSE) # 0.646
### RESULTS: Good fit

# Accuracy
concordance(coxInterest, newdata=datCredit_valid_TFD)# Concordance= 0.5638 se= 0.00442
### RESULTS: Improved accuracy, but not particularly good.

timedROC(datCredit_valid_TFD, coxInterest_valid, month_Start=0, month_End=36,
         fld_ID="LoanID", fld_Event="Default_Ind",fld_StartTime="Start",
         fld_EndTime="End", numDigits=0, Graph=FALSE)
# AUC: 0.5214215

# House keeping
rm(coxInterest)
#===========================================================================================

modelVar <- c(modelVar,"InterestRate_Nom", "InterestRate_Margin_Aggr_Med",
              "InterestRate_Margin_Aggr_Med_3")

#============================================================================================
# ------ 4. General
varlist <- data.table(vars=c("Balance","Instalment",
                             "Principal","Undrawn_Amt",
                             "LN_TPE","Term"),
                      vartypes=c("int", "fin", "fin","int","cat", "int"))

#=========================================================================================
# ------ 4.1 Which variables are highly correlated in the group?

# - Correlation analysis
corrAnalysis(datCredit_train_TFD, varlist[vartypes!="cat"]$vars, corrThresh = 0.6, method = 'spearman') # Obtain correlation groups

### RESULTS: 1) [Balance], [Installment] and [Principal] are highly correlated.

# ------ 4.2 Which variable in the correlated group should be kept in the model?

# Initialize variable
vars <- c("Balance", "Instalment", "Principal")

# Goodness of fit test of variables
csTable(datCredit_train_TFD, vars, seedVal=NA)
#     Variables Iteration_1 Iteration_2 Iteration_3 Iteration_4 Iteration_5 Average
# 1: Instalment      0.6469      0.6492      0.6462      0.6506      0.6491 0.64840
# 2:    Balance      0.6485      0.6453      0.6463      0.6501      0.6482 0.64768
# 3:  Principal      0.6510      0.6434      0.6471      0.6472      0.6478 0.64730
# 4:      Range      0.0041      0.0058      0.0009      0.0034      0.0013 0.00310

### RESULTS:  Difficult to make a conclusion with varying B-statistics.

# Obtain the concordance for different variables.
concTable(datCredit_train_TFD, datCredit_valid_TFD, vars)
#     Variable Concordance          SD LR_Statistic
# 1:  Principal   0.5827355 0.004098595          335
# 2: Instalment   0.5636894 0.004286208           84
# 3:    Balance   0.5544691 0.004145147          101

## RESULTS: [Principal] has the highest concordance

### CONCOLUSION:  Considering Principal has the best concordance, but the worse fit
###               complicates the decision.

# ------ 4.2.1 Can [Balance] and [Instalment] be replaced with [InstalmentToBalance_Aggr_Prop]?

# Initialize variable
vars <- c("Balance", "Instalment", "InstalmentToBalance_Aggr_Prop")

# Goodness of fit test of variables
csTable(datCredit_train_TFD, vars, seedVal=NA)
#                         Variables Iteration_1 Iteration_2 Iteration_3 Iteration_4 Iteration_5 Average
# 1:                    Instalment      0.6478      0.6507      0.6512      0.6501      0.6495 0.64986
# 2: InstalmentToBalance_Aggr_Prop      0.6510      0.6501      0.6489      0.6497      0.6482 0.64958
# 3:                       Balance      0.6464      0.6495      0.6457      0.6483      0.6478 0.64754
# 4:                         Range      0.0046      0.0012      0.0055      0.0018      0.0017 0.00296

### RESULTS:  From the 5 iterations, [InstalmentToBalance_Aggr_Prop] consistently seems to outfit [Balance]

# Accuracy of variables
concTable(datCredit_train_TFD, datCredit_valid_TFD, vars)
#                         Variable Concordance          SD LR_Statistic
# 1:                    Instalment   0.5636894 0.004286208           84
# 2:                       Balance   0.5544691 0.004145147          101
# 3: InstalmentToBalance_Aggr_Prop   0.5146318 0.004341557           18

### RESULTS:  [InstalmentToBalance_Aggr_Prop] has a worse concordance than [Instalment] and [Balance].

### CONCLUSION: [InstalmentToBalance_Aggr_Prop] cannot replace [Instalment] and [Balance].

# ------ 4.2.2 Can [Balance] be replaced with [ArrearsToBalance_Aggr_Prop]?

# Initialize variable
vars <- c("Balance", "ArrearsToBalance_Aggr_Prop")

# Goodness of fit test of variables
csTable(datCredit_train_TFD, vars, seedVal=NA)
#                     Variables Iteration_1 Iteration_2 Iteration_3 Iteration_4 Iteration_5 Average
# 1:                    Balance      0.6475      0.6519      0.6479      0.6482      0.6483 0.64876
# 2: ArrearsToBalance_Aggr_Prop      0.6457      0.6440      0.6480      0.6446      0.6470 0.64586
# 3:                      Range      0.0018      0.0079      0.0001      0.0036      0.0013 0.00294

### RESULTS:  From the 5 iterations, [Balance] seems to have the best fit.

# Test for accuracy
concTable(datCredit_train_TFD, datCredit_valid_TFD, vars)
#                     Variable Concordance          SD LR_Statistic
# 1:                    Balance   0.5544691 0.004145147          101
# 2: ArrearsToBalance_Aggr_Prop   0.5367480 0.004126016          118

## RESULTS: [Balance] has the highest concordance.

### CONCOLUSION:  [ArrearsToBalance_Aggr_Prop] cannot reaplace [Balance]

### CONCLUSION: Remove [Instalment] and [Balance].

varlist <- vecChange(varlist,Remove=c("Instalment","Balance"))

# ------ 4.3 How does [Term] compare to [AgeToTerm_Aggr_Mean]?

# Initialize variables
vars <- c("Term", "AgeToTerm_Aggr_Mean")

# Goodness of fit
csTable(datCredit_train_TFD,vars,seedVal = NA)
#               Variables Iteration_1 Iteration_2 Iteration_3 Iteration_4 Iteration_5 Average
# 1: AgeToTerm_Aggr_Mean      0.6462      0.6514      0.6487      0.6471      0.6509 0.64886
# 2:                Term      0.6455      0.6458      0.6486      0.6455      0.6465 0.64638
# 3:               Range      0.0007      0.0056      0.0001      0.0016      0.0044 0.00248

### RESULTS: [AgeToTerm_Aggr_Mean] seems to consistently outfit [Term].

# Accuracy
concTable(datCredit_train_TFD, datCredit_valid_TFD,vars)
#               Variable Concordance          SD LR_Statistic
# 1: AgeToTerm_Aggr_Mean   0.5091579 0.003670705           24
# 2:                Term   0.5066915 0.002188255            2

### RESULTS: Both have poor concordance.

### CONCLUSION: Remove [Term].

varlist <- vecChange(varlist,Remove=c("Term"))

# ------ 4.4 What is the predictive power of the variables in univariate models?

vars <- c("Principal", "Undrawn_Amt", "LN_TPE")

# Goodness of fit
csTable(datCredit_train_TFD,vars, seedVal=NA)
#       Variables Iteration_1 Iteration_2 Iteration_3 Iteration_4 Iteration_5 Average
# 1: Undrawn_Amt      0.6474      0.6479      0.6505      0.6484      0.6453 0.64790
# 2:      LN_TPE      0.6466      0.6449      0.6460      0.6451      0.6481 0.64614
# 3:   Principal      0.6439      0.6465      0.6455      0.6455      0.6478 0.64584
# 4:       Range      0.0035      0.0030      0.0050      0.0033      0.0028 0.00352

### RESULTS: [Undrawn_Amt] appears to have the best goodness of fit. The other two variables'
###           vary quite a lot.

# Accuracy
concTable(datCredit_train_TFD, datCredit_valid_TFD,vars)
#       Variable Concordance           SD LR_Statistic
# 1: Undrawn_Amt   0.6883511 0.0006726495         2126
# 2:   Principal   0.5827355 0.0040985952          335
# 3:      LN_TPE   0.5121296 0.0018064070           29

### RESULTS:  [Undrawn_Amt] appears to have significant more concordance than the other variables.
###           [LN_TPE] could possibly be removed.

# ------ 4.5 What is the predictive power of the variables in a single model?

# Goodness of fit
coxGen <- coxph(as.formula(paste0("Surv(Start,End,Default_Ind) ~ ",
                                  paste(vars,collapse=" + "))), id=LoanID,
                      datCredit_train_TFD)

summary(coxGen)# All varaibles are significant (p-value < 0.001)

GoF_CoxSnell_KS(coxGen,datCredit_train_TFD,GraphInd = FALSE) # 0.6485
### RESULTS: Improved goodness of fit.

# Accuracy
concordance(coxGen, newdata=datCredit_valid_TFD) # Concordance= 0.6765 se= 0.003312
### RESULTS: Good Concordance

# House keeping
rm(coxGen)

### CONCLUSION: Leave all variables in the model.

#============================================================================================

modelVar <- c(modelVar, "Principal", "Undrawn_Amt", "LN_TPE")

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

# Goodness of fit test of variables
csTable(datCredit_train_TFD, vars, seedVal=NA)
#                 Variables Iteration_1 Iteration_2 Iteration_3 Iteration_4 Iteration_5 Average
# 1:     NewLoans_Aggr_Prop      0.6471      0.6513      0.6494      0.6486      0.6494 0.64916
# 2: CuringEvents_Aggr_Prop      0.6512      0.6474      0.6494      0.6455      0.6475 0.64820
# 3:                  Range      0.0041      0.0039      0.0000      0.0031      0.0019 0.00260

### RESULTS:  Both seem to have a good goodness of fit.

# Obtain the concordance for different variables.
concTable(datCredit_train_TFD, datCredit_valid_TFD, vars)
#                 Variable Concordance          SD LR_Statistic
# 1:     NewLoans_Aggr_Prop   0.5571220 0.004159613          194
# 2: CuringEvents_Aggr_Prop   0.5359976 0.004123115           82

## RESULTS: The concordances are not particularly good.

### CONCOLUSION:  Keep both variables in the model

# ------ 5.2 What is the performance of the variables in a single model?

# Goodness of fit
coxPort <- coxph(as.formula(paste0("Surv(Start,End,Default_Ind) ~ ",
                                        paste(vars,collapse=" + "))), id=LoanID,
                      datCredit_train_TFD)

summary(coxPort)
### RESULS: [CuringEvents_Aggr_Prop] has a high se(coef) ((22)realative to coef (42)
###         and does not have a significant p-value. [NewLoans_Aggr_Prop] is significant.

### CONCLUSION: Remove [CuringEvents_Aggr_Prop]

#============================================================================================

modelVar <- c(modelVar,"NewLoans_Aggr_Prop")

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

# Obtain the B-statistics for goodness of fit tests.
csTable(datCredit_train_TFD, vars, seedVal=NA)
#             Variables Iteration_1 Iteration_2 Iteration_3 Iteration_4 Iteration_5 Average
# 1:        M_Emp_Growth      0.6469      0.6484      0.6468      0.6498      0.6474 0.64786
# 2:    M_RealGDP_Growth      0.6473      0.6463      0.6516      0.6449      0.6489 0.64780
# 3: M_RealIncome_Growth      0.6445      0.6467      0.6448      0.6509      0.6518 0.64774
# 4:               Range      0.0028      0.0021      0.0068      0.0060      0.0044 0.00442

### RESULTS:  Non-conclusive results as rankings fluctuated too much over 5 iterations.

# Obtain the concordance for different variables.
concTable(datCredit_train_TFD, datCredit_valid_TFD, vars)
#               Variable Concordance          SD LR_Statistic
# 1: M_RealIncome_Growth   0.5341913 0.003951199           68
# 2:        M_Emp_Growth   0.5267075 0.003999905           45
# 3:    M_RealGDP_Growth   0.5211591 0.003910739           18

### RESULTS:  [M_RealIncome_Growth] has the highest concordance, although the concordances
###           are not that high.

### CONCLUSION: Keep M_RealIncome_Growth in the model.

varlist = vecChange(varlist, Remove=c("M_Emp_Growth", "M_RealGDP_Growth"))



# ------ 6.3 What are the predictive powers of the current variables in the models?

# Initialize variables
vars <- c("M_DTI_Growth", "M_Inflation_Growth", "M_RealIncome_Growth")

# Compare goodness of fit of different variables
csTable(datCredit_train_TFD,vars,seedVal = NA)
#               Variables Iteration_1 Iteration_2 Iteration_3 Iteration_4 Iteration_5 Average
# 1:  M_Inflation_Growth      0.6479      0.6488      0.6452      0.6466      0.6503 0.64776
# 2: M_RealIncome_Growth      0.6448      0.6474      0.6470      0.6473      0.6471 0.64672
# 3:        M_DTI_Growth      0.6454      0.6441      0.6466      0.6484      0.6471 0.64632
# 4:               Range      0.0031      0.0060      0.0037      0.0018      0.0054 0.00400

### RESULTS:  After 5 iterations, there is too much variability in B-statistice, therefore no conclusions can be made.

# Compare concordance of different variables
concTable(datCredit_train_TFD, datCredit_valid_TFD,vars)
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

# Obtain the B-statistics for goodness of fit tests.
csTable(datCredit_train_TFD, vars, seedVal=NA)
#           Variables Iteration_1 Iteration_2 Iteration_3 Iteration_4 Iteration_5 Average
# 1:  M_DTI_Growth_6      0.6481      0.6506      0.6489      0.6497      0.6484 0.64914
# 2:  M_DTI_Growth_1      0.6526      0.6474      0.6466      0.6499      0.6445 0.64820
# 3:  M_DTI_Growth_9      0.6466      0.6472      0.6482      0.6496      0.6464 0.64760
# 4: M_DTI_Growth_12      0.6451      0.6504      0.6480      0.6465      0.6475 0.64750
# 5:  M_DTI_Growth_3      0.6467      0.6496      0.6456      0.6461      0.6452 0.64664
# 6:    M_DTI_Growth      0.6482      0.6474      0.6485      0.6446      0.6439 0.64652
# 7:           Range      0.0075      0.0034      0.0033      0.0053      0.0045 0.00480

### RESULTS: Only [M_DTI_Growth_6], [M_DTI_Growth_1] and [M_DTI_Growth_9] appears to have the best fit.

# Obtain the concordance for different variables.
concTable(datCredit_train_TFD, datCredit_valid_TFD, vars)
#           Variable Concordance          SD LR_Statistic
# 1: M_DTI_Growth_12   0.5660056 0.004260414         1044
# 2:  M_DTI_Growth_9   0.5563671 0.004283938         1220
# 3:  M_DTI_Growth_6   0.5447379 0.004238437         1284
# 4:  M_DTI_Growth_3   0.5361471 0.004117214         1276
# 5:  M_DTI_Growth_1   0.5360005 0.004008607         1387
# 6:    M_DTI_Growth   0.5350355 0.003955393         1450

### RESULTS: Only [M_DTI_Growth_12], [M_DTI_Growth_9] and [M_DTI_Growth_6] appears to have the best fit.

### CONCLUSION: Include the following lags in the model: [M_DTI_Growth_6] and [M_DTI_Growth_9]

varlist <- vecChange(varlist,Add=data.table(vars=c("M_DTI_Growth_6", "M_DTI_Growth_9"),
                                            vartypes=c("prc", "prc")))

# ------ 6.4.2 Should lags for [M_Inflation_Growth] be included in the models?
# Initialize variables
vars <- c("M_Inflation_Growth_1","M_Inflation_Growth_12","M_Inflation_Growth_3",
          "M_Inflation_Growth_6","M_Inflation_Growth_9","M_Inflation_Growth")

# Obtain the B-statistics for goodness of fit tests.
csTable(datCredit_train_TFD, vars, seedVal=NA)
#                 Variables Iteration_1 Iteration_2 Iteration_3 Iteration_4 Iteration_5 Average
# 1:  M_Inflation_Growth_6      0.6484      0.6488      0.6489      0.6482      0.6503 0.64892
# 2:  M_Inflation_Growth_1      0.6497      0.6467      0.6494      0.6467      0.6498 0.64846
# 3:  M_Inflation_Growth_3      0.6468      0.6471      0.6485      0.6482      0.6480 0.64772
# 4: M_Inflation_Growth_12      0.6476      0.6479      0.6495      0.6452      0.6469 0.64742
# 5:    M_Inflation_Growth      0.6479      0.6478      0.6412      0.6471      0.6463 0.64606
# 6:  M_Inflation_Growth_9      0.6475      0.6496      0.6469      0.6437      0.6423 0.64600
# 7:                 Range      0.0029      0.0029      0.0083      0.0045      0.0080 0.00532

### RESULTS: [M_Inflation_Growth_6], [M_Inflation_Growth_1] and [M_Inflation_Growth_3] appear to have the best fits.

# Obtain the concordance for different variables.
concTable(datCredit_train_TFD, datCredit_valid_TFD, vars)
#                 Variable Concordance          SD LR_Statistic
# 1: M_Inflation_Growth_12   0.5582303 0.004239198          693
# 2:  M_Inflation_Growth_6   0.5529680 0.004358212         1190
# 3:  M_Inflation_Growth_9   0.5517474 0.004360946         1009
# 4:  M_Inflation_Growth_3   0.5432397 0.004326945         1155
# 5:  M_Inflation_Growth_1   0.5275290 0.004264846          976
# 6:    M_Inflation_Growth   0.5202318 0.004253188          915

### RESULTS: [M_Inflation_Growth_12], [M_Inflation_Growth_6] and [M_Inflation_Growth_9] appear to be the most accurate.

### CONCLUSION: Add [M_Inflation_Growth_6] to the model since it seems
###             to have the first best goodness of fit and second best accuracy

varlist <- vecChange(varlist,Add=data.table(vars=c("M_Inflation_Growth_6"),
                                            vartypes=c("prc")))

# ------ 6.4.3 Should lags for [M_RealIncome_Growth] be included in the models?
vars <- c("M_RealIncome_Growth_1","M_RealIncome_Growth_12","M_RealIncome_Growth_3",
          "M_RealIncome_Growth_6","M_RealIncome_Growth_9","M_RealIncome_Growth")

# Obtain the B-statistics for goodness of fit tests.
csTable(datCredit_train_TFD, vars, seedVal=NA)
#                   Variables Iteration_1 Iteration_2 Iteration_3 Iteration_4 Iteration_5 Average
# 1:  M_RealIncome_Growth_6      0.6502      0.6457      0.6486      0.6506      0.6488 0.64878
# 2:  M_RealIncome_Growth_9      0.6495      0.6436      0.6481      0.6487      0.6512 0.64822
# 3:  M_RealIncome_Growth_3      0.6484      0.6458      0.6470      0.6512      0.6470 0.64788
# 4:  M_RealIncome_Growth_1      0.6477      0.6474      0.6450      0.6438      0.6491 0.64660
# 5: M_RealIncome_Growth_12      0.6448      0.6451      0.6472      0.6494      0.6461 0.64652
# 6:    M_RealIncome_Growth      0.6445      0.6447      0.6441      0.6466      0.6477 0.64552
# 7:                  Range      0.0057      0.0038      0.0045      0.0074      0.0051 0.00530

### RESULTS: [M_RealIncome_Growth_6], [M_RealIncome_Growth_9] and [M_RealIncome_Growth_3] appear to have the best fits.

# Obtain the concordance for different variables.
concTable(datCredit_train_TFD, datCredit_valid_TFD, vars)
#     Variable Concordance          SD LR_Statistic
# 1: M_RealIncome_Growth_12   0.5131676 0.004061662          744
# 2:  M_RealIncome_Growth_9   0.4960102 0.004007988          456
# 3:  M_RealIncome_Growth_6   0.4844719 0.003940162          273
# 4:  M_RealIncome_Growth_3   0.4754604 0.003902888          141
# 5:  M_RealIncome_Growth_1   0.4681097 0.003925243           62
# 6:    M_RealIncome_Growth   0.4658087 0.003951199           35

### RESULTS: [M_RealIncome_Growth_12], [M_RealIncome_Growth_9] and [M_RealIncome_Growth_6] appear to be the most accurate,
###           althoug all of the have extremely weak concordances.

### CONCLUSION: Do not include the lags for [M_RealIncome_Growth].

# ------ 6.5 What is the predictive performance of the current univariate thematic models?

# Initialize macro-economic thematic variables
vars <- c("M_DTI_Growth", "M_Inflation_Growth","M_DTI_Growth_6", "M_DTI_Growth_9",
          "M_Inflation_Growth_6", "M_Inflation_Growth_12", "M_RealIncome_Growth")

# Test Goodness of fit
csTable(datCredit_train_TFD,vars)
#               Variable B_Statistic
# 5 M_Inflation_Growth_6      0.6500
# 6  M_RealIncome_Growth      0.6481
# 1         M_DTI_Growth      0.6479
# 2   M_Inflation_Growth      0.6477
# 4       M_DTI_Growth_9      0.6461
# 3       M_DTI_Growth_6      0.6449

### RESULTS: All variables seem to be relatively good fits with the lowest being 0.6449.

# Test Accuracy
concTable(datCredit_train_TFD, datCredit_valid_TFD, vars)
#                 Variable Concordance          SD LR_Statistic
# 1: M_Inflation_Growth_12   0.5582303 0.004239198          693
# 2:        M_DTI_Growth_9   0.5563671 0.004283938         1220
# 3:  M_Inflation_Growth_6   0.5529680 0.004358212         1190
# 4:        M_DTI_Growth_6   0.5447379 0.004238437         1284
# 5:          M_DTI_Growth   0.5350355 0.003955393         1450
# 6:    M_Inflation_Growth   0.5202318 0.004253188          915
# 7:   M_RealIncome_Growth   0.4658087 0.003951199           35

### RESULTS: All variables should be kept in the model. Although [M_RealIncome_Growth]
###           has a low concordance, it is kept in the model as a proxy for M_Emp_Growth.

# ------ 6.6 What is the predictive performance of the current thematic model?
# Initialize macro-economic thematic variables
vars <- c("M_DTI_Growth", "M_Inflation_Growth","M_DTI_Growth_6", "M_DTI_Growth_9",
          "M_Inflation_Growth_6", "M_RealIncome_Growth", "M_Inflation_Growth_12")

# Build model based on macro economic variables
coxMacro <- coxph(as.formula(paste0("Surv(Start,End,Default_Ind) ~ ",
                                    paste(vars,collapse=" + "))), id=LoanID,
                  datCredit_train_TFD)
summary(coxMacro)

### RESULTS: [M_DTI_Growth_6] and [M_DTI_Growth_9] have high p-values and seemingly 
###           the lowest Goodness of Fit (also opposite signs), therefore they should be removed.
###           [M_Inflation_Growth_12] has a high p-values and an
###           opposite signs to M_Inflation_Growth and M_Inflation_Growth_6, therefore they should be removed.

vars <- c("M_DTI_Growth", "M_Inflation_Growth", "M_Inflation_Growth_6", "M_RealIncome_Growth",
          "M_Inflation_Growth_12")

# Build model based on macro economic variables
coxMacro <- coxph(as.formula(paste0("Surv(Start,End,Default_Ind) ~ ",
                                    paste(vars,collapse=" + "))), id=LoanID,
                  datCredit_train_TFD)
summary(coxMacro)

### RESULTS:  All values have significant p-value.

# Test Goodness of fit of the single model based on all variables.
GoF_CoxSnell_KS(coxMacro,datCredit_train_TFD,GraphInd = F) # 0.6499
# Goodness of fit is equal to greatest Goodness of fit univariate model.

# Test accuracy of the single model based on all variables.
concordance(coxMacro, newdata=datCredit_valid_TFD)# 0.5509
# Concordance is much better than that of individual models.

#============================================================================================

modelVar <- c(modelVar,"M_DTI_Growth", "M_Inflation_Growth", "M_Inflation_Growth_6",
              "M_RealIncome_Growth")

#============================================================================================
# ------ 7. Final Model

# - Confirm prepared datasets are loaded into memory
if (!exists('datCredit_train_TFD')) unpack.ffdf(paste0(genPath,"creditdata_train_TFD"), tempPath);gc()
if (!exists('datCredit_valid_TFD')) unpack.ffdf(paste0(genPath,"creditdata_valid_TFD"), tempPath);gc()
### AB: If the changes are made to the process flow (i.e., moving data prep to where they belong), then the above
# alone should work before fitting the final model. Naturally, I am not going to run through the above, I am only
# interested in the final model, and will trust that you have now learned the essentials of the 
# "thematic variable selection process". ;-) Delete this comment (as with all my "### AB"-comments) after addressing/reading

# ------ 7.1 How does the univariate models compare with one another?
# Initialize variables
vars <- c("PerfSpell_g0_Delinq_Num", "Arrears", "g0_Delinq_Ave", "TimeInDelinqState_Lag_1",      
          "slc_acct_arr_dir_3_Change_Ind", "slc_acct_pre_lim_perc_imputed_med", 
          "pmnt_method_grp", "InterestRate_Margin_Aggr_Med", "InterestRate_Margin_Aggr_Med_3", 
          "Principal", "Undrawn_Amt", "LN_TPE", "NewLoans_Aggr_Prop", "M_DTI_Growth",
          "M_Inflation_Growth", "M_Inflation_Growth_6", "M_RealIncome_Growth")

# Test Goodness of fit
csTable(datCredit_train_TFD,vars)
#                              Variable B_Statistic
# 12                            LN_TPE      0.6520
# 9     InterestRate_Margin_Aggr_Med_3      0.6504
# 14                      M_DTI_Growth      0.6484
# 11                       Undrawn_Amt      0.6483
# 6  slc_acct_pre_lim_perc_imputed_med      0.6481
# 10                         Principal      0.6474
# 7                    pmnt_method_grp      0.6473
# 2                            Arrears      0.6472
# 15                M_Inflation_Growth      0.6468
# 8       InterestRate_Margin_Aggr_Med      0.6464
# 13                NewLoans_Aggr_Prop      0.6460
# 1            PerfSpell_g0_Delinq_Num      0.6458
# 17               M_RealIncome_Growth      0.6457
# 16              M_Inflation_Growth_6      0.6455
# 3                      g0_Delinq_Ave      0.6449
# 5      slc_acct_arr_dir_3_Change_Ind      0.6446
# 4            TimeInDelinqState_Lag_1      0.6331

### RESULTS: All variables seem to have comparible Goodness of fits, although
###           [TimeInDelinqState_Lag_1] have a particularly bad one.

# Test accuracy
concTable(datCredit_train_TFD, datCredit_valid_TFD, vars)
#                             Variable Concordance           SD LR_Statistic
# 1:                           Arrears   0.9573776 0.0016578273        13229
# 2:           PerfSpell_g0_Delinq_Num   0.9474644 0.0006456311         7954
# 3:           TimeInDelinqState_Lag_1   0.9064346 0.0031587395        35701
# 4:     slc_acct_arr_dir_3_Change_Ind   0.8433911 0.0019345216        33655
# 5:                   pmnt_method_grp   0.7363669 0.0036080644         8299
# 6:                       Undrawn_Amt   0.6883511 0.0006726495         2941
# 7: slc_acct_pre_lim_perc_imputed_med   0.6360895 0.0008227257         5572
# 8:                         Principal   0.5827355 0.0040985952          516
# 9:    InterestRate_Margin_Aggr_Med_3   0.5603674 0.0041996844         1174
# 10:                NewLoans_Aggr_Prop   0.5571220 0.0041596133          107
# 11:              M_Inflation_Growth_6   0.5529680 0.0043582119         1190
# 12:      InterestRate_Margin_Aggr_Med   0.5506784 0.0041951744         1198
# 13:                     g0_Delinq_Ave   0.5397435 0.0040305498         1623
# 14:                      M_DTI_Growth   0.5350355 0.0039553926         1450
# 15:                M_Inflation_Growth   0.5202318 0.0042531877          915
# 16:                            LN_TPE   0.5121296 0.0018064070          154
# 17:               M_RealIncome_Growth   0.4658087 0.0039511993           35

# ------ 7.2 What is the predictive performance of the final model containing the selected variables?
# Initialize variables
vars <- c("PerfSpell_g0_Delinq_Num", "Arrears", "g0_Delinq_Ave", "TimeInDelinqState_Lag_1",      
          "slc_acct_arr_dir_3_Change_Ind", "slc_acct_pre_lim_perc_imputed_med", 
          "pmnt_method_grp", "InterestRate_Margin_Aggr_Med", "InterestRate_Margin_Aggr_Med_3", 
          "Principal", "LN_TPE", "M_DTI_Growth","M_Inflation_Growth",
          "M_Inflation_Growth_6", "M_RealIncome_Growth")
### AB: I removed [Undrawn_Amt] since it is known to be a faulty variable. This removal was already actioned 
# in script 2f. Kindly remove any reliance in the themes above on this variable

# - Initialize variables | AB-variant
vars2 <- c("PerfSpell_g0_Delinq_Num", "Arrears", "g0_Delinq_Ave", "TimeInDelinqState_Lag_1",      
           "slc_acct_arr_dir_3_Change_Ind", "slc_acct_pre_lim_perc_imputed_med", 
           "InterestRate_Margin_Aggr_Med", "InterestRate_Margin_Aggr_Med_3", 
           "Principal", "LN_TPE", "M_DTI_Growth","M_Inflation_Growth",
           "M_Inflation_Growth_6", "M_RealIncome_Growth")
 ### AB: I received "exp overflow due to covariates" error when trying to fit this model.
# This is likely because the exponentiated form of certain variables and the 
# coxPH-function then struggles to obtain a good fit. So I simply tweaked the variable list a bit just 
# so that I can at least get some fitted model to test. Problem 'children': [pmnt_method_grp]

# Test Goodness of fit
csTable_TFD <- csTable(datCredit_train_TFD,vars2)
# 12                M_Inflation_Growth      0.6528
# 9                          Principal      0.6525
# 14               M_RealIncome_Growth      0.6508
# 7       InterestRate_Margin_Aggr_Med      0.6502
# 5      slc_acct_arr_dir_3_Change_Ind      0.6501
# 6  slc_acct_pre_lim_perc_imputed_med      0.6499
# 11                      M_DTI_Growth      0.6496
# 10                            LN_TPE      0.6489
# 2                            Arrears      0.6481
# 8     InterestRate_Margin_Aggr_Med_3      0.6477
# 1            PerfSpell_g0_Delinq_Num      0.6475
# 13              M_Inflation_Growth_6      0.6467
# 3                      g0_Delinq_Ave      0.6463
# 4            TimeInDelinqState_Lag_1      0.6385

# Test accuracy
concTable_TFD <- concTable(datCredit_train_TFD, datCredit_valid_TFD, vars2)
#                             Variable Concordance           SD LR_Statistic
# 1:                           Arrears   0.9568891 0.0017938453        13502
# 2:           PerfSpell_g0_Delinq_Num   0.9439437 0.0006710513         6307
# 3:           TimeInDelinqState_Lag_1   0.8933628 0.0031888664        31583
# 4:     slc_acct_arr_dir_3_Change_Ind   0.8722179 0.0014858895        28431
# 5: slc_acct_pre_lim_perc_imputed_med   0.6496281 0.0005780128         5558
# 6:                         Principal   0.6137518 0.0043225132          494
# 7:    InterestRate_Margin_Aggr_Med_3   0.5869687 0.0043844515         1148
# 8:                     g0_Delinq_Ave   0.5862387 0.0042293435         1872
# 9:      InterestRate_Margin_Aggr_Med   0.5858681 0.0044509364         1214
# 10:                      M_DTI_Growth   0.5813828 0.0042020200         1696
# 11:              M_Inflation_Growth_6   0.5743351 0.0045314545         1353
# 12:                M_Inflation_Growth   0.5552019 0.0045393073         1098
# 13:                            LN_TPE   0.5209834 0.0019699539          100
# 14:               M_RealIncome_Growth   0.4892308 0.0041952801           26

# Create table object
Table_TFD <- left_join(csTable_TFD$Table,concTable_TFD,by="Variable")

# 
# Build model based on variables
cox_TFD <- coxph(as.formula(paste0("Surv(Start,End,Default_Ind) ~ ",
                               paste(vars2,collapse=" + "))), id=LoanID, datCredit_train_TFD)
c <- coefficients(cox_TFD)
(c <- data.table(Variable=names(c),Coefficient=c))
#                             Variable        Coefficient
# 1:           PerfSpell_g0_Delinq_Num    0.0033543172679
# 2:                           Arrears    0.0000705150839
# 3:                     g0_Delinq_Ave   -5.5525678389239
# 4:           TimeInDelinqState_Lag_1   -0.1139113572196
# 5:     slc_acct_arr_dir_3_Change_Ind    3.3294073869448
# 6: slc_acct_pre_lim_perc_imputed_med  -33.3990668446264
# 7:       pmnt_method_grpMISSING_DATA   -1.4980174597352
# 8:    pmnt_method_grpSalary/Suspense    0.9390684205780
# 9:          pmnt_method_grpStatement    0.3823914797478
# 10:      InterestRate_Margin_Aggr_Med  445.5030569360081
# 11:    InterestRate_Margin_Aggr_Med_3 -356.2167460966670
# 12:                         Principal   -0.0000014990546
# 13:                       Undrawn_Amt    0.0000002719867
# 14:                         LN_TPEWHL    0.3001587384029
# 15:                      M_DTI_Growth    4.1911996761251
# 16:                M_Inflation_Growth    4.4736472467573
# 17:              M_Inflation_Growth_6    1.4355263441157
# 18:               M_RealIncome_Growth   -2.3490754410016

# Goodnes of fit
GoF_CoxSnell_KS(cox_TFD,datCredit_train_TFD, GraphInd=TRUE, legPos=c(0.6,0.4)) # 0.6307

### RESULTS: Goodness of fit for the model seems to be a bit low.

# Save objects
pack.ffdf(paste0(genObjPath,"TFD_Univariate_Models"), Table_TFD)
pack.ffdf(paste0(genObjPath,"TFD_Cox_Model"), cox_TFD)
