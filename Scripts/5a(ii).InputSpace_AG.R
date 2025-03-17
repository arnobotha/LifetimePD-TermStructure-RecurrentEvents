# ============================================ INPUT SPACE: AG =========================================
# Divide data into thematic groups and perform data analysis on them to compile an input space for the AG model
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
#   - datCredit_train_AG | Prepared from script 3c
#   - datCredit_valid_AG | Prepared from script 3c

#
# -- Outputs:
#   - Input_Space
# ------------------------------------------------------------------------------------------------------

# ------ 1. Preliminaries
# - Confirm prepared datasets are loaded into memory
if (!exists('datCredit_train_AG')) unpack.ffdf(paste0(genPath,"creditdata_train_AG"), tempPath);gc()
if (!exists('datCredit_valid_AG')) unpack.ffdf(paste0(genPath,"creditdata_final_AG"), tempPath);gc()


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
csTable(datCredit_train_AG,vars)
#         Variable     B
# Variable B_Statistic
# 1  g0_Delinq_SD_4      0.6186
# 5 g0_Delinq_SD_12      0.6157
# 2  g0_Delinq_SD_5      0.6101
# 4  g0_Delinq_SD_9      0.5955
# 3  g0_Delinq_SD_6      0.5901

### RESULTS:  [g0_Delinq_SD_4] fits the data the best, slightly better than [g0_Delinq_SD_5]

# Accuracy test
concTable(datCredit_train_AG, datCredit_valid_AG, vars)
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
corrAnalysis(datCredit_train_AG, varlist[vartypes!="cat"]$vars, corrThresh = 0.6, method = 'spearman') # Obtain correlation groups

### RESULTS:  1) [g0_Delinq_Any_Aggr_Prop] and [g0_Delinq_Ave] with a correlation of 1
###           2) [g0_Delinq] and [Arrears]

### CONCLUSION: A single variable from each group must be retained while the rest are removed.







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
