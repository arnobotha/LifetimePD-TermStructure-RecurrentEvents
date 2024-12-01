# ================================ DATA FUSION - TFD ====================================
# Subsample the full TFD credit dataset before fusing it to the wider "input space" data.
# Finally, analyse and fix basic deficiencies within the input space, followed by
# some basic feature engineering towards cultivating a comprehensive input space, suitable
# for Cox regression modelling.
# ---------------------------------------------------------------------------------------
# PROJECT TITLE: Default survival modelling
# SCRIPT AUTHOR(S): Dr Arno Botha, Marcel Muller, Roland Breedt, Bernard Scheepers
# ---------------------------------------------------------------------------------------
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

# -- Inputs:
#   - datInput.raw | raw input space imported in script 1
#   - datCredit_TFD | prepared credit data from script 3b
#   - various parameters set in the setup script 0
#   - datMV | prepared feature engineered macroeconomic data from script 2a
#
# -- Outputs:
#   - datCredit_final_c_TFD | enriched credit dataset, fused with macroeconomic data and
#                        various input fields, with feature engineering applied
#   - datCredit_smp_TFD   | subsampled dataset
# ---------------------------------------------------------------------------------------
# NOTE: This script predominantly comes from another project (Kasmeer).
# =======================================================================================





# ------ 1. Preliminaries
# --- 1.1 Load in Dataset
if (!exists('datCredit_TFD')) unpack.ffdf(paste0(genPath,"creditdata_final_TFD"), tempPath)




# ------- 2. Subsampling


# --- 2.1 Preliminaries and parameter definitions

# - Field Names
# Required field names
resolPerf_targetVar <- "Defaulted" # Reference level in the performance spell resolution type (as specified by [resolPerf_start]) for the target variable
clusVar <- "LoanID"
clusVar_Perf <- "PerfSpell_Key"
timeVar <- "Date"
counter <- "Counter"
PerfSpell_counter <- "PerfSpell_Counter"
smp_size <- 90000 # fixed size of downsampled set in terms of the number of unique performance spells
cat(smp_size, " is ", sprintf("%.4f", smp_size/length(unique(datCredit_TFD[,get(clusVar)]))*100), "% of all performance spells.")
smp_frac <- 0.7 # sampling fraction for resampling scheme
minStrata_size_Perf <- 20 # Minimum size for performance spells
# Optional field named (stratification)
stratifiers <- c("Date_First", "Event_Type") # First variable should be of type "date" | Assign "NA" for no stratifiers | Other good stratifier candidates are "Event_Type", "LN_TPE", and HasDefaulted_Ever
# Facet specification field names (for graphing purposes of the resolution rates)
resolPerf <- "PerfSpellResol_Type_Hist2" # Field name of performance spell resolution types - first level should be the target event (default)
resolPerf_stop <- "PerfSpellResol_Type_Hist3" # Field name for performance spell resolution rate; specific for the stopping time cohort | Assign [resolPerf] if not interested in controlling the resolution rate facet for stop dates
resolPerf_stop2 <- "Defaulted" # Name of the main resolution type (typically default) used to obtain a single facet (this level needs to be within the resolPerf_stop variable) | Set to NA if not interested in creating additional facets for the performance spells using stopping time

# - Implied sampling fraction for downsampling step
smp_perc <- smp_size/length(unique(datCredit_TFD[,get(clusVar)]))

# - Minimum strata requirements
minStrata_size <- 0 # Minimum strata size specified for subsample


# --- 2.2 Subsampling (simple clustered )
# - Set seed
set.seed(1, kind="Mersenne-Twister")
# - Conditional loop for stratifiers
if (all(is.na(stratifiers))){ # No stratifiers
  # Get unique loan account IDs from the full dataset
  dat_keys <- unique(datCredit_TFD[, mget(c(clusVar))])
  # Use simple random sampling to select the loan IDs that ought to be in the subsampled dataset
  dat_smp_keys <- dat_keys %>% slice_sample(prop=smp_perc) %>% as.data.table()
} else { # Stratifiers
  # Get unique loan account IDs from the full dataset
  dat_keys <- unique(datCredit_TFD[, mget(c(clusVar, stratifiers))])
  # Use simple random sampling with the stratifiers to select the loan IDs that ought to be in the subsampled dataset
  dat_smp_keys <- dat_keys %>% group_by(across(all_of(stratifiers))) %>% slice_sample(prop=smp_perc) %>% as.data.table()
}

# - Obtain the associated loan records as to create the subsampled dataset
datCredit_smp <- copy(datCredit_TFD[get(clusVar) %in% dat_smp_keys[, get(clusVar),]])

# House keeping
rm(datCredit_TFD)


# --- 2.3 Minimum stratum analysis and subsequent exclusions to ensure adherence to specified threshold
# - Obtaining the stratum that are below the minimum
if (all(!is.na(stratifiers))){ # - Conditional loop for strata
  selectionVar_smp <- c(clusVar, timeVar, stratifiers)
  datStrata_smp_min <- datCredit_smp[get(counter)==1, mget(selectionVar_smp)][, list(Freq = .N), by=stratifiers][Freq<minStrata_size,]
  cat(sum(datStrata_smp_min[,Freq]), "accounts of ", datCredit_smp[get(counter)==1,.N], "(", sprintf("%.4f", sum(datStrata_smp_min[,Freq])/datCredit_smp[get(counter)==1,.N]*100), "%) need to be excluded to ensure a minimum strata size of ", minStrata_size)
  
  # - Conditionally applying the exclusions
  if (sum(datStrata_smp_min[,Freq]) > 0){
    # Saving the number of records and the prior probability, in the subsampled dataset, for reporting
    datCredit_smp_old_n <- datCredit_smp[,.N]
    datCredit_smp_prior <- datCredit_smp[get(timeVar)==timeVar_Perf_Min, get(resolPerf)] %>% table() %>% prop.table(); datCredit_smp_prior <- datCredit_smp_prior[names(datCredit_smp_prior)[names(datCredit_smp_prior) == resolPerf_targetVar]][[1]] # Computing the prior probabilities of the performance spell resolution outcomes
    # Initiating a vector which will contain the exclusion IDs
    dat_keys_exc <- NA
    # Looping through the minimum strata dataset and building an exclusion condition (filter) for each row therein
    for (i in 1:datStrata_smp_min[,.N]){
      class_type <- sapply(datStrata_smp_min[,1:length(stratifiers)], function(x) {class(x[[1]])}) # Getting the type of class of each stratifier (used for building the ith condition)
      
      excCond <- datStrata_smp_min[i,1:length(stratifiers)] # Getting the values of the ith minimum strata
      excCond <- data.table(Stratifier = colnames(excCond), # Building a dataset
                            Value = unname(t(excCond)), # Ensure that the column name is Value instead of Value.V1
                            Class = class_type)
      excCond[, Value.V1 := ifelse(Class %in% c("numeric", "Date"), paste0("as.",Class,"(",'"',Value.V1,'"',")"), paste0('"', Value.V1, '"'))]
      excCond[, Condition := paste0(Stratifier, " == ", Value.V1, " & ")] # Adding an "and" operator to enable multiple conditions
      excCond2 <- parse(text = paste0(paste0(excCond$Condition, collapse = ""), counter,"==1")) # Compiling the ith condition
      
      dat_keys_exc <- c(dat_keys_exc, as.vector(datCredit_smp[eval(excCond2), get(clusVar)]))
    }
    dat_keys_exc <- dat_keys_exc[-1] # Removing the first value (as it is a missing value stemming from the vector's creation)
    
    # Applying the exclusions to the subsampled dataset
    datCredit_smp <- copy(datCredit_smp[!(get(clusVar) %in% dat_keys_exc),])
    
    cat(datCredit_smp_old_n-datCredit_smp[,.N], " observations removed (", sprintf("%.4f", (datCredit_smp_old_n-datCredit_smp[,.N])/datCredit_smp_old_n*100), "% ) \n",
        "Prior probability = ", sprintf("%.4f", datCredit_smp_prior*100), "% comapred to ", sprintf("%.4f", (datCredit_smp[get(timeVar)==timeVar_Perf_Min, get(resolPerf)] %>% table() %>% prop.table())[[2]]*100), "%")
  }
  # - Obtaining the stratum that are below the minimum
  datStrata_smp_min <- datCredit_smp[get(counter)==1, mget(selectionVar_smp)][, list(Freq = .N), by=stratifiers][Freq<minStrata_size,]
  cat(sum(datStrata_smp_min[,Freq]), "accounts of ", datCredit_smp[get(counter)==1,.N], "(", sprintf("%.4f", sum(datStrata_smp_min[,Freq])/datCredit_smp[get(counter)==1,.N]*100), "%) need to be excluded to ensure a minimum strata size of ", minStrata_size)
}




# ------- 3. Fusing credit dataset with additional input fields
# - Confirm if the input space data is loaded into memory
if (!exists('datInput.raw')) unpack.ffdf(paste0(genPath,"creditdata_input1"), tempPath)

# [SANITY CHECK] Prevalence of overlapping fields in the input space and the main credit dataset
# Find intersection between fields in input space and those perhaps already in the main credit dataset
overlap_flds <- intersect(colnames(datCredit_smp), colnames(datInput.raw))
check.fuse1 <- length(overlap_flds) == 0 # FALSE; duplicate columns exists.
cat(check.fuse1 %?% 'SAFE: No overlapping fields in the input space and the main credit dataset' %:%
      'WARNING: Overlapping field(s) detected in the input space and the main credit dataset.')
# Conditional reporting
if (check.fuse1 == 0) {cat('NOTE: The following fields overlap: ', overlap_flds,"\n",sep="\t")}
### RESULTS: slc_past_due_amt overlap

# - Remove any additional variables that are not going to be used
suppressWarnings( datInput.raw[, `:=`(slc_status_final_pred7 = NULL, slc_status_final = NULL, 
                                      slc_curing_ind = NULL, datex = NULL)])
# - Ensure variables are not present in dataset before fusion (useful during debugging)
suppressWarnings( datCredit_smp[, `:=`(slc_pmnt_method = NULL, slc_past_due_amt = NULL, slc_days_excess = NULL,
                                       slc_status_final_pred7 = NULL, slc_status_final = NULL, slc_curing_ind = NULL,
                                       slc_acct_pre_lim_perc = NULL, slc_acct_roll_ever_24 = NULL,
                                       slc_acct_arr_dir_3 = NULL, slc_acct_prepaid_perc_dir_12 = NULL, 
                                       ccm_ute_lvl_40_cnt_24m = NULL, ccm_worst_arrears_6m = NULL, ccm_worst_arrears_24m = NULL)])

# - Format the date in the correct format for merging
datInput.raw[, date := as.Date(date, format="%Y-%m-%d")]
# - Rename the datasets for merging
colnames(datInput.raw)[colnames(datInput.raw) %in% c("date", "acct_no")] <- c("Date", "LoanID")
# - Check the data grain
data_grain_check <- datInput.raw[, list(Freq = .N), by=list(LoanID, Date)][Freq>1,]
sum(is.na(data_grain_check$LoanID)); gc()
# the data grain is broken in the cases where a Loan_ID does not exist - we are not interested in these accounts in any case
# - Merge on LoanID and Date by performing a left-join
datCredit_smp <- merge(datCredit_smp, datInput.raw, by=c("Date", "LoanID"), all.x=T); gc()
# - Check the data grain
NROW(data_grain_check_merge <- datCredit_smp[, list(Freq = .N), by=list(LoanID, Date)][Freq>1,])==0
# success, the data grain check is passed

# - Save intermediary snapshot to disk (zip) for quick disk-based retrieval later
pack.ffdf(paste0(genPath,"creditdata_final_TFD_smp1a"), datCredit_smp)
# - Clean-up
rm(datInput.raw, data_grain_check, data_grain_check_merge); gc()




# ------- 4. Feature engineering for modelling purposes
if (!exists('datCredit_smp')) unpack.ffdf(paste0(genPath,"creditdata_final_TFD_smp1a"), tempPath)


# --- 4.1 Missing value diagnostics & treatments
# - Diagnostics of missing values in the additional engineered "SLC" input space | If missingness > 50% missing remove variable
# Categorical variables
table(is.na(datCredit_smp$slc_pmnt_method)) %>% prop.table()              # missingness: 11.84% - keep variable
table(is.na(datCredit_smp$slc_acct_arr_dir_3)) %>% prop.table()           # missingness: 11.84%  - keep variable
# Numerical variables
table(is.na(datCredit_smp$slc_past_due_amt)) %>% prop.table()             # missingness: 11.84% - keep variable
table(is.na(datCredit_smp$slc_days_excess)) %>% prop.table()              # missingness: 74.83% - discard variable
table(is.na(datCredit_smp$slc_acct_pre_lim_perc)) %>% prop.table()        # missingness: 11.84% - keep variable
table(is.na(datCredit_smp$slc_acct_prepaid_perc_dir_12)) %>% prop.table() # missingness: 11.84% - keep variable

# - Categorical variables
# [slc_pmnt_method] - Missing value indicators
describe(datCredit_smp$slc_pmnt_method)
### RESULTS: [slc_pmnt_method] has 7 levels and 5555434 observations and 7 missing values. Bin the missing values with the already existent "unknown" bin.
# [TREATMENT] Binning "Unknown" values and missing values into one level
datCredit_smp[, slc_pmnt_method := 
                ifelse(is.na(slc_pmnt_method) | slc_pmnt_method == "" | slc_pmnt_method == "Unknown",
                       "MISSING_DATA", slc_pmnt_method)]
describe(datCredit_smp$slc_pmnt_method)
### RESULTS: [slc_pmnt_method] has 7 levels and 0 observations have missing values.
# [TREATMENT] Apply factor transformation
datCredit_smp[,slc_pmnt_method:=factor(slc_pmnt_method)]
### RESULTS: Missing values imputed and factorization applied to the levels of the variable.

# [slc_acct_arr_dir_3] - Missing value indicators
describe(datCredit_smp$slc_acct_arr_dir_3)
### RESULTS: [slc_acct_arr_dir_3] has 4 levels and 5555434 observations have missing values.
### [TREATMENT] Binning "N/A" values and missing values into one level
datCredit_smp[, slc_acct_arr_dir_3 := 
                ifelse(is.na(slc_acct_arr_dir_3) | slc_acct_arr_dir_3 == "" | slc_acct_arr_dir_3 == "N/A",
                       "MISSING_DATA", slc_acct_arr_dir_3)]
describe(datCredit_smp$slc_acct_arr_dir_3)
### RESULTS: [slc_acct_arr_dir_3] has 4 levels and 0 missing values after binning the missing values with the already existent "N/A" bin.
# [TREATMENT] Apply factor transformation
datCredit_smp[,slc_acct_arr_dir_3:=factor(slc_acct_arr_dir_3)]
### RESULTS: Missing values imputed and facorisation applied to the levels of the variable.

# - Numerical variables
# [slc_past_due_amt] - Missing value indicators
datCredit_smp[, value_ind_slc_past_due_amt := ifelse(is.na(slc_past_due_amt) | slc_past_due_amt == "", 0, 1)]
describe(datCredit_smp$value_ind_slc_past_due_amt)
### RESULTS: No missing values, variable created successfully.

# [slc_past_due_amt] - Missing value imputation
describe(datCredit_smp$slc_past_due_amt)
### RESULTS:    [slc_past_due_amt] has scale [0;279141] and 746130 observations have missing values; with 0 at 50% quantile and 0 at 75% quantile and mean of 322.9; note a VERY large single outlier in the right hand-side tail
### CONCLUSION: Median imputation proceeds from the distribution of values which constitutes performance spells only (all associated values from default spells are excluded)
# [TREATMENT] Median imputation for missing values (due to non-normality, i.e. high skewness)
datCredit_smp[, slc_past_due_amt_imputed_med := 
                ifelse(is.na(slc_past_due_amt) | slc_past_due_amt == "", 
                       median(datCredit_smp[!is.na(PerfSpell_Key) & !(is.na(slc_past_due_amt) | slc_past_due_amt == ""), slc_past_due_amt], na.rm=TRUE), slc_past_due_amt)]
describe(datCredit_smp$slc_past_due_amt_imputed_med)
### RESULTS: [slc_past_due_amt_imputed_med] has scale [0;279141] and 0 observations have missing values; with 0 at 50% quantile and 0 at 75% quantile and mean of 284.7
hist(datCredit_smp$slc_past_due_amt_imputed_med, breaks=500); skewness(datCredit_smp$slc_past_due_amt_imputed, na.rm = T); datCredit_smp[slc_past_due_amt_imputed_med==0,.N]/datCredit_smp[,.N]
### RESULTS:    [slc_past_due_amt_imputed_med] is skewed to the right; Skewness = 18.58226; 95.12% of variables have zero values.

# [slc_acct_pre_lim_perc] - Missing value indicator
datCredit_smp[, value_ind_slc_acct_pre_lim_perc := ifelse(is.na(slc_acct_pre_lim_perc) | slc_acct_pre_lim_perc == "", 0, 1)]
describe(datCredit_smp$value_ind_slc_acct_pre_lim_perc)
### RESTULS: Binning successful, no missing values in [value_ind_slc_acct_pre_lim_perc]

# [slc_acct_pre_lim_perc] - Missing value imputation
describe(datCredit_smp$slc_acct_pre_lim_perc)
### RESULTS:    [slc_acct_pre_lim_perc] has scale [0;1] and 529862 observations with 746130 missing values; with 0 at 50% quantile and 0 at 75% quantile and mean of 0.09851
### CONCLUSION: Median imputation proceeds from the distribution of values which constitutes performance spells only (all associated values from default spells are excluded)
# [TREATMENT] Median imputation for missing values (due to non-normality, i.e. high skewness)
datCredit_smp[, slc_acct_pre_lim_perc_imputed_med := 
                ifelse(is.na(slc_acct_pre_lim_perc) | slc_acct_pre_lim_perc == "", 
                       median(datCredit_smp[!is.na(PerfSpell_Key) & !(is.na(slc_acct_pre_lim_perc) | slc_acct_pre_lim_perc == ""), slc_acct_pre_lim_perc], na.rm=TRUE), slc_acct_pre_lim_perc)]
describe(datCredit_smp$slc_acct_pre_lim_perc_imputed_med)
### RESTULS: [slc_acct_pre_lim_perc_imputed_med] has scale [0;1] and 0 observations have missing values; with 0 at 50% quantile and 0.009435 at 75% quantile and mean of 0.08684
hist(datCredit_smp$slc_acct_pre_lim_perc_imputed_med, breaks=500); skewness(datCredit_smp$slc_acct_pre_lim_perc_imputed_med, na.rm = T); datCredit_smp[slc_acct_pre_lim_perc_imputed_med==0,.N]/datCredit_smp[,.N]
### RESULTS: Distribution is skewed to the right; Skewness = 2.978002; 66.74% of variables have zero values.
### SUMMARY: Median imputation used to treat missing values.
###          No outliers present, the variable has a reasonable scale and hence no scaling is applied.
###          [slc_acct_pre_lim_perc_imputed] has scale [0;1] and 529862 observations with 746130 missing values; with 0 at 50% quantile and 0 at 75% quantile and mean of 0.09851
###          No further treatment necessary

# [slc_acct_prepaid_perc_dir_12] - Missing value indicator
datCredit_smp[, value_ind_slc_acct_prepaid_perc_dir_12 := ifelse(is.na(slc_acct_prepaid_perc_dir_12) | slc_acct_prepaid_perc_dir_12 == "", 0, 1)]
describe(datCredit_smp$value_ind_slc_acct_prepaid_perc_dir_12)
### RESTULS: Binning successful, no missing values in [value_ind_slc_acct_prepaid_perc_dir_12]

# [slc_acct_prepaid_perc_dir_12] - Missing value imputation
# describe(datCredit_smp$slc_acct_prepaid_perc_dir_12) ### AB: describe() stalled the process. disabled for now
### RESULTS:    [slc_acct_prepaid_perc_dir_12] has scale [0;103696000000000] and 746130 observations have missing values; with 0 at 50% quantile and 0.2061 at 75% quantile and mean of 19686957
### CONCLUSION: Median imputation proceeds from the distribution of values which constitutes performance spells only (all associated values from default spells are excluded)
# [TREATMENT] Median imputation for missing values
datCredit_smp[, slc_acct_prepaid_perc_dir_12_imputed_med := 
                ifelse(is.na(slc_acct_prepaid_perc_dir_12) | slc_acct_prepaid_perc_dir_12 == "", 
                       median(datCredit_smp[!is.na(PerfSpell_Key) & !(is.na(slc_acct_prepaid_perc_dir_12) | slc_acct_prepaid_perc_dir_12 == ""), slc_acct_prepaid_perc_dir_12], na.rm=TRUE), slc_acct_prepaid_perc_dir_12)]
#describe(datCredit_smp$slc_acct_prepaid_perc_dir_12_imputed_med) ### AB: describe() stalled the process. disabled for now
### RESULTS: [slc_acct_prepaid_perc_dir_12_imputed] has scale [0;103696000000000] and 0 observations have missing values; with 0 at 50% quantile and 0.0002421 at 75% quantile and mean of 17355944
hist(datCredit_smp$slc_acct_prepaid_perc_dir_12_imputed_med, breaks=500); skewness(datCredit_smp$slc_acct_prepaid_perc_dir_12_imputed_med, na.rm = T); datCredit_smp[slc_acct_prepaid_perc_dir_12_imputed_med==0,.N]/datCredit_smp[,.N]
### RESULTS:    Distribution is skewed to the right; Skewness = 2509.185; 74.8% of variables have zero values.

# [slc_acct_roll_ever_24] - Missing value indicator
datCredit_smp[, value_ind_slc_acct_roll_ever_24 := ifelse(is.na(slc_acct_roll_ever_24) | slc_acct_roll_ever_24 == "", 0, 1)]
describe(datCredit_smp$value_ind_slc_acct_roll_ever_24)
### RESULTS: Binning successful, no missing values in [value_ind_slc_acct_roll_ever_24]

# [slc_acct_roll_ever_24] - Missing value imputation
describe(datCredit_smp$slc_acct_roll_ever_24)
### RESULTS:    [slc_acct_roll_ever_24] has scale [0;4] and 746700 observations have missing values and mean of 0.304
### CONCLUSION: Median imputation proceeds from the distribution of values which constitutes performance spells only (all associated values from default spells are excluded)
# [TREATMENT] Median imputation for missing values
datCredit_smp[, slc_acct_roll_ever_24_imputed_med := 
                ifelse(is.na(slc_acct_roll_ever_24) | slc_acct_roll_ever_24 == "" | slc_acct_roll_ever_24 == "Unknown", 
                       median(datCredit_smp[!is.na(PerfSpell_Key) & !(is.na(slc_acct_roll_ever_24) | slc_acct_roll_ever_24 == "") ,slc_acct_roll_ever_24], na.rm=TRUE), slc_acct_roll_ever_24)]
describe(datCredit_smp$slc_acct_roll_ever_24_imputed_med)
### RESULTS: [slc_acct_roll_ever_24_imputed_med] has scale [0;4] and 0 observations have missing values and mean of 0.1815
hist(datCredit_smp[,slc_acct_roll_ever_24_imputed_med], breaks=500); skewness(datCredit_smp$slc_acct_roll_ever_24_imputed_med, na.rm = T); datCredit_smp[slc_acct_roll_ever_24_imputed_med==0,.N]/datCredit_smp[,.N]
### RESULTS: Distribution is skewed to the right; Skewness = 3.203737; 86.18% of variables have zero values.
### SUMMARY: Median imputation used to treat missing values.
###          No outliers present, since the variable represents counts and only has 4 levels.
###          No further treatment necessary

# - Remove the variables that have missingness > 50%
suppressWarnings( datCredit_smp[, `:=`(value_ind_slc_days_excess = NULL, slc_days_excess = NULL, 
                                        value_ind_ccm_ute_lvl_40_cnt_24m = NULL, ccm_ute_lvl_40_cnt_24m = NULL,
                                        value_ind_ccm_worst_arrears_6m = NULL, ccm_worst_arrears_6m = NULL,
                                        value_ind_ccm_worst_arrears_24m = NULL, ccm_worst_arrears_24m = NULL)]); gc()



# --- 4.2 Analysis and Treatments of Numeric variables
# - [Principal]
describe(datCredit_smp$Principal)
### RESULTS: [Principal] has scale [0.17;22370000], with 512428    at 50% quantile and 850000 at 75% quantile and mean of 652573.
hist(datCredit_smp$Principal, breaks='FD'); skewness(datCredit_smp$Principal, na.rm = T); datCredit_smp[Principal==0,.N]/datCredit_smp[,.N]
### RESULTS:    Distribution is skewed to the right; Skewness = 4.120052; 0% of observations have zero values.

# - [Instalment]
describe(datCredit_smp$Instalment)
### RESULTS: [Instalment] has scale [0;14244400], with 4728.9 at 50% quantile and 7976.7 at 75% quantile and mean of 5998
hist(datCredit_smp$Instalment, breaks ='FD'); skewness(datCredit_smp$Instalment, na.rm = T); datCredit_smp[Instalment==0,.N]/datCredit_smp[,.N]
### RESULTS:    Distribution is skewed to the right; Skewness = 434.135; 0.31% of observations have zero values.

# - [Balance]
describe(datCredit_smp$Balance)
### RESULTS: [Balance] has scale [-146.71;21738700], with 362199.88 at 50% quantile and 693429.75 at 75% quantile and mean of 491971
hist(datCredit_smp$Balance, breaks = 'FD'); skewness(datCredit_smp$Balance, na.rm = T); datCredit_smp[Balance==0,.N]/datCredit_smp[,.N]
### RESULTS:    Distribution is skewed to the right; Skewness = 3.385451; 3.77% of variables have zero values.

# - [InterestRate_Margin] (incorporating risk-based pricing info)
# Create [InterestRate_Margin] using the repo rate + 3.5% (Prime Rate's definition in South Africa)
datCredit_smp <- datCredit_smp %>% mutate(InterestRate_Margin = round(InterestRate_Nom - (M_Repo_Rate+0.035), digits=4)) %>%
  relocate(InterestRate_Margin, .after=InterestRate_Nom)
# Distributional analysis
#describe(datCredit_smp$InterestRate_Margin); ### AB: describe() stalled the process. disabled for now
hist(datCredit_smp$InterestRate_Margin, breaks="FD")
datCredit_smp[is.na(InterestRate_Margin), .N] / datCredit_smp[,.N] * 100
### RESULTS:    Highly right-skewed distribution (as expected), with mean of -0.007401 vs median of -0.008, 
###             bounded by [-0.02, 0.012] for 5%-95% percentiles; some negative outliers distort shape of distribution
### CONCLUSION: Use median imputation, given 0.53% missingness degree
datCredit_smp[, InterestRate_Margin_imputed_mean := 
                ifelse(is.na(InterestRate_Margin) | InterestRate_Margin == "", 
                       median(InterestRate_Margin, na.rm=TRUE), InterestRate_Margin)]
# [SANITY CHECK] Confirm treatment success
cat( ( datCredit_smp[is.na(InterestRate_Margin_imputed_mean), .N] == 0) %?% 
       'SAFE: Treatment successful for [InterestRate_Margin_imputed_mean].\n' %:% 
       'ERROR: Treatment failed for [InterestRate_Margin_imputed_mean] \n' )
#describe(datCredit_smp$InterestRate_Margin_imputed_mean); ### AB: describe() stalled the process. disabled for now
hist(datCredit_smp$InterestRate_Margin_imputed_mean, breaks="FD")
### RESULTS: Highly right-skewed distribution (as expected), with mean of -0.007405 vs median of -0.008, bounded by [-0.02, 0.12] for 5%-95% percentiles.

# - Save intermediary snapshot to disk (zip) for quick disk-based retrieval later
pack.ffdf(paste0(genPath,"creditdata_final_TFD_smp1b"), datCredit_smp)



# --- 4.3 Binning and factorisation
# - Confirm prepared datasets are loaded into memory
if (!exists('datCredit_smp')) unpack.ffdf(paste0(genPath,"creditdata_final_TFD_smp1b"), tempPath)

# - Condense the payment group
datCredit_smp[, pmnt_method_grp := 
                case_when(slc_pmnt_method == "Debit Order FNB account" | slc_pmnt_method == "Debit Order other bank" ~ "Debit Order",
                          slc_pmnt_method == "Salary" | slc_pmnt_method == "Suspense" ~ "Salary/Suspense",
                          TRUE ~ slc_pmnt_method)]
# [SANITY CHECK] Check new feature for illogical values
cat( ( datCredit_smp[is.na(pmnt_method_grp), .N] == 0) %?% 
       'SAFE: New feature [pmnt_method_grp] has logical values.\n' %:% 
       'WARNING: New feature [pmnt_method_grp] has illogical values \n' )
describe(datCredit_smp$pmnt_method_grp)
### RESULTS: Bins grouped logically such that each bin now has sufficient observations

# - [PerfSpellResol_Type_Hist]
describe(datCredit_smp$PerfSpellResol_Type_Hist)
datCredit_smp$PerfSpellResol_Type_Hist %>% table() %>% prop.table()
barplot(table(datCredit_smp$PerfSpellResol_Type_Hist))
### RESULTS [PerfSpellResol_Type_Hist] has 5 levels and no missing values. 12.1% of performance spells defaulted; 41.71% of performance spells have been right-censored.
# [TREATMENT] Apply factor transformation
datCredit_smp[,PerfSpellResol_Type_Hist:=factor(PerfSpellResol_Type_Hist)]
### RESULTS: No further treatment necessary

# - Factorised [g0_Delinq] variable
datCredit_smp[,g0_Delinq_fac := as.factor(g0_Delinq)]
describe(datCredit_smp$g0_Delinq_fac)
### RESULTS: proportion in 0 equals 94%, proportion in 1 equals 5.1%, proportion in 2 equals 0.6%, proportion in 3 equals 0.3%

# - Bin [InterestRate_Margin_imputed] | Binning the variable into three equally sized bins
datCredit_smp[, InterestRate_Margin_imputed_bin := factor(ntile(InterestRate_Margin_imputed_mean, n=3))]
describe(datCredit_smp$InterestRate_Margin_imputed_bin)



# --- 4.4 Feature Engineering: ratio-type variables (Period-level)
# - [AgeToTerm]
# Create [AgeToTerm] variable
datCredit_smp <- datCredit_smp %>% mutate(AgeToTerm = Age_Adj/Term)
#describe(datCredit_smp$AgeToTerm) ### AB: describe() stalled the process. disabled for now
### RESULTS: [AgeToTerm] has scale [0.00185185 ;58.3333], with 0.2875 at 50% quantile and 0.5333 at 75% quantile and mean of 0.3621
hist(datCredit_smp$AgeToTerm, breaks='FD'); skewness(datCredit_smp$AgeToTerm, na.rm = T); datCredit_smp[AgeToTerm==0,.N]/datCredit_smp[,.N]
### RESULTS:    Distribution is skewed to the right; Skewness = 16.15954; 0% of variables have zero values.

# - [BalanceToTerm]
# Create [BalanceToTerm] variable
datCredit_smp <- datCredit_smp %>% mutate(BalanceToTerm = Balance/Term)
# Distributional analysis
describe(datCredit_smp$BalanceToTerm) ### AB: describe() stalled the process. disabled for now
## RESULTS: [BalanceToTerm] has scale [-0.611292;90577.9], with 1531.3769 at 50% quantile and 2919.7919 at 75% quantile and mean of 2081
hist(datCredit_smp$BalanceToTerm, breaks='FD'); skewness(datCredit_smp$BalanceToTerm, na.rm = T); datCredit_smp[BalanceToTerm==0,.N]/datCredit_smp[,.N]
### RESULTS:    Distribution is skewed to the right; Skewness = 3.769301; 3.8% of variables have zero values.

# - Save intermediary snapshot to disk (zip) for quick disk-based retrieval later
pack.ffdf(paste0(genPath,"creditdata_final_TFD_smp1c"), datCredit_smp)



# --- 4.5 Featuring Engineering: Portfolio-level information
# - Confirm prepared datasets are loaded into memory
if (!exists('datCredit_smp')) unpack.ffdf(paste0(genPath,"creditdata_final_TFD_smp1c"), tempPath)

# - Pre default delinquency rate
#Note: Creating an aggregated dataset with which to fuse to the full dataset
dat_g0_Delinq_Aggr <- data.table(datCredit_smp[!is.na(PerfSpell_Key) & Default_Ind==0, list(sum(g0_Delinq>0, na.rm=T)/.N), by=list(Date)])
colnames(dat_g0_Delinq_Aggr) <- c("Date", "g0_Delinq_Any_Aggr_Prop")
# Applying various lags
lags <- c(1,2,3,4,5,6,9,12) # Lags
ColNames <- colnames(dat_g0_Delinq_Aggr)[-1] # Names of the columns
for (i in seq_along(lags)){ # Looping over the specified lags and applying each to each of the specified columns
  for (j in seq_along(ColNames)){
    dat_g0_Delinq_Aggr[, (paste0(ColNames[j],"_Lag_",lags[i])) := fcoalesce(shift(get(ColNames[j]), n=lags[i], type="lag"),get(ColNames[j]))] # Impute NA's with the non lagged value
  }
}
# [Sanity Check] Check for any missing values before merging the dat_g0_Delinq_Aggr dataset to datCredit_smp
cat((anyNA(dat_g0_Delinq_Aggr)) %?% 'WARNING: One of the new [g0_Delinq_Any_Aggr_Prop] features has missing values. \n' %:%
      'SAFE: New [g0_Delinq_Any_Aggr_Prop] features created sucessfully without any missing values. \n')
### RESULTS: [g0_Delinq_Any_Aggr_Prop] variables created successfully without any missingness
# Fusing the aggregated variable with its various lags to the full dataset
datCredit_smp <- merge(datCredit_smp, dat_g0_Delinq_Aggr, by="Date", all.x=T)
# [SANITY CHECK] Check new feature for illogical values
cat( ( sum(datCredit_smp[DefaultStatus1==0, sum(g0_Delinq_Any_Aggr_Prop + sum(g0_Delinq==0)/.N, na.rm=T), by=Date][,2])==sum(datCredit_smp[DefaultStatus1==0,.N,by=Date][,2]) & (sum(is.na(datCredit_smp$g0_Delinq_Any_Aggr_Prop))==0)) %?% 
       'SAFE: New feature [g0_Delinq_Any_Aggr_Prop] has logical values.\n' %:% 
       'WARNING: New feature [g0_Delinq_Any_Aggr_Prop] has illogical values \n' )
describe(datCredit_smp$g0_Delinq_Any_Aggr_Prop); plot(unique(datCredit_smp$g0_Delinq_Any_Aggr_Prop), type="b")
### RESULTS: Variable has a logical trend, with mean of 0.05828 vs median of 0.04720, 
# bounded by [0.03926, 0.11895] for 5%-95% percentiles; no large outliers
#PROBLEM

# [SANITY CHECK] Check new feature for missingness after fusion
cat((anyNA(datCredit_smp$g0_Delinq_Any_Aggr_Prop)) %?% 'WARNING: New feature [g0_Delinq_Any_Aggr_Prop] has missing values. \n' %:%
      'SAFE: New feature [g0_Delinq_Any_Aggr_Prop] has no missing values. \n')
### RESULTS: [g0_Delinq_Any_Aggr_Prop] created without any missingness


# - Average pre-default delinquency level
datCredit_smp[,g0_Delinq_Ave:=mean(ifelse(DefaultStatus1==0,g0_Delinq,0), na.rm=T), by=Date]
# [SANITY CHECK] Check new feature for illogical values
cat( (sum(datCredit_smp[, sum(is.na(g0_Delinq_Ave)), by=Date][,2])==0) %?% 
       'SAFE: New feature [g0_Delinq_Ave] has logical values.\n' %:% 
       'WARNING: New feature [g0_Delinq_Ave] has illogical values \n' )
#describe(datCredit_smp$g0_Delinq_Ave);  ### AB: describe() stalled the process. disabled for now
hist(datCredit_smp$g0_Delinq_Ave, breaks="FD")
### RESULTS: Follows a logical trend, with mean of 0.06397 vs median of 0.05269,
# bounded by [0.04357, 0.12998] for 5%-95% percentiles; no outliers


# - Ratio type variables (portfolio-level) during performance spells
# (Total) Arrears to (Total) Balance; (Total) Instalments to (Total) Balance
# NOTE: These portfolio-level aggregated variables are engineered to capture/ aggregate information only for accounts that are in a performance spell
# The resulting aggregated dataset can be fused to the full dataset
dat_Aggr <- data.table(datCredit_smp[DefaultStatus1==0, list(sum(Arrears, na.rm=T)/sum(Balance, na.rm=T)), by=list(Date)], # [ArrearsToBalance_Aggr]
                       datCredit_smp[DefaultStatus1==0, list(sum(Instalment, na.rm=T)/sum(Balance)), by=list(Date)][,2]) # [InstalmentToBalance_Aggr]
colnames(dat_Aggr) <- c("Date", "ArrearsToBalance_Aggr_Prop", "InstalmentToBalance_Aggr_Prop")
# Fusing the aggregated dataset to the full dataset
datCredit_smp <- merge(datCredit_smp, dat_Aggr, by="Date", all.x=T)
# [SANITY CHECK] Check new feature for illogical values
cat( (sum(datCredit_smp[, sum(is.na(ArrearsToBalance_Aggr_Prop)), by=Date][,2])==0) %?% 
       'SAFE: New feature [ArrearsToBalance_Aggr_Prop] has logical values.\n' %:% 
       'WARNING: New feature [ArrearsToBalance_Aggr_Prop] has illogical values \n' )
cat( (sum(datCredit_smp[, sum(is.na(InstalmentToBalance_Aggr_Prop)), by=Date][,2])==0) %?% 
       'SAFE: New feature [InstalmentToBalance_Aggr_Prop] has logical values.\n' %:% 
       'WARNING: New feature [InstalmentToBalance_Aggr_Prop] has illogical values \n' )
describe(datCredit_smp$InstalmentToBalance_Aggr_Prop); plot(unique(datCredit_smp$Date),unique(datCredit_smp$InstalmentToBalance_Aggr_Prop), type="b")
#describe(datCredit_smp$ArrearsToBalance_Aggr_Prop); ### AB: describe() stalled the process. disabled for now
plot(unique(datCredit_smp$ArrearsToBalance_Aggr_Prop), type="b")
### RESULTS [InstalmentToBalance_Aggr_Prop]: Variable has high volatility around 2010 as seen through the graphical plot. Mean of 0.01228 vs median of 0.01207,
#            bounded by [0.01088, 0.01422] for 5%-95% percentiles; no outliers
#                   [ArrearsToBalance_Aggr_Prop]: Variable has  mean of 0.0006134 vs median of 0.0004895,
#                    bounded by [0.0003790, 0.0015335] for 5%-95% percentiles


# - Proportion of curing loans across performing/default spell type
datCredit_smp[, CuringEvents_Aggr_Prop := sum(PerfSpell_Counter==1 & PerfSpell_Num>=2, na.rm=T)/.N, by=list(Date)]
cat( (sum(datCredit_smp[, sum(is.na(CuringEvents_Aggr_Prop)), by=Date][,2])==0) %?% 
       'SAFE: New feature [CuringEvents_Aggr_Prop] has logical values.\n' %:% 
       'WARNING: New feature [CuringEvents_Aggr_Prop] has illogical values \n' )
#describe(datCredit_smp$CuringEvents_Aggr_Prop); ### AB: describe() stalled the process. disabled for now
plot(unique(datCredit_smp$CuringEvents_Aggr_Prop), type="b")
### RESULTS: Variable has mean of 0.001373 vs median of 0.0012615,
# bounded by [0.0006567, 0.0026194] for 5%-95% percentiles; no outliers

# Clean up
rm(lags, ColNames, dat_g0_Delinq_Aggr, dat_Aggr)


# - Aggregated age-to-term of portfolio over time, i.e., percentage-based maturity
datCredit_smp[, AgeToTerm_Aggr_Mean := mean(Age_Adj/Term, na.rm=T), by=Date]
cat( (sum(datCredit_smp[, sum(is.na(AgeToTerm_Aggr_Mean)), by=Date][,2])==0) %?% 
       'SAFE: New feature [AgeToTerm_Aggr_Mean] has logical values.\n' %:% 
       'WARNING: New feature [AgeToTerm_Aggr_Mean] has illogical values \n' )
describe(datCredit_smp$AgeToTerm_Aggr_Mean); plot(unique(datCredit_smp$AgeToTerm_Aggr_Mean), type="b")
### RESULTS: Variable behaves as expected, i.e., increases as the loan portfolio matures. Has mean 0.3621 and median 0.3878
# bounded by [0.2568, 0.4006] for 5%-95% percentiles; no outliers


# - Aggregate maturity of performance spell ages over time
datCredit_smp[, PerfSpell_Maturity_Aggr_Mean := mean(PerfSpell_Age, na.rm=T), by=Date]
cat( (sum(datCredit_smp[, sum(is.na(PerfSpell_Maturity_Aggr_Mean)), by=Date][,2])==0) %?% 
       'SAFE: New feature [PerfSpell_Maturity_Aggr_Mean] has logical values.\n' %:% 
       'WARNING: New feature [Perf_SpellMaturity_Aggr_Mean] has illogical values \n' )
describe(datCredit_smp$PerfSpell_Maturity_Aggr_Mean); plot(unique(datCredit_smp$PerfSpell_Maturity_Aggr_Mean), type="b")
### RESULTS: Mean performance spell age seem to decrease over time. Has mean 135.1 and median 140.68;
# bounded by [93.54, 152.44] for 5%-95% percentiles; no outliers


# - Median-aggregated interest rate margin
# NOTE: The median is preferred over the mean since it resulted in a superior model, as investigated in the experimental script 3c(v)
# Creating an aggregated dataset
dat_IRM_Aggr <- datCredit_smp[, list(InterestRate_Margin_Aggr_Med = median(InterestRate_Margin_imputed_mean, na.rm=T)), by=list(Date)]
gc()
# Checking the time series of this variable
plot(dat_IRM_Aggr$InterestRate_Margin_Aggr_Med, type="b")
# Applying various lags
lags <- c(1,2,3,4,5,6,9,12) 
dat_IRM_Aggr_Check1 <- data.table(Variable = NULL, # Dataset for conducting sanity checks
                                  Check = NULL)
ColNames <- colnames(dat_IRM_Aggr)[-1] # Names of the columns
for (i in seq_along(lags)){ # Looping over the specified lags and applying each to each of the specified columns
  for (j in seq_along(ColNames)){
    dat_IRM_Aggr[, (paste0(ColNames[j],"_",lags[i])) := fcoalesce(shift(get(ColNames[j]), n=lags[i], type="lag"),get(ColNames[j]))] # Impute NA's with non-lagged version of variable
  }
}
# [SANITY CHECK] Check whether the lags were created correctly
cat((anyNA(dat_IRM_Aggr[,InterestRate_Margin_Aggr_Med_1]) | anyNA(dat_IRM_Aggr[,InterestRate_Margin_Aggr_Med_2]) | anyNA(dat_IRM_Aggr[,InterestRate_Margin_Aggr_Med_3]) | anyNA(dat_IRM_Aggr[,InterestRate_Margin_Aggr_Med_9])) %?%
      "WARNING: Missingness detected, [InterestRate_Margin_Aggr_Med_1], [InterestRate_Margin_Aggr_Med_2] and/or [InterestRate_Margin_Aggr_Med_3] compromised.\n" %:%
      "SAFE: No missingness, [InterestRate_Margin_Aggr_Med_1], [InterestRate_Margin_Aggr_Med_2] and [InterestRate_Margin_Aggr_Med_3] created successfully.\n")
### RESULTS: Safe, no missingness, hence continue with merge

# Merging the credit dataset with the aggregated dataset
datCredit_smp <- merge(datCredit_smp, dat_IRM_Aggr, by="Date", all.x=T)
# Validate merging success )by checking for missingness (should be zero)
list_merge_variables <- list(colnames(dat_IRM_Aggr))
results_missingness <- list()
for (i in 1:length(list_merge_variables)){
  output <- sum(is.na(datCredit_smp$list_merge_variables[i]))
  results_missingness[[i]] <- output
}
cat( (length(which(results_missingness > 0)) == 0) %?% "SAFE: No missingness, fusion with aggregated data is successful.\n" %:%
       "WARNING: Missingness in certain aggregated fields detected, fusion compromised.\n")
describe(datCredit_smp$InterestRate_Margin_Aggr_Med); plot(datCredit_smp[!duplicated(Date),InterestRate_Margin_Aggr_Med], type="b") # Only saving the base variable's descriptive statistics
### RESULTS: Variable follows a logical trend over time. Has mean -0.008077 and median -0.0085;
# bounded by [-0.012, -0.0040] for 5%-95% percentiles; no outliers

# - Save final snapshot to disk (zip) for quick disk-based retrieval later
pack.ffdf(paste0(genPath,"creditdata_final_TFD_smp2"), datCredit_smp)

# Clean up
rm(dat_IRM_Aggr, dat_IRM_Aggr_Check1, list_merge_variables, results_missingness, output, lags, ColNames,varSLC_Info_Cat, varSLC_Info_Num, varCredit_Info_Cat, varCredit_Info_Num, check.fuse1, check.fuse3, check.fuse4, lookup_IDs,
   Covariate_Info, lookup, lookup2);gc()



# ------ 5. Implementing a simple (loan-level) resampling scheme
# --- 5.1 Apply resampling
# - Set seed
set.seed(1, kind="Mersenne-Twister")

# - Use simple random sampling with the stratifiers to select the loan IDs that ought to be in the training dataset
if (all(!is.na(stratifiers))){ # Stratifiers
  dat_train_keys <- dat_smp_keys %>% group_by(across(all_of(stratifiers))) %>% slice_sample(prop=smp_frac) %>% as.data.table() 
} else { # No stratifiers
  dat_train_keys <- dat_smp_keys %>% slice_sample(prop=smp_frac) %>% as.data.table()
}

# - Obtain the associated loan records as to create the training dataset
datCredit_train_TFD <- copy(datCredit_smp[get(clusVar) %in% dat_train_keys[, get(clusVar)],]);gc()
# - Ensure that the training set only contains the first performance spell
datCredit_train_TFD <- datCredit_train_TFD[PerfSpell_Num == 1,]

# - Obtain the associated loan records of the validation dataset
datCredit_valid_TFD <- copy(datCredit_smp[!(get(clusVar) %in% dat_train_keys[, get(clusVar)]),]);gc()

# - [SANITY CHECK] Reconciling the cardinalities of the training and the validation dataset to the subsampled dataset
check.2 <- datCredit_smp[,.N] == datCredit_train_TFD[,.N] + datCredit_valid_TFD[,.N] # Should be TRUE
ifelse(check.2, print('SAFE: Training and validation datasets reconstitue subsampled dataset'),
                print('WARNING: Training and validation datasets does not reconstitue subsampled dataset \n' ))

# - [SANITY CHECK] Checking if all account-specific observations are clustered together
check.3 <- length(unique(datCredit_smp$LoanID)) == length(unique(datCredit_train_TFD$LoanID)) + length(unique(datCredit_valid_TFD$LoanID))
ifelse(check.3, print('SAFE: All associated loan observations are clustered together in the training and validation datasets; no cross-contamination detected '),
                print('WARNING: Not all associated loan observations are clustered together in the training and validation datasets; cross-contamination may exist \n' ))

# - Clean up
suppressWarnings(rm(smp_perc, dat_keys, dat_smp_keys, dat_train_keys, check.2, check.3, datCredit_smp_old_n, datCredit_smp_prior, dat_keys_exc, class_type, excCond, excCond2, datCredit_smp))


# --- 5.2 Saving the cross-validation scheme
# - Training dataset
pack.ffdf(paste0(genPath,"creditdata_train_TFD"), datCredit_train_TFD)

# - Validation dataset
pack.ffdf(paste0(genPath,"creditdata_valid_TFD"), datCredit_valid_TFD)


# --- 5.3 Clean up
suppressWarnings(rm(dat_keys_smp_perf, dat_keys_smp_perf,  dat_train_keys_perf, dat_train_keys_def, datCredit_train_perf, datCredit_train_def,  datCredit_valid_perf, datCredit_valid_def,
                    check.4_a, check.4_b, check.4_c, check.5_a, check.5_b, datCredit_smp, datStrata_smp_min, datCredit_train_TFD, datCredit_valid_TFD));gc()
