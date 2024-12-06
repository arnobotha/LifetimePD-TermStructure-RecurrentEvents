# ========================= INVESTIGATING SURVIVALROC() ==========================
# As a poof of concept, we shall compare various functions from various packages
# in conducting time-dependent ROC-analyses on the same fitted Cox regression model
# --------------------------------------------------------------------------------
#  PROJECT TITLE: Default Survival Modelling
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
# ================================================================================





# --------- 1. Preliminaries
# -- Confirm prepared database is loaded into memory
if (!exists('datCredit_real')) unpack.ffdf(paste0(genPath,"creditdata_final4a"), tempPath)





# --------- 2. Prepare data according to time definitions


# --- 2.1a Time to First Default time definition
# NOTE: The first spell condition is only enforced during the resampling of data into the training set, since we purposefully
# would like the validation set to include multiple spells in order to validate certain assumptions.
# Dplyr-actions include in order the following tasks:
# 1) Rename and adjust start and ending times for survival modelling
# 2) Subset for performance spells only
datCredit_TFD1 <- datCredit_real %>% mutate( Start = ifelse(is.na(PerfSpell_Counter),NA,PerfSpell_Counter-1), # Records the start of time interval
                                             End = ifelse(is.na(PerfSpell_Counter),NA,PerfSpell_Counter), # Records the end of time interval
                                             Default_Ind = ifelse(!is.na(PerfSpell_Num) & DefaultStatus1==1,1,
                                                                  ifelse(!is.na(PerfSpell_Num),0,NA))) %>%
  filter(!is.na(PerfSpell_Num)) # Filter for only performance spells
# Sanity check - Should be TRUE
datCredit_TFD1[is.na(PerfSpell_Num),.N] == 0  # TRUE, field created successfully




# --- 2.1b Time to First Default time definition | Alternate version from AB
# NOTE: The first spell condition is only enforced during the resampling of data into the training set, since we purposefully
# would like the validation set to include multiple spells in order to validate certain assumptions.
# - Dplyr-actions include in order the following tasks:
# 1) Subset for performance spells only
# 2) Rename and adjust start and ending times for survival modelling
datCredit_TFD2 <- subset(datCredit_real, !is.na(PerfSpell_Num)) %>%
  mutate(Start = PerfSpell_Counter-1, End = PerfSpell_Counter,
         Default_Ind = DefaultStatus1)
### AB: Using the "_Counter" variety is incorrect since it refers simply to the row index within a spell, i.e.,
# 1, ..., <Spell_Length>. Instead, one should rather use the "timeInPerfSpell" as the spell period since it
# accounts for left-truncated spells, which are quite prevalent in our dataset.
# We should fix this after the close-out and rerun all subsequent analyses/models downstream from this point.
# It is an important update and the results should shift on some axis, though will perhaps not change in
# distributional shape.
# Sanity check - Should be TRUE
datCredit_TFD2[is.na(PerfSpell_Num),.N] == 0  # TRUE, field created successfully




# --- 2.1c Time to First Default time definition | Alternate and corrected version from AB
# NOTE: The first spell condition is only enforced during the resampling of data into the training set, since we purposefully
# would like the validation set to include multiple spells in order to validate certain assumptions.
# - Dplyr-actions include in order the following tasks:
# 1) Subset for performance spells only
# 2) Rename and adjust start and ending times for survival modelling
datCredit_TFD3 <- subset(datCredit_real, !is.na(PerfSpell_Num)) %>%
  mutate(Start = TimeInPerfSpell-1, End = TimeInPerfSpell,
         Default_Ind = DefaultStatus1)
datCredit_TFD3[is.na(PerfSpell_Num),.N] == 0  # TRUE, field created successfully



# --------- 3. Investigate similarities between datasets.
# -- Compare similarity between datasets 1) datCredit_TFD1 and 2) datCredit_TFD2.
all.equal(datCredit_TFD1, datCredit_TFD2) # TRUE
# datCredit_TFD1 and datCredit_TFD2 are indeed identical and therefore the order of operation did not change final form of table.

# -- Compare similarity between datasets 1) datCredit_TFD1 and 2) datCredit_TFD3.
all.equal(datCredit_TFD1, datCredit_TFD3) # Not TRUE
# datCredit_TFD1 and datCredit_TFD2 are not identical and the differences appear to be in the Start column

# - Investigate whether differences are due to truncated loans.
all.equal(datCredit_TFD1[PerfSpell_LeftTrunc==0,], datCredit_TFD3[PerfSpell_LeftTrunc==0,]) # TRUE
# datCredit_TFD1 and datCredit_TFD2 are identical if truncated performance spells are filtered out.





# --------- 3. Determine whether the difference between datCredit_TFD1 and datCredit_TFD2 are significant
# Build coxph model based on datCredit_TFD1 and fit on random variable but compare base hazard rates.
cox_TFD1 <- coxph(Surv(Start,End,Default_Ind) ~ g0_Delinq_SD_4, id=LoanID, data=datCredit_TFD1)

# Build coxph model based on datCredit_TFD1 and fit on random variable but compare base hazard rates.
cox_TFD3 <- coxph(Surv(Start,End,Default_Ind) ~ g0_Delinq_SD_4, id=LoanID, data=datCredit_TFD3)



### AB: Compare the differences among TFD1-3 and establish whether or not these 
# are significant by building different Cox regression models from them.
# Method of comparison: time-dependent ROC-analysis.
