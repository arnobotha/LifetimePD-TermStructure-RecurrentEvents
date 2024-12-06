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





# --------- 1.1 Preliminaries
# -- Confirm prepared database is loaded into memory
if (!exists('datCredit_real')) unpack.ffdf(paste0(genPath,"creditdata_final4a"), tempPath)



# --------- 1.2 Prepare data according to time definitions

# --- 1.2a Time to First Default time definition
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

# --- 1.2b Time to First Default time definition | Alternate version from AB
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

# --- 1.3c Time to First Default time definition
# NOTE: The first spell condition is only enforced during the resampling of data into the training set, since we purposefully
# would like the validation set to include multiple spells in order to validate certain assumptions.
# - Dplyr-actions include in order the following tasks:
# 1) Subset for performance spells only
# 2) Rename and adjust start and ending times for survival modelling
datCredit_TFD3 <- subset(datCredit_real, !is.na(PerfSpell_Num)) %>%
  mutate(Start = TimeInPerfSpell-1, End = TimeInPerfSpell,
         Default_Ind = DefaultStatus1)
datCredit_TFD3[is.na(PerfSpell_Num),.N] == 0  # TRUE, field created successfully





# --------- 2. How similar are the 3 data tables?
# -- Compare similarity between datasets 1) datCredit_TFD1 and 2) datCredit_TFD2.
all.equal(datCredit_TFD1, datCredit_TFD2) # TRUE
# datCredit_TFD1 and datCredit_TFD2 are indeed identical and therefore the order of operation did not change final form of table.

# -- Compare similarity between datasets 1) datCredit_TFD1 and 2) datCredit_TFD3.
all.equal(datCredit_TFD1, datCredit_TFD3) # Not TRUE
# datCredit_TFD1 and datCredit_TFD3 are not identical and the differences appear to be in the Start column.

# --- 2.1 Is TFD1 and TFD3 different due to left-truncated spells?
all.equal(datCredit_TFD1[PerfSpell_LeftTrunc==0,], datCredit_TFD3[PerfSpell_LeftTrunc==0,]) # TRUE
### RESULTS: datCredit_TFD1 and datCredit_TFD3 are identical if truncated performance spells are filtered out.

# --- 2.2 How different are TFD1 and TFD3?
all.equal(datCredit_TFD1[,PerfSpell_Key], datCredit_TFD3[,PerfSpell_Key]) # TRUE
### RESUTLS: The left-trucated spells has no influence on the [PerfSpell_Key].

all.equal(datCredit_TFD1[,Start], datCredit_TFD3[,Start]) # FALSE
### RESUTLS: The left-trucated spells has an influence on the [PerfSpell_Key].

# -- Investigate distribution of differences
# Perform distributional analysis on difference in [Start] variables.
describe(datCredit_TFD3$Start - datCredit_TFD2$Start)
# Mean of 23.69 with a percentile distribution [0.05, 0.95] of [0; 138]

hist(datCredit_TFD3$Start - datCredit_TFD2$Start, breaks='FD')
# Extremely skewed distribution with significant extreme values (a maximum of 1216).

# -- Investigate distribution of differences given the loans are left-truncated.
# Perform distributional analysis on difference in [Start] variables.
describe(datCredit_TFD3[PerfSpell_LeftTrunc==1, Start] - datCredit_TFD1[PerfSpell_LeftTrunc==1, Start])
# Mean of 52.34 with a percentile distribution [0.05, 0.95] of [2; 167], with a minimum of 1

hist(datCredit_TFD3[PerfSpell_LeftTrunc==1, Start] - datCredit_TFD1[PerfSpell_LeftTrunc==1, Start], breaks='FD')
# Skewed distribution with significant extreme values

### CONCLUSION: There is a significant difference between the two time definitions, specifically as they relate to left-truncated performance spells.

# --------- 3. How prevalent is left-truncated performance spells
datCredit_TFD1[PerfSpell_LeftTrunc==1 & !duplicated(PerfSpell_Key), .N]/datCredit_TFD1[ !duplicated(PerfSpell_Key), .N]
# 0.373624

### CONCLUSION: Due to the high prevalence of left truncated performance spells, TFD3 should be implemented in favour of TFD1/TFD2.









