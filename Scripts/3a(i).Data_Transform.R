# =========================================== Data enrichment ===========================================
# Prepare extractions from the base credit dataset towards creating various time definitions, themselves
#  used in subsequent Cox regression-based modelling of the default term-structure
# ------------------------------------------------------------------------------------------------------
# PROJECT TITLE: Default survival modelling
# SCRIPT AUTHOR(S): Bernard Scheepers (BS), Dr Arno Botha (AB)

# DESCRIPTION:
# This script performs the following high-level tasks:
#   1) Extract a smaller subset of the main credit dataset
#   2) Engineer additional features in assisting with analyses and modelling later
#   3) Extract, prepare, and save subsets of the dataset according to the following time definitions: 
#       a) Time to first default; b) Anderson-Gill (AG) method; c) Prentice-Williams-Peterson (PWP) method 
#         based on Total-time approach; d) PWP-method based on Spell-time approach.
#
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

# -- Inputs:
#   - datCredit_real | Prepared from script 2b.
#
# -- Outputs:
#   - datCredit_TFD
#   - datCredit_AG
#   - datCredit_PWP_TT
#   - datCredit_PWP_ST
# ------------------------------------------------------------------------------------------------------




# --------- 1. Preliminaries
# -- Confirm prepared database is loaded into memory
if (!exists('datCredit_real')) unpack.ffdf(paste0(genPath,"creditdata_final4a"), tempPath)

Defcol <- c("DefSpell_Num", "TimeInDefSpell") # Default variables to remove
ncol <- ncol(datCredit_real)
datCredit_real[, (Defcol) := NULL]

# Sanity Check - TRUE
ncol(datCredit_real) + length(Defcol) == ncol




# --------- 2. Additional (light) feature engineering towards implementing various time definitions
# - Max date of each performance spell (used as a stratifier)
datCredit_real[!is.na(PerfSpell_Key), PerfSpell_Max_Date := max(Date, na.rm=T), by=list(PerfSpell_Key)]
datCredit_real[is.na(PerfSpell_Key), PerfSpell_Max_Date:= NA]

# - Min date of each performance spell (used as a stratifier)
datCredit_real[!is.na(PerfSpell_Key), PerfSpell_Min_Date := min(Date, na.rm=T), by=list(PerfSpell_Key)]
datCredit_real[is.na(PerfSpell_Key), PerfSpell_Min_Date:= NA]

# - Create an indicator variable variable for when a loan exits a performance spell
datCredit_real[,PerfSpell_Exit_Ind := ifelse(Date==PerfSpell_Max_Date,1,0)]
#datCredit_real[,PerfSpell_Exit_Ind2 := with(datCredit_real, ave(seq_along(PerfSpell_Key), PerfSpell_Key, FUN = function(x) x == max(x)))]
#all.equal(datCredit_real[!is.na(PerfSpell_Key), PerfSpell_Exit_Ind], datCredit_real[!is.na(PerfSpell_Key), PerfSpell_Exit_Ind2])
### RESULTS: TRUE

# - Create new spell resolution types by grouping competing risks together & relocating variable
# Performance spells
datCredit_real <- datCredit_real %>% mutate(PerfSpellResol_Type_Hist2 = case_when(PerfSpellResol_Type_Hist=="Defaulted" ~ "Defaulted",
                                                                                  PerfSpellResol_Type_Hist=="Censored" ~ "Censored",
                                                                                  TRUE ~ "Settled & Other")) %>%
  relocate(PerfSpell_Exit_Ind, PerfSpell_Max_Date, PerfSpell_Min_Date,  .after=PerfSpellResol_Type_Hist)

# - Checking the proportions of the newly created variable
describe(datCredit_real$PerfSpellResol_Type_Hist2)




# --------- 3. Prepare data according to time definitions

# --- 3.1 Time to First Default time definition
# Select performance spells only and create timing variables
datCredit_TFD <- subset(datCredit_real, !is.na(PerfSpell_Num)) %>%
  mutate(Start = TimeInPerfSpell-1, End = TimeInPerfSpell,
         Default_Ind = DefaultStatus1)
datCredit_TFD[is.na(PerfSpell_Num),.N] == 0  # TRUE, field created successfully
# - Save snapshots to disk (zip) for quick disk-based retrieval later
pack.ffdf(paste0(genPath,"creditdata_final_TFD"), datCredit_TFD)
# - Remove from memory as an expedient
rm(datCredit_TFD); gc()


# --- 3.2 Anderson-Gill (AG) time definition
# Select performance spells only and create timing variables
datCredit_AG <- subset(datCredit_real, !is.na(PerfSpell_Num)) %>%
  mutate(Start = Age_Adj-1, End = Age_Adj,
         Default_Ind = DefaultStatus1)
datCredit_AG[is.na(PerfSpell_Num),.N] == 0  # TRUE, field created successfully
# - Save snapshots to disk (zip) for quick disk-based retrieval later
pack.ffdf(paste0(genPath,"creditdata_final_AG"), datCredit_AG)
# - Remove from memory as an expedient
rm(datCredit_AG); gc()


# - 3.3 Prentice-Williams-Peterson (PWP) Spell-time definition
datCredit_PWPST <- subset(datCredit_real, !is.na(PerfSpell_Num)) %>%
  mutate(Start = TimeInPerfSpell-1, End = TimeInPerfSpell,
         Default_Ind = DefaultStatus1)
datCredit_PWPST[is.na(PerfSpell_Num),.N] == 0  # TRUE, field created successfully
# - Save snapshots to disk (zip) for quick disk-based retrieval later
pack.ffdf(paste0(genPath,"creditdata_final_PWPST"), datCredit_PWPST)
# - Remove from memory as an expedient
rm(datCredit_PWPST); gc()




# --------- 4. Housekeeping
suppressWarnings(rm(datCredit_real, exclusions_all, exclusions_credit)); gc()
