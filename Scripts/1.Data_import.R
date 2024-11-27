# ===================================== DATA IMPORT =====================================
# Import credit data and macroeconomic data
# ---------------------------------------------------------------------------------------
# PROJECT TITLE: Default survival modelling
# SCRIPT AUTHOR(S): Dr Arno Botha, Roelinde Bester

# DESCRIPTION:
# This script imports the loan performance credit dataset (SAS) into R, as well
# as macroeconomic datasets, and a large dataset containing various input fields
# ---------------------------------------------------------------------------------------
# -- Script dependencies:
#   - 0.Setup.R

# -- Inputs:
#   - creditdata_final_sub | monthly loan performance credit data (FNB Mortgages) (subject-period format)
#   - macro_data_monthly.sas7bdat | monthly macroeconomic data
#   - macro_data_quarterly.sas7bdat | quarterly macroeconomic data
#   - creditdata_input.sas7bdat | input fields associated with credit data (FNB Mortgages)

# -- Outputs:
#   - dat.raw
#   - macro_data_m
#   - macro_data_q
#   - datInput.raw
# =======================================================================================



# --------------------------------- IMPORTS ------------------------------------

# --------- Importing local data from various sources (using the haven package)
# Note: local sources should ideally be uncompressed SAS files, or at least 
# compressed using the COMPRESS=CHAR SAS-option.

# --- 1. Macroeconomic history + forecasts from FNB Group Economics | monthly data
# Import, then recast data into a more pliable data.table object for greater memory efficiency, during which the
# key is set based on preliminary analysis on the data grain (itself tested later again)
macro_data_m <- as.data.table(read_sas(paste0(genRawPath,"macro_data_monthly.sas7bdat")), stringsAsFactors=T,
                              key=c("EffectiveDate", "Scenario"))

# --- 2. Macroeconomic history + forecasts from FNB Group Economics | quarterly data
# Import, then recast data into a more pliable data.table object for greater memory efficiency, during which the
# key is set based on preliminary analysis on the data grain (itself tested later again)
macro_data_q <- as.data.table(read_sas(paste0(genRawPath,"macro_data_quarterly.sas7bdat")), stringsAsFactors=T,
                              key=c("EffectiveDate", "Scenario"))

# --- 3. Mortgage credit dataset
ptm <- proc.time() # for runtime calculations (ignore)
# Import, then recast data into a more pliable data.table object for greater memory efficiency

# - User-dependent dataset extraction (subsampled vs full)
if (Sys.getenv("USERNAME") == "Arno Botha") { # Dr Arno Botha | Kralkatorrik-machine
  dat.raw <- as.data.table(read_sas(paste0(genRawPath, "creditdata_final.sas7bdat")), stringsAsFactors=T) 
} else if (Sys.getenv("USERNAME") == "R8873885") { # Bernard Scheepers
  dat.raw <- as.data.table(read_sas(paste0(genRawPath, "creditdata_final_sub.sas7bdat")), stringsAsFactors=T) 
}
proc.time() - ptm # IGNORE: elapsed runtime

# - Save to disk(zip) for quick disk-based retrieval later
pack.ffdf(paste0(genPath, "creditdata_final1"), dat.raw); gc()


# --- 4. Input fields associated with mortgage credit dataset

ptm <- proc.time()# for runtime calculations (ignore)
# - User-dependent dataset extraction (subsampled vs full)
# NOTE: Import the dataset, but then recast it into a more pliable data.table object for greater memory efficiency
if (Sys.getenv("USERNAME") == "Arno Botha") { # Dr Arno Botha | Kralkatorrik-machine
  datInput.raw <- as.data.table(read_sas(paste0(genRawPath, "creditdata_input.sas7bdat")), stringsAsFactors=T) 
} else if (Sys.getenv("USERNAME") == "R8873885") { # Bernard Scheepers
  datInput.raw <- as.data.table(read_sas(paste0(genRawPath, "creditdata_input_sub.sas7bdat")), stringsAsFactors=T) 
}
proc.time() - ptm # IGNORE: elapsed runtime

# - Save to disk (zip) for quick disk-based retrieval later
pack.ffdf(paste0(genPath, "creditdata_input1"), datInput.raw); gc()
