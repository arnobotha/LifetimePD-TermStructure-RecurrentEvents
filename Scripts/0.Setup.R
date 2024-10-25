# =================================== SETUP =============================================
# Setting up R environment, parameters, and function definitions
# ---------------------------------------------------------------------------------------
# PROJECT TITLE: Default survival modelling
# SCRIPT AUTHOR(S): Dr Arno Botha, Roundel Bester, Marcel Muller, Roland Breedt, Bernard Scheepers

# DESCRIPTION: 
# This script installs and loads various libraries and packages, compiles all
# custom functions, and set requisite parameters.
# ---------------------------------------------------------------------------------------
# -- Inputs:
#   - DelinqM.R | Delinquency measures and related functions 
# =======================================================================================



# ================ 0. Library setup

# ------ Install and load packages
# - data access and big data management
require(haven) # for SAS imports
require(ETLUtils)
require(ffbase)
require(ff)
tempPath <- "C:/TempData"; options("fftempdir"=tempPath)

# for data wrangling
require(tidyr)
require(dplyr)
require(data.table)
require(lubridate)
require(readr)
require(bit64) # for very big numeric values
require(stringr) # common string operations, e.g, str_pad
require(purrr) # mapping functions from tidyverse in working with matrices, lists

# for advanced looping functionality in simulation tasks
require(doBy)
require(foreach)
require(doParallel)

# for analyses
require(Hmisc)
require(moments) # for using skewness() function
require(regclass) # for VIF

# for modelling
require(survival) # for survival modelling
require(zoo)
require(car)
require(survivalROC) # for time-dependent ROC-analysis from Heagerty et al.
#require(survAUC) # for time-dependent ROC-analysis (alternative from Potapov et al.)
#require(tdROC) # for time-dependent ROC-analysis ([outdated?] alternative from Li et al.)
#require(timeROC) # for time-dependent ROC-analysis ([outdated?] alternative from Blanche)
require(pROC); require(ROCR) # both for cross-sectional ROC-analysis (main:pROC)
require(discSurv)
require(MASS)
require(ldatools)

#for plots
require(ggplot2)
require(ggpp) # Extensions to ggplot2, particularly geom_table
require(scales)
require(ggthemes)
require(RColorBrewer)
require(extrafont) #remotes::install_version("Rttf2pt1", version = "1.3.8"); Sys.setenv(R_GSCMD="C:/Program Files/gs/gs9.55.0/bin/gswin32c.exe"); font_import(); loadfonts(); loadfonts(device="win")
require(survminer)
require(gridExtra)
require(corrplot)
require(Metrics)
require(treem)
require(treemapify)



# ================ 1. Parametrisation

# - general R options
options(scipen=999) # Suppress showing scientific notation

# - Parameters used in calculating delinquency measures
sc.Thres <- 0.9; # repayment ratio - g1
d <- 3 # default threshold for g0/g1-measures of delinquency (payments in arrears)
k <- 6 # Probation period

# -- Path variables | User-dependent

if (Sys.getenv("USERNAME") == "Arno Botha") {
  path_cust <- "E:/WorkLife/Analytix/Research/LifetimePD-TermStructure-RecurrentEvents/Scripts/"
  
  # - Common path for storing important R-objects as back-up
  genObjPath <- "E:/WorkLife/Analytix/Research/LifetimePD-TermStructure-RecurrentEvents/Objects/"
  
  # - Common path for saving important analytics (e.g., sampling)
  genFigPath <- "E:/WorkLife/Analytix/Research/LifetimePD-TermStructure-RecurrentEvents/Figures/"
  
  # - Common path for saving big data objects
  genPath <- "E:/DataDump/FNB SLC/LifetimePD-TermStructure-RecurrentEvents_Data/"
  
  # - Common path for importing raw data
  genRawPath <- "E:/DataDump/FNB SLC/"
  
} else if (Sys.getenv("USERNAME") == "R8873885") { # Bernard
  
  # - Common path for saving big data objects
  genPath <- "C:/BMI Data/LifetimePD-TermStructure-RecurrentEvents_Data/"
  
  # - Common path for importing raw data
  genRawPath <- "C:/BMI Data/"
  
  # - Custom path where R-scripts are saved
  path_cust <- "C:/Users/R8873885/OneDrive - FRG/Documents/LifetimePD-TermStructure-RecurrentEvents/Scripts/"
  
  # - Common path for storing important R-objects as back-up
  genObjPath <- "C:/Users/R8873885/OneDrive - FRG/Documents/LifetimePD-TermStructure-RecurrentEvents/Objects/"
  
  # - Common path for saving important analytics (e.g., sampling)
  genFigPath <- "C:/Users/R8873885/OneDrive - FRG/Documents/LifetimePD-TermStructure-RecurrentEvents/Figures/"
  
}



# ================ 2. Custom functions

# ------ Custom function definitions
# - Load all custom functions defined in a separate R-script
source(paste0(path_cust,"0a(i).CustomFunctions.R"))

# - True End procedure functions defined in a separate R-script
source(paste0(path_cust,"0a(ii).TruEnd.R"))

# - Compile Delinquency Calculation Functions (CD, MD/DoD)
source(paste0(path_cust,'0a(iii).DelinqM.R'))

# - Survival functions
source(paste0(path_cust,'0a(iv).SurvFunc.R'))

# - Survival functions - Recurrent events
source(paste0(path_cust,'0a(v).SurvFunc_RecurrentEvents.R'))


