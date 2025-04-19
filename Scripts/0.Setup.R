# =================================== SETUP =============================================
# Setting up R environment, parameters, and function definitions
# ---------------------------------------------------------------------------------------
# PROJECT TITLE: Default survival modelling
# SCRIPT AUTHOR(S): Dr Arno Botha, Roundel Bester, Marcel Muller, Roland Breedt, 
#                   Bernard Scheepers

# DESCRIPTION: 
# This script installs and  loads various libraries and packages, compiles all
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
#require(timeROC) # for time-dependent ROC-analysis from Blanche2013 (disavowed in script 0b(iii)). DO NOT USE IN CREDIT DOMAIN
#require(risksetROC) # for time-dependent ROC-analysis (I/D Cox regression method from Heagerty, P.J., Zheng Y. (2005))
require(pROC); require(ROCR) # both for cross-sectional ROC-analysis (main:pROC)
require(discSurv)
require(MASS)
#require(ldatools) ### AB: ??? why necessary?

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



# ================ 1. Parametrisation

# - general R options
options(scipen=999) # Suppress showing scientific notation

# - Parameters used in calculating delinquency measures
sc.Thres <- 0.9; # repayment ratio - g1
d <- 3 # default threshold for g0/g1-measures of delinquency (payments in arrears)
k <- 6 # Probation period

# -- Path variables | User-dependent

if (Sys.getenv("USERNAME") == "Arno Botha") { # Dr Arno Botha | Kralkatorrik-machine
  # - Common path for saving large R-objects as back-up and/or as reusable checkpoints
  genPath <- "E:/DataDump/RetailMortgages-FNB/LifetimePD-TermStructure-RecurrentEvents_Data/"
  # - Common path from which raw big datasets are imported
  genRawPath <- "E:/DataDump/RetailMortgages-FNB/"
  # - Common path for sourcing R-scripts in main codebase
  path_cust <- "E:/Backupz/Google Drive/WorkLife/Analytix/R&D Codebases/LifetimePD-TermStructure-RecurrentEvents/Scripts/"
  # - Common path for storing important (but small!) R-objects as back-up
  genObjPath <- "E:/Backupz/Google Drive/WorkLife/Analytix/R&D Codebases/LifetimePD-TermStructure-RecurrentEvents/Objects/"
  # - Common path for saving important analytics and figures
  genFigPath <- "E:/Backupz/Google Drive/WorkLife/Analytix/R&D Codebases/LifetimePD-TermStructure-RecurrentEvents/Figures/"
  
} else if (Sys.getenv("USERNAME") == "w8873885") { # Bernard Scheepers
  # - Common path for saving large R-objects as back-up and/or as reusable checkpoints
  genPath <- "C:/Data/LifetimePD-TermStructure-RecurrentEvents_Data/"
  # - Common path from which raw big datasets are imported
  genRawPath <- "C:/Data/"
  # - Common path for sourcing R-scripts in main codebase
  path_cust <- "C:/Users/w8873885/Desktop/Article/LifetimePD-TermStructure-RecurrentEvents/Scripts/"
  # - Common path for storing important (but small!) R-objects as back-up
  genObjPath <- "C:/Users/w8873885/Desktop/Article/LifetimePD-TermStructure-RecurrentEvents/Objects/"
  # - Common path for saving important analytics and figures
  genFigPath <- "C:/Users/w8873885/Desktop/Article/LifetimePD-TermStructure-RecurrentEvents/Figures/"
}



# ================ 2. Custom functions

# ------ Custom function definitions
# - Load all custom functions defined in a separate R-script
source(paste0(path_cust,"0a.CustomFunctions.R"))

# - True End procedure functions defined in a separate R-script
source(paste0(path_cust,"TruEnd.R"))

# - Compile Delinquency Calculation Functions (CD, MD/DoD)
source(paste0(path_cust,'DelinqM.R'))

# - Custom survival-related functions - generic; 
### AB: Needs to be dissected a bit and collapsed into other scripts, but only after closeout!
# I'm not longer sure of its utility if its primary function (Schoenfeld) is implemented in the next/your script
source(paste0(path_cust,'0b(i).FunkySurv.R'))

# - Custom survival-related functions - Residuals (Cox-Snell, Schoenfeld)
source(paste0(path_cust,'0b(ii).FunkySurv_Residuals.R'))

# - Custom survival-related functions - time-dependent ROC-analyses and unit tests
source(paste0(path_cust,'0b(iii).FunkySurv_tROCkit.R'))
