# =============================== FINITE-SAMPLE BEHAVIOUR OF tROC ======================================
# A simulation study of the finite-sample behaviour of the classical NN-based tROC vs the
# clustered tROC-extension. Performance measures include: 1) bias; 2) empirical standard deviation; and
# 3) RMSE.
# ------------------------------------------------------------------------------------------------------
# PROJECT TITLE: Default survival modelling
# SCRIPT AUTHOR(S): Dr Arno Botha (AB)
# ------------------------------------------------------------------------------------------------------
# -- Script dependencies:
#   - 0.Setup.R

# -- Inputs:
#   - Various simulation parameters
#
# -- Outputs:
#   - <Analytics>
# ------------------------------------------------------------------------------------------------------


# ------ 0. Preliminaries  ------ 
# Set various simulation parameters

B <- 500 # number of Monte Carlo replications
vecN <- c(2500, 5000, 10000, 20000) # Number of independent loans per simulated portfolio
sT <- 36 # Number of observations per typical spell
maxNumSpells <- 5 # maximum number of spells that can be generated per loan
t_eval <- 12 # Time horizon(s) at which tROC/tAUC is to be evaluated

# - Strength of within-loan and within-spell dependencies by specifying variance components
# The first parameter (variance component) controls the variance contribution of the 
# loan-level random effect. Larger values imply stronger dependence amongst observations 
# belonging to different spells of the same loan.
vecRhoLoan <- c(0, 0.25, 0.5, 0.75)
# The second parameter (variance component) controls the variance contribution of the 
# spell-level random effect. Larger values imply stronger dependence amongst observations 
# belonging to the same spell of a loan.
vecRhoSpell <- c(0, 0.25, 0.5, 0.75)

# - Parameters related to Cox-regression
# Note: these parameters are found interactively later in this script such that 
# the resulting credit dataset has realistic characteristics.
beta_Cox <- 0.75 # Cox coefficient determining the association between generated marker and event hazard
lambda_0_Cox <- 0.005 # Baseline hazard h_0




# ------ 1a. Data-generating process ------ 
simulateData <- function(nObs, nSpellObs, maxNumSpells, rhoLoan=0.25, rhoSpell=0.25,
                         beta_Cox=0.75, lambda_0_Cox=0.001, seedVal = NA) {
  # --- Testing conditions
  # nObs <- vecN[2]; nSpellObs <- sT; rhoLoan <- vecRhoLoan[2]; rhoSpell <- vecRhoSpell[2]; seedVal <- 1234
  
  # - Set seed
  if (!is.na(seedVal)) set.seed(seedVal, kind="Mersenne-Twister")
  
  
  # --- Generate loan-level random effects U_i from a specified distribution
  # U_i induces dependence amongst different spells belonging to the same loan,
  # including observations from different spells of that loan.
  # rhoLoan controls the variance contribution of this loan-level effect.
  datLoanEffects <- data.table(ObsID = seq_len(nObs),
                               U = rnorm(nObs, mean=0, sd=sqrt(rhoLoan)))
  
  
  # --- Generate recurrent-event spells sequentially
  # A subsequent spell is only generated if the preceding spell experienced the event of interest
  # (e.g., loan default), after which the entity is assumed to have cured and become subject to the
  # main risk event again.
  
  datSpells <- list() # object in which sequentially-generated spells will be stored
  # Initially, all entities are eligible for their first spell
  datActive <- data.table(ObsID = seq_len(nObs))
  
  # - Iterate over maximum number of spells
  for (j in seq_len(maxNumSpells)) {
    # Testing conditions
    # j <- 1
    
    # - Stop the process should there be no more active loans
    if (NROW(datActive) == 0) break;
    
    
    # --- Create Observation-spell data structure, e.g., denoted by (i,j) for 
    # loan i and spell j, assuming 1 observation per such pair 
    dt <- data.table(ObsID = datActive$ObsID, SpellNum = j)
    
    # - Expand this structure such that each loan-spell pair has the designated number of 
    # (starting) observations
    dt <- dt[rep(seq_len(.N), each=nSpellObs)]
    dt[, TimeInSpell := 1:.N, by=.(ObsID, SpellNum)]
    
    # - Attach loan-level random effects U_i
    dt <- merge(dt, datLoanEffects, by="ObsID", all.x=T)
    
    
    # --- Generate spell-level random effects V_{ij}
    # V_{ij} induces dependence amongst different observations belonging to the same spell.
    # rhoSpell controls the variance contribution of this spell-level effect
    
    # First, retrieve the Loan-Spell structure (i,j)
    datSpellEffects <- unique(dt[,.(ObsID, SpellNum)])
    # Then, generate the spell-level effects
    datSpellEffects[, V := rnorm(.N, mean=0, sd=sqrt(rhoSpell))]
    # Lastly, attach these effects to the main data object
    dt <- merge(dt, datSpellEffects, by=c("ObsID", "SpellNum"), all.x=T)
    
    
    # --- Generate observation-level random effects \epsilon_{ijt} from standard normal distribution
    dt[, epsilon := rnorm(.N, mean=0, sd=1)]
    
    
    # --- Calculate overall marker and associated hazard quantities
    # Model is given by M = U_i + V_{ij} + \epsilon_{ijt}
    dt[, M := U + V + epsilon]
    # Calculate hazard using Cox-regression model: h_{ijt} = lambda_0 * exp( \beta * M_{ijt})
    dt[, hazard := lambda_0_Cox * exp(beta_Cox*M)]
    # Calculate cumulative hazard over time t
    dt[, Cumul_hazard := cumsum(hazard * 1), by=.(ObsID,SpellNum)] # assumes an interval length of 1 (discrete-time)
    
    
    # --- Generate Event times 
    # First generate a latent random variable that represents a default intensity 'threshold' 
    # for each loan-spell. Default is not an absorbing state at the loan-level, hence a new
    # threshold is generated independently for each spell.
    # This threshold is generated as E_{ij} ~ Exp(1) by definition.
    # Loan default is then induced whenever the cumulative default hazard breaches this threshold.
    # We generate this threshold from the exponential distribution given the fundamental relationship
    # of survival probability: S(t) = P(T>t) = exp(-H(t)), where H(t) is the cumulative hazard function
    # We shall generate such thresholds for each loan-spell construct (i,j)
    
    # - Generate thresholds
    # First, retrieve the Loan-Spell structure (i,j)
    datSpellThresholds <- unique(dt[,.(ObsID, SpellNum)])
    # Generate thresholds from specified distribution
    datSpellThresholds[, Threshold := rexp(.N, 1)]
    # Attach these thresholds to the main data object
    dt <- merge(dt, datSpellThresholds, by=c("ObsID", "SpellNum"), all.x=T)  
    
    # - Calculate event indicators
    # Default is induced whenever H_{ij}(t) >= E_{ij}
    dt[, Event_Ind := as.integer(Cumul_hazard >= Threshold)]
    
    # - Retain generated history only up to the first instance of default, if it exists; otherwise all history
    dt[, Cumul_Event_Ind := cumsum(Event_Ind), by=.(ObsID, SpellNum)]
    dt <- dt[Cumul_Event_Ind <= 1, ]
    
    # - Store completed spell
    datSpells[[j]] <- dt
    
    
    # --- Determine which loans are eligible for the next spell
    # A subsequent spell can only exist if the current spell has experienced the event of interest
    datActive <- dt[, .(HasEvent = max(Event_Ind, na.rm=T)), by=.(ObsID)][HasEvent==1, .(ObsID)]
    
  }
  
  
  # --- Combine all generated spells and perform miscellaneous tasks
  dt <- rbindlist(datSpells, use.names=T, fill=T)
  
  # - Create miscellaneous fields
  # Calculate spell age
  dt[, SpellAge := .N, by=.(ObsID, SpellNum)]
  # Obtain maximum spell number per loan 
  dt[, MaxSpellCount := max(SpellNum), by=.(ObsID)]
  # Set start column
  dt[, Start := TimeInSpell - 1]
  # Set Obs-spell id (e.g., LoanID-SpellNum)
  dt[, KeyID := paste0(ObsID, "_", SpellNum)]
  # Summarise period-specific event indicators to spell-level quantities
  dt[, Max_Event_Ind := max(Event_Ind, na.rm=T), by=.(ObsID, SpellNum)]
  
  # - Create flag for denoting illegal spell progression, i.e., subsequent spells that 
  # occur without the first spell defaulting
  dt[, IllegalSpellProgression_Ind := ifelse(SpellAge==sT & Max_Event_Ind==0 & SpellNum < MaxSpellCount, 1, 0), 
     by=.(ObsID,SpellNum)]
  
  # - Remove intermediary variables
  dt[, c("Cumul_Event_Ind", "Max_Event_Ind") := NULL]
  
  return(dt)
}



# ------ 1b. Find realistic values for parameters ------  

# -- Simulate credit data
datCreditSim <- simulateData(nObs=vecN[2], nSpellObs=sT, maxNumSpells=maxNumSpells, 
                             rhoLoan=vecRhoLoan[2], rhoSpell=vecRhoSpell[2],
                             beta_Cox=0.75, lambda_0_Cox=0.005, seedVal=1234)

# -- High-level analysis

# max spell number distribution
hist(unique(datCreditSim[, .(ObsID, MaxSpellCount)])$MaxSpellCount, breaks="FD")
table(unique(datCreditSim[, .(ObsID, MaxSpellCount)])$MaxSpellCount) %>% prop.table()
# Should roughly correspond to given probability distribution
#      1      2      3      4      5 
# 0.7718 0.1696 0.0440 0.0106 0.0040 

# Aggregate to loan-spell level and retrieve and plot the distribution of spell ages
hist(unique(datCreditSim[, .(ObsID, SpellNum, SpellAge)])$SpellAge, breaks="FD")

# proportion of spells ending in default
datCreditSim[, .(HasDefault = max(Event_Ind)), by=.(ObsID, SpellNum)][,sum(HasDefault)] / 
  unique(datCreditSim[, .(ObsID,SpellNum)])[,.N] * 100 # 23.5177%

# Diagnostic: Illegal spell progression prevalence
datCreditSim[TimeInSpell==1 & IllegalSpellProgression_Ind == 1,.N] / datCreditSim[TimeInSpell==1,.N] * 100



# ------ 1c. Different tROC-calculations

# -- Calculate the clustered tAUC (compact) by KeyID
ptm <- proc.time() #IGNORE: for computation time calculation;
objROC_0 <- tROC.Clustered_compact(datCreditSim, t_eval = 12, sLambda = 0.05, numDigits = 2, numThreads = 12)
objROC_0$AUC; proc.time() - ptm
### RESULTS: 0.6434665 (when not ordering vStartTimes)
### RESULTS: 0.7050944 (when ordering vStartTimes; now the default implementation)


# -- Calculate the clustered tAUC (verbose) by KeyID
ptm <- proc.time() #IGNORE: for computation time calculation;
objROC <- tROC.multi(datGiven=datCreditSim, Markers=datCreditSim$M, month_End=t_eval, sLambda=0.05, 
                     numDigits=2, eventVal=1, fld_ID="KeyID", fld_Event="Event_Ind", 
                     fld_StartTime="Start", fld_EndTime="TimeInSpell", 
                     graphName="SimStudy_Truth", genFigPath=paste0(genFigPath, "Proof of Concepts/"),
                     caseStudyName=paste0("SimStudy_", t_eval), numThreads=12, logPath=genPath)
objROC$AUC; objROC$ROC_graph
proc.time() - ptm
### RESULTS: 0.6434665 (when not ordering vStartTimes)
### RESULTS: 0.7050944 (when ordering vStartTimes; now the default implementation)


# -- Calculate the classical tAUC (compact)
ptm <- proc.time() #IGNORE: for computation time calculation;
tROC1.classic <- tROC.Classical_compact(datCreditSim, t_eval=t_eval, sLambda = 0.05, 
                                        numDigits = 2, eventVal = 1, numThreads = 12)
tROC1.classic$AUC
proc.time() - ptm
### RESULTS: 0.6421406  (when not ordering vStartTimes; ImplementFix=F)
### RESULTS: 0.7052204 (when ordering vStartTimes; now the default implementation)


# -- Calculate the classical tAUC (verbose)
ptm <- proc.time() #IGNORE: for computation time calculation;
tROC2.classic <- tROC.multi(datGiven=datCreditSim, Markers=datCreditSim$M, month_End=t_eval, sLambda=0.05, 
                     numDigits=2, eventVal=1, fld_Event="Event_Ind", fld_StartTime="Start", fld_EndTime="TimeInSpell", 
                     graphName="SimStudy_Truth", genFigPath=paste0(genFigPath, "Proof of Concepts/"),
                     caseStudyName=paste0("SimStudy_", t_eval), numThreads=12, logPath=genPath)
tROC2.classic$AUC; tROC2.classic$ROC_graph
proc.time() - ptm
### RESULTS: 0.6421406 (when not ordering vStartTimes)
### RESULTS: 0.7052204 (when ordering vStartTimes; now the default implementation)


# -- Calculate the classical tAUC (Heagerty version from survivalROC package)
ptm <- proc.time() #IGNORE: for computation time calculation;
objROC_Heagerty <- survivalROC::survivalROC(Stime=datCreditSim$TimeInSpell, status=datCreditSim$Event_Ind, 
                                            entry=datCreditSim$Start, method = "NNE",  span=0.05, 
                                            predict.time=t_eval, marker=round(datCreditSim$M,2))
objROC_Heagerty$AUC
proc.time() - ptm
### RESULTS: 0.6421406 (vStartTimes is unordered by ranks of end times, which is a latent bug in this function)







# ------ 1c. Generate a "truth" dataset for testing purposes

# --- Simulate the "truth" using simulateData()
# First calculate a baseline truth given some parametrisation, purely for testing purposes
datTruth <- simulateData(nObs=1e05, nSpellObs=sT, maxNumSpells=maxNumSpells, 
                         rhoLoan=0.25, rhoSpell=0.25,
                         beta_Cox=beta_Cox, lambda_0_Cox=lambda_0_Cox, seedVal=1234)





# ------ 2. Calculate the true tROC-quantities on a generated truth ------ 

# - Calculate "reference tAUC-value using very large datTruth object using a bespoke function from script 0b(iii)
ptm <- proc.time() #IGNORE: for computation time calculation;
objROC_0 <- tROC.Clustered_compact(datTruth, t_eval = 12, sLambda = 0.05, numDigits = 2, numThreads = 12)
objROC_0$AUC; proc.time() - ptm
### RESULTS: AUC up to t: 0.7518179 (numDigits=3, lambda_0_Cox=0.01), achieved in 5384.03  secs ( 89.7 mins)
### RESULTS: AUC up to t: 0.751558 (numDigits=2, lambda_0_Cox=0.01), achieved in 612.06  secs ( 10.2 mins)
### RESULTS: AUC up to t: 0.7536163 (numDigits=2, lambda_0_Cox=0.001), achieved in 748.22  secs (12.5  mins)






# ------ 3. Monte Carlo simulation on testing finite-sample behaviour ------ 
# The finite-sample behaviour is tested of the MAA-estimator (mean-adjusted Akritas)

# This testing is based on a Monte Carlo simulation, and includes the bias and 
# variance of \hat{tAUC}_b for the tAUC-value within the bth Monte Carlo run. 
# This includes the following steps:
#   1. For each parametrisation of \rho_loan and \rho_spell, generate the "truth" according to a known
#       data-generating process (DGP) using a very large sample size
#   2. For each "truth", estiamte both the clustered tROC-extension and the classical tROC objects
#   3. For each "truth", simulate a credit dataset according to the same DGP across varying portfolio size, but
#      do so within a Monte Carlo setup, whereby (B) credit datasets are generated
#   4. For each "bth" Monte Carlo run, calculate the tAUC at t_eval using the clustered tROC-extension
#   5. Collect the (B) tAUC-values across the Monte Carlo runs into a single data object
#   6. Calculate the bias and variance for the current parameter scenario, given (B) Monte Carlo runs

# - Reporting
cat("New Job: Estimating the tAUC",
    file=paste0(logPath,"assesslog_SimStudy_MAA_Control.txt"), append=F)

# - Iteration parameter
counter <- 1;
numThreads <- 13

# - Load up partial simulation results file, if it exists

if (file.exists(paste0(genPath,"SimStudy_Results"))) {
  suppressWarnings(rm(datResults_study))
  unpack.ffdf(paste0(genPath,"SimStudy_Results"), tempPath);gc()
} 

# - Monte Carlo setup across a parameter scenario grid of
ptm <- proc.time() #IGNORE: for computation time calculation;
for (k in seq_along(vecRhoLoan)) { 
  for (j in seq_along(vecRhoSpell)) {
    
    datTruth <- NULL # Prepare truth object for current parameter scenario (rhoLoan, rhoSpell)
    
    for (i in seq_along(vecN)) { # Runtime for B=500: 62.95 mins
      
      # - Testing conditions
      # k <- 1; j <- 3; i <- 1
      
      # - Reporting
      cat(paste0("\n - [", Sys.time(), "] Running parameter scenario (rhoLoan=",vecRhoLoan[k],
                 ", rhoSpell=", vecRhoSpell[j],", nLoans=", vecN[i],") .. "), 
          file=paste0(logPath,"assesslog_SimStudy_MAA_Control.txt"), append=T);
      
      # - Determine if current parameter scenario (rhoLoan, rhoSpell) and portfolio size has been run before
      if (exists('datResults_study')) {
        dFiltered <- subset(datResults_study, rhoLoan == vecRhoLoan[k] & rhoSpell == vecRhoSpell[j] & 
                              nLoans==vecN[i])
        if (NROW(dFiltered) > 0) bSkip <- T else bSkip <- F
      } else bSkip <- F
      
      # - Generate truth only once for the current i, but only if current parameter scenario 
      # (rhoLoan, rhoSpell) has not been run before
      if (is.null(datTruth) & bSkip == F) {
        
        # - Reporting
        cat(paste0("\n\t -- 1. [", Sys.time(), "] Generating the 'truth' for parameter scenario (rhoLoan=",vecRhoLoan[k],
                   ", rhoSpell=", vecRhoSpell[j],"), and estimating the reference tAUC-value (clustered) .. "), 
            file=paste0(logPath,"assesslog_SimStudy_MAA_Control.txt"), append=T);
        
        # - Generate the "truth" (very large dataset) using current parameter scenario
        datTruth <- simulateData(nObs=1e5, nSpellObs=sT, maxNumSpells=maxNumSpells, 
                                 rhoLoan=vecRhoLoan[k], rhoSpell=vecRhoSpell[j],
                                 beta_Cox=beta_Cox, lambda_0_Cox=lambda_0_Cox, seedVal=1000000*k + 10000*j)
        
        # - Calculate the large-sample reference tAUC-value using very large datTruth object
        # This uses the MAA-estimator with the NN-estimator for the marker-conditional survival probability S(t|M=m)
        # See script 0b(iii) for the compact function definition
        objROC_0 <- tROC.Clustered_compact(datTruth, t_eval=t_eval, sLambda=0.05, numDigits=2, numThreads=numThreads)
        
        # - Reporting
        cat(paste0("\n\t -- 2. [", Sys.time(), "] Estimating another reference tAUC-value (classical) .. "), 
            file=paste0(logPath,"assesslog_SimStudy_MAA_Control.txt"), append=T);
        
        # - Recalculate the classical tAUC assuming independence amongst observations, i.e., 
        # classical tROC from Heagerty et al. (2000)
        objROC_1 <- tROC.Classical_compact(datTruth, t_eval=t_eval, sLambda = 0.05, 
                                           numDigits = 2, eventVal = 1, numThreads = numThreads)
      }
      
      # - Initiate Monte Carlo setup, but only if the current parameter scenario and portfolio size 
      # does not exist in the result set
      if (bSkip == F) {
        # - Preparations
        datResults <- NULL
        
        # - Reporting
        cat(paste0("\n\t -- 3. [", Sys.time(), "] Running parameter scenario (rhoLoan=",vecRhoLoan[k],
                   ", rhoSpell=", vecRhoSpell[j], ", nLoans=", vecN[i], ") across MC-iterations .. "), 
            file=paste0(logPath,"assesslog_SimStudy_MAA_Control.txt"), append=T);
        cat(paste0("Running MC-iterations for parametrisation (rhoLoan=",vecRhoLoan[k],
                   ", rhoSpell=", vecRhoSpell[j], ", nLoans=", vecN[i], ") .. "), 
            file=paste0(logPath,"assesslog_SimStudy_MAA_MonteCarlo.txt"), append=F);
        
        # - Implement credit data generation within a Monte Carlo setup across B replications
        for (b in seq_len(B)) {
          # - Testing conditions
          # k <- 2; j <- 2; i <- 1; b <- 1
          
          # - Reporting
          cat(paste0("\n - Running MC-iteration [", b, " of ", B, "] for parametrisation (rhoLoan=",vecRhoLoan[k],
                     ", rhoSpell=", vecRhoSpell[j], ", nLoans=", vecN[i], ") .. "), 
              file=paste0(logPath,"assesslog_SimStudy_MAA_MonteCarlo.txt"), append=T);
          
          # - Generate dataset according to the same data-generating process (DGP) as used in the "truth"
          datCreditSim <- simulateData(nObs=vecN[i], nSpellObs=sT, maxNumSpells=maxNumSpells, 
                                       rhoLoan=vecRhoLoan[k], rhoSpell=vecRhoSpell[j],
                                       beta_Cox=beta_Cox, lambda_0_Cox=lambda_0_Cox, 
                                       seedVal=1000000*k + 10000*j + 100*i + b) 
          
          # - Calculate tAUC using the clustered tROC-extension
          objROC <- tROC.multi(datGiven=datCreditSim, Markers=datCreditSim$M, month_End=t_eval, sLambda=0.05, 
                               numDigits=2, fld_ID="KeyID", fld_Event="Event_Ind", 
                               eventVal=1, fld_StartTime="Start", fld_EndTime="TimeInSpell", 
                               Graph = F, caseStudyName=paste0("SimStudy_", t_eval), 
                               numThreads=numThreads, logPath=genPath)
          
          # - Calculate the classical tAUC assuming independence amongst observations, i.e., 
          # the classical tROC-estimator from Heagerty et al. (2000)
          objROC_classical <-  tROC.Classical_compact(datCreditSim, t_eval=t_eval, sLambda = 0.05, 
                                                      numDigits = 2, eventVal = 1, numThreads = numThreads)
          
          # - Collect results into data object
          datResults_indiv <- data.table(rhoLoan=vecRhoLoan[k], rhoSpell=vecRhoSpell[j], nLoans=vecN[i],
                                         B=b, tAUC=objROC$AUC, tAUC_Classical = objROC_classical$AUC)
          if (b == 1) datResults <- datResults_indiv else datResults <- rbind(datResults, datResults_indiv)
        }
        
        # - Bias and mean and variance calculation
        # - For the current parameter scenario, calculate the bias and mean and variance across (B) number of 
        # Monte Carlo iterations. Specifically, the bias is calculated as 
        # \hat{Bias} = 1/B \sum_b (tAUC(b) - tAUC_truth), where tAUC_truth comes from objROC_0 or objROC_1.
        # The sample variance of tAUC is simply calculated across the (B) tAUC-values
        sBias <- mean(datResults$tAUC, na.rm=T) - objROC_0$AUC
        sBias_1 <- mean(datResults$tAUC_Classical, na.rm=T) - objROC_1$AUC
        sMean <- mean(datResults$tAUC, na.rm=T)
        sMean_1 <- mean(datResults$tAUC_Classical, na.rm=T)
        sVar <- var(datResults$tAUC, na.rm=T)
        sVar_1 <- var(datResults$tAUC_Classical, na.rm=T)
        
        # - Calculate mean paired difference between classical and clustered tROC-estimators.
        Mean_Diff = mean(datResults$tAUC - datResults$tAUC_Classical, na.rm=T)
        Variance_Diff = var(datResults$tAUC - datResults$tAUC_Classical, na.rm=T)
        
        # - Collect results into data object
        datResults_indiv2 <- data.table(rhoLoan=vecRhoLoan[k], rhoSpell=vecRhoSpell[j], nLoans=vecN[i],
                                        tAUC_0=objROC_0$AUC, tAUC_classical_0=objROC_1$AUC,
                                        tAUC_Mean=sMean, tAUC_classical_Mean=sMean_1,
                                        tAUC_Bias=sBias, tAUC_classical_Bias=sBias_1, 
                                        tAUC_Variance=sVar, tAUC_classical_Variance=sVar_1, 
                                        tAUC_meanDiff=Mean_Diff, tAUC_varDiff=Variance_Diff)
        if (counter == 1) datResults_study <- datResults_indiv2 else
          datResults_study <- rbind(datResults_study,datResults_indiv2)
        
        # - commit results so far to disk
        pack.ffdf(paste0(genPath,"SimStudy_datasets/SimStudy_Results_rhoLoan",vecRhoLoan[k],
                         "_rhoSpell", vecRhoSpell[j],"_nLoans", vecN[i],".zip"), datResults)
        pack.ffdf(paste0(genPath,"SimStudy_Results"), datResults_study)
      }
      
      # - Simply increment the overall loop counter
      counter <- counter + 1 
    }
  }
}
proc.time() - ptm





# ------ 4. Graphing results from the simulation study ------ 


# --- Graph: Clustered tAUC's Bias = E(tAUC) - tAUC_ref ---

# Graphing data preparation
datResults_study[, rhoLoan_fact := factor(sprintf("%.2f",rhoLoan))]
datResults_study[, rhoSpell_fact := factor(sprintf("%.2f",rhoSpell))]

alpha <- 0.05
# Create 95% confidence intervals using z-scores
datResults_study[, ErrMargin := sqrt(tAUC_Variance/500) * qnorm(1-alpha / 2)] # sample size at each point is 500
datResults_study[, tAUC_Bias_lower := tAUC_Bias - ErrMargin]
datResults_study[, tAUC_Bias_upper := tAUC_Bias + ErrMargin]

# Aesthetic engineering
chosenFont <- "Cambria"
vCol <- brewer.pal(n=8, name="Set1")
vName <- bquote("Spell-level VC "*italic(rho)["s"])
vFacets <- as_labeller(c("0.00"="'Loan-level VC '*italic(rho)['l']*' = 0'", 
                         "0.25"="'Loan-level VC '*italic(rho)['l']*' = 0.25'",
                         "0.50"="'Loan-level VC '*italic(rho)['l']*' = 0.50'",
                         "0.75"="'Loan-level VC '*italic(rho)['l']*' = 0.75'"),
                       default=label_parsed)

# - Main plot
(ggObj <- ggplot(datResults_study, aes(x=nLoans, y=tAUC_Bias)) + theme_minimal() + 
    labs(y=bquote("Bias of clustered tAUC at "*italic(t)*" = "*.(t_eval)), x=bquote("Number of loans "*italic(N[p]))) + 
    theme(text=element_text(family=chosenFont),legend.position = "bottom",
          strip.background=element_rect(fill="snow2", colour="snow2"),
          strip.text=element_text(size=8, colour="gray50"), strip.text.y.right=element_text(angle=90)) + 
    geom_line(aes(x=nLoans, y=tAUC_Bias, colour=rhoSpell_fact, linetype=rhoSpell_fact)) + 
    geom_point(aes(x=nLoans, y=tAUC_Bias, colour=rhoSpell_fact, shape=rhoSpell_fact)) + 
    geom_ribbon(aes(x=nLoans, ymin=tAUC_Bias_lower, ymax=tAUC_Bias_upper, fill=rhoSpell_fact), alpha=0.15) + 
    facet_wrap(.~rhoLoan_fact, labeller = labeller(rhoLoan_fact=vFacets), ncol=2) + 
    scale_color_manual(name=vName, values=vCol) + 
    scale_fill_manual(name=vName, values=vCol) + 
    scale_linetype_discrete(name=vName) + 
    scale_shape_discrete(name=vName) + 
    scale_x_continuous(labels=comma) + 
    scale_y_continuous(labels=comma))

# - Save plot
dpi <- 200
ggsave(ggObj, file=paste0(paste0(genFigPath, "SimStudy_Bias.png")), 
       width=1600/dpi, height=1200/dpi, dpi=dpi, bg="white")



# --- Graph: Mean paired difference = E(tAUC_clustered) - E(tAUC_classical) ---

alpha <- 0.05
# Create 95% confidence intervals using z-scores
datResults_study[, ErrMargin := sqrt(tAUC_varDiff/500) * qnorm(1-alpha / 2)] # sample size at each point is 500
datResults_study[, tAUC_meanDiff_lower := tAUC_meanDiff - ErrMargin]
datResults_study[, tAUC_meanDiff_upper := tAUC_meanDiff + ErrMargin]

# - Main plot
(ggObj2 <- ggplot(datResults_study, aes(x=nLoans, y=tAUC_meanDiff)) + theme_minimal() + 
    labs(y=bquote("Mean paired difference between clustered and classical tAUC at "*italic(t)*" = "*.(t_eval)), x=bquote("Number of loans "*italic(N[p]))) + 
    theme(text=element_text(family=chosenFont),legend.position = "bottom",
          strip.background=element_rect(fill="snow2", colour="snow2"),
          strip.text=element_text(size=8, colour="gray50"), strip.text.y.right=element_text(angle=90)) + 
    geom_line(aes(x=nLoans, y=tAUC_meanDiff, colour=rhoSpell_fact, linetype=rhoSpell_fact)) + 
    geom_point(aes(x=nLoans, y=tAUC_meanDiff, colour=rhoSpell_fact, shape=rhoSpell_fact)) + 
    geom_ribbon(aes(x=nLoans, ymin=tAUC_meanDiff_lower, ymax=tAUC_meanDiff_upper, fill=rhoSpell_fact), alpha=0.15) + 
    facet_wrap(.~rhoLoan_fact, labeller = labeller(rhoLoan_fact=vFacets), ncol=2) + 
    scale_color_manual(name=vName, values=vCol) + 
    scale_fill_manual(name=vName, values=vCol) + 
    scale_linetype_discrete(name=vName) + 
    scale_shape_discrete(name=vName) + 
    scale_x_continuous(labels=comma) + 
    scale_y_continuous(labels=comma))

# - Save plot
dpi <- 200
ggsave(ggObj2, file=paste0(paste0(genFigPath, "SimStudy_MeanPairedDiff.png")), 
       width=1600/dpi, height=1200/dpi, dpi=dpi, bg="white")