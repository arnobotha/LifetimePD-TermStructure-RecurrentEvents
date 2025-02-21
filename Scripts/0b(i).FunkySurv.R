# ============================== SURVIVAL FUNCTIONS ==============================
# Defining bespoke customs relating to various generic aspects of survival analysis
# --------------------------------------------------------------------------------
# PROJECT TITLE: Default Survival Modelling
# SCRIPT AUTHOR(S): Marcel Muller, Dr Arno Botha, Bernard Scheepers

# VERSION: 1.0 (July-2023)
# DESCRIPTION: 
# This script defines various functions specific to survival modelling
# that are used elsewhere in this project or, indeed, used across other projects.
# Functions are grouped thematically.
# ================================================================================




# ----------------- 1. Functions related to the modelling process------------------

# --- Function to ensure that the column is indeed part of the given dataset.
# Input:  [cols], Vector containing the column names to be checked if they exist in the dataset.
#         [dataset], Dataset to check on if the columns exist.
# Output: print on whether the columns are in the dataset or which columns are not in the dataset.
colCheck <- function(cols, dataset) {
  # Check if all columns exist in the dataset
  missing_cols <- cols[!(cols %in% colnames(dataset))]
  
  if (length(missing_cols) == 0) { # All columns are in the dataset.
    paste0("SAFE: All columns are in the dataset")
  } else { # List columns not in the dataset.
    warning("WARNING: Some columns are not in the dataset.")
    paste0("Columns not in the dataset: ", paste(missing_cols, collapse = ", "))
  }
}



# --- Function to add and/or remove certain values from a vector
# Input:  [mat], (2 x n) Matrix on which changes will be made. The first column contains the variable names (vars) and the second their respective variable type (vartypes)
#         [Remove], Vector containing the entries to be removed.
# Output:        [m], Updated matrix.
vecChange <- function(mat,Remove=FALSE,Add=FALSE){
  m <- mat
  if (all(Remove != FALSE)) {m <- mat[!(mat$vars %in% Remove)]} # Remove values in [Remove]
  if (all(Add != FALSE)) {m <- rbind(m,Add[!(Add$vars %in% mat$vars)])} # Add values in [Add]
  return(m)
}

# --- Function to detect significant correlations (abs(cor) > 0.6) between vectors.
# Input:  [data], Dataset containing the variables in varlist to dertermine correlations.
#         [varlist], list of variable names to determine correlations.
#         [corrThresh], the absolute correlation threshold above which correlations are deemed significant.
#         [method], the method by which correlations are determined.
#         [Remove], Vector containing the entries to be removed.
# Output: <graph>  Upper half of correlation matrix.
#         <print> Text indicating the variable pairs with high correlation.
corrAnalysis <- function(data, varlist, corrThresh = 0.6, method = 'spearman') {
  # Compute the correlation matrix
  corrMat <- as.data.table(data) %>% subset(select = varlist) %>% cor(method = method)
  
  # Visualize the correlation matrix
  corrplot(corrMat, type = 'upper', addCoef.col = 'black', tl.col = 'black', diag=FALSE,
           tl.srt = 45)
  
  # Find correlation coordinates exceeding the threshold
  corrCoordinates <- which(abs(corrMat) > corrThresh & abs(corrMat) < 1 & upper.tri(corrMat), arr.ind = TRUE)
  
  if(nrow(corrCoordinates) != 0){
    # Create a data table with correlation pairs
    corrProbs <- data.table(x = rownames(corrMat)[corrCoordinates[, 1]], y = colnames(corrMat)[corrCoordinates[, 2]])
    
    # Print the identified correlations
    for (i in 1:nrow(corrProbs)) {
      cat("Absolute correlations of ",percent(corrMat[corrProbs[i, x], corrProbs[i, y]]),
          " found for ", corrProbs[i, x], " and ", corrProbs[i, y],"\n")
    }
  }else{
    cat("No significant correlations were detected")
  }
}



# --- Function to return the appropriate formula object based on the time definition.
#         [TimeDef], Time definition on which the cox ph models are based on.
#         [var], Single variable name
TimeDef_Form <- function(TimeDef="TFD", vars){
  # Create formula based on time definition of the dataset.
  if(TimeDef=="TFD"){# Formula for time to first default time defintion (containing only the fist performance spell).
    formula <- as.formula(paste0("Surv(Start,End,Default_Ind) ~ ",
                                 paste(vars,collapse=" + ")))
  }else if(TimeDef=="PWP_ST"){# Formula for Prentice-Williams-Peterson Spell time definition (containing only the fist performance spell).
    formula <- as.formula(paste0("Surv(Start,End,Default_Ind) ~ PerfSpell_Grp + ",
                                 paste(vars,collapse=" + ")))
  }
}



# --- Function to fit a given formula within a Cox regression model towards extracting Harrell's C-statistic and related quantities
calc_HarrellC <- function(formula, data_train, data_valid, variable="", it=NA, logPath="") {
  model <- coxph(formula,id=get(fldSpellID), data = data_train) # Fit Cox model
  if (!is.na(it)) {
    cat(paste0("\n\t ", it,") Single-factor survival model built. "),
        file=paste0(logPath,"HarrelsC_log.txt"), append=T)
  }
  c <- concordance(model, newdata=data_valid) # Calculate concordance of the model based on the validation set.
  conc <- as.numeric(c[1])# Extract concordance
  sd <- sqrt(c$var)# Extract concordance variability
  lr_stat <- round(2 * (model$loglik[2] - model$loglik[1]),0)# Extract LRT from the model's log-likelihood
  # Return results as a data.table
  return(data.table(Variable = variable, Concordance = conc, SD = sd, LR_Statistic = lr_stat))
}



# --- Function to extract the concordances (Harrell's C) from single-factor models
# Input:  [data_train], Training dataset on which the models are built on.
#         [data_valid], Validation dataset on which the models' concordance is validated on.
#         [variables], Vector containing a list of the variables for the models.
#         [TimeDef], Time definition on which the cox ph models are based on.
# Output: [matResults]  Table containing the concordance, se(concordance) and log ratio of the singular models.
concTable <- function(data_train, data_valid, variables, fldSpellID="PerfSpell_Key",
                      TimeDef="TFD", numThreads=6, genPath) {
  # - Testing conditions
  # data_valid <- datCredit_valid_TFD; TimeDef="TFD"' numThreads=6
  
  # - Iterate across loan space using a multi-threaded setup
  ptm <- proc.time() #IGNORE: for computation time calculation
  cl.port <- makeCluster(round(numThreads)); registerDoParallel(cl.port) # multi-threading setup
  cat("New Job: Estimating B-statistic (1-KS) for each variable as a single-factor survival model ..",
      file=paste0(genPath,"HarrelsC_log.txt"), append=F)
  
  results <- foreach(j=1:length(variables), .combine='rbind', .verbose=F, .inorder=T,
                                 .packages=c('data.table', 'survival'), .export=c('calc_HarrellC')) %dopar%
    { # ----------------- Start of Inner Loop -----------------
      # - Testing conditions
      # j <- 1
      calc_HarrellC(formula=TimeDef_Form(TimeDef,variables[j]), variable=variables[j],
                    data_train=data_train, data_valid=data_valid, it=j, logPath=genPath)
    } # ----------------- End of Inner Loop -----------------
  stopCluster(cl.port); proc.time() - ptm  
  
  # Sort by concordance in descending order.
  setorder(results, -Concordance)
  
  # Return resulting table.
  return(results)
}


# --- Function to calcualte the complement of the KS test statistic "B-statistic"
# Inputs: [formula]: Cox regression formula object; [data_train]: training data
#         [fldSpellID]: Field name of spell-level ID; [vEvents]: spell-level vector of event indicators
#         [seedVal]: Seed value for random number generation; [it]: optional iteration parameter for logging purposes;
#         [logPath]: Optional path for log file for logging purposes
# Outputs: b-statistic (single value)
calcBStat <- function(formula, data_train, fldSpellID="PerfSpell_Key", vEvents, seedVal, it=NA, logPath=NA) {
  # Fit Model
  model <- coxph(formula, data = data_train, id=get(fldSpellID))
  
  if (!is.na(it)) {
    cat(paste0("\n\t ", it,") Single-factor survival model built. "),
        file=paste0(logPath,"BStat_log.txt"), append=T)
  }
  
  # Calculate Cox-Snell (adjusted) residuals
  vCS <- calc_CoxSnell_Adj(model, vIDs=data_train[[fldSpellID]], vEvents=vEvents)
  # Initialize a unit exponential distribution
  set.seed(seedVal, kind = "Mersenne-Twister")
  vExp <- rexp(length(vCS),1)
  # Perform the two-sample Kolmogorov-Smirnov test of distribution equality
  #   H_0: vCS and vExp originates from the same distribution
  #   NOTE: We only desire the KS test statistic in measuring distributional dissimilarity
  #   Then, we subtract this from 1 in creating a coherent statistic; greater is better
  bStat <- 1 - round(suppressWarnings(ks.test(vCS,vExp))$statistic,4)
  return(bStat)
}


# --- Function to extract the B-statistic from a range of models built on a list of variables based on a time definition.
# Input:  [data_train], Training dataset on which the models are built on.
#         [seedVal], Seed value to ensure results are reproducible.
#         [numIt], Number of simulations to calculate the B-statistic.
#         [TimeDef], Time definition on which the cox ph models are based on.
# Output: [Results]  Table containing the concordance, se(concordance) and log ratio of the singular models.
csTable <- function(data_train, variables, TimeDef="TFD", seedVal=1, numIt=5, 
                    fldSpellID="PerfSpell_Key", fldLstRowInd="PerfSpell_Exit_Ind", fldEventInd="Default_Ind",
                    numThreads=6, genPath=NA){
  
  # - Testing conditions
  # data_train <- datCredit_train_TFD; variables<-vars2; TimeDef<-"TFD"; seedVal<-1; numIt<-5; 
  # fldLstRowInd="PerfSpell_Exit_Ind";  fldSpellID="PerfSpell_Key"; fldEventInd="Default_Ind"; numThreads=6
  
  # - Initialize results
  results <- data.frame(Variable = variables, B_Statistic = NA_real_)
  
  # - Data preparation
  # Subset last row per performing spell for Goodness-of-Fit (GoF) purposes
  datLstRow <- copy(data_train[get(fldLstRowInd)==1,])
  vLstRow_Events <- datLstRow[, get(fldEventInd)]
  
  # - Simulate null distribution if seedVal is not NA
  if (!is.na(seedVal)) {
    
    
    # - Iterate across loan space using a multi-threaded setup
    ptm <- proc.time() #IGNORE: for computation time calculation
    cl.port <- makeCluster(round(numThreads)); registerDoParallel(cl.port) # multi-threading setup
    cat("New Job: Estimating B-statistic (1-KS) for each variable as a single-factor survival model ..",
        file=paste0(genPath,"BStat_log.txt"), append=F)
    
    results$B_Statistic <- foreach(j=1:length(variables), .combine='rbind', .verbose=F, .inorder=T,
                      .packages=c('data.table', 'survival'), .export=c('calc_CoxSnell_Adj', 'calcBStat', 'TimeDef_Form')) %dopar%
      
      { # ----------------- Start of Inner Loop -----------------
        # - Testing conditions
        # var <- variables[1]
        calcBStat(formula=TimeDef_Form(TimeDef,variables[j]), data_train=data_train, fldSpellID=fldSpellID, vEvents=vLstRow_Events,
                  seedVal=seedVal, it=j, logPath=genPath)
        
      } # ----------------- End of Inner Loop -----------------
    stopCluster(cl.port); proc.time() - ptm
    
    # Sort results by B statistic in descending order
    results <- results[order(-results$B_Statistic, na.last = TRUE), ]
    
    # Return results and range of B statistics
    return(list(Results = results, Range = diff(range(results$B_Statistic, na.rm = TRUE))))
    
  } else {
    # Perform iterative B calculation when seedVal is NA
    # Initialize Results matrix to contain the number of interations
    matResults <- matrix(NA, nrow = length(variables), ncol = numIt,
                         dimnames = list(variables,
                                         paste0("Iteration_", 1:numIt))) %>%
                          as.data.table()
    
    # - Iterate across loan space using a multi-threaded setup
    ptm <- proc.time() #IGNORE: for computation time calculation
    cl.port <- makeCluster(round(numThreads)); registerDoParallel(cl.port) # multi-threading setup
    cat("New Job: Estimating B-statistics (1-KS) ..",
        file=paste0(genPath,"BStat_log.txt"), append=F)
    
    for (it in seq_len(numIt)) {
      
      cat(paste0("\n Estimating B-statistic (1-KS) for each variable as a single-factor survival model for iteration ", it, " .."),
          file=paste0(genPath,"BStat_log.txt"), append=T)
      
      matResults[, it] <- foreach(j=1:length(variables), .combine='rbind', .verbose=F, .inorder=T,
                                     .packages=c('data.table', 'survival'), .export=c('calc_CoxSnell_Adj', 'calcBStat', 'TimeDef_Form')) %dopar%
        
        { # ----------------- Start of Inner Loop -----------------
          # - Testing conditions
          # var <- variables[1]
          calcBStat(formula=TimeDef_Form(TimeDef,variables[j]), data_train=data_train, fldSpellID=fldSpellID, vEvents=vLstRow_Events,
                    seedVal=seedVal*it, it=j, logPath=genPath)
          
        } # ----------------- End of Inner Loop -----------------
    }
    stopCluster(cl.port); proc.time() - ptm
    
    # Compute additional statistics for the results matrix
    colRanges <- matResults[, lapply(.SD, function(x) diff(range(x, na.rm = TRUE)))] # Calculate the Range of B-statistic value for each iteration
    matResults[, Average := rowMeans(.SD, na.rm = TRUE), .SDcols = patterns("^Iteration_")]# Calculate the average B-statistic for each variable
    matResults[, Variable := variables] 
    matResults <- matResults %>% relocate(Variable, .before=Iteration_1)
    
    #matResults <- cbind(Variables = c(variables,"Range"),matResults)# Add a column to the Results matrix to cross reference the variables with their respective B-statistics
    setorder(matResults,-Average)# Arrange matrix according to average
    
    # Return matrix of B statistics
    return(list(Results=matResults, IterationRanges=colRanges))
  }
}



# HW BS: Comment better the spline function
# -- Cubic spline function
spline_estimation <- function(times, hazard, nknots, degree) {
  n <- length(times)  # Number of observations
  
  # Calculate quantiles for knot placement
  qtiles <- (1:nknots + 1) / (nknots + 2)  # Generate nknots evenly spaced quantiles
  knots <- quantile(times, probs = qtiles, na.rm = TRUE)
  
  # Construct the T matrix efficiently
  matT <- cbind(
    sapply(0:degree, function(j) times^j),                           # Polynomial terms
    sapply(knots, function(k) pmax(0, (times - k)^degree))           # Truncated power basis
  )
  
  # Ordinary Least Squares: Compute the parameter estimates
  coef <- solvet(crossprod(matT), tol=1e-15) %*% crossprod(matT, hazard)  # More efficient than t(T) %*% T and t(T) %*% y
  
  # Polynomial terms: ∑_{j=0}^{d} α_j * t^j
  poly_terms <- sapply(0:degree, function(j) coef[j + 1] * times^j)
  
  # Truncated power basis terms: ∑_{p=1}^{r} α_{p+d} * (t - K_p)^d_+
  truncated_terms <- sapply(1:nknots, function(p) {
    coef[1 + degree + p] * pmax(0, (times - knots[p])^degree)
  })
  
  # Combine polynomial and truncated power terms
  y <- rowSums(poly_terms) + rowSums(truncated_terms)
  y <- pmax(y,0)
  return(y)
}







### AB: Given the work of Bernard, I'm no longer sure of the utility of the below.

# --- function to compute the (unscaled) Schoenfeld residuals for a Cox PH model
#   1) Schoenfeld residuals are computed for each specified variable in the training dataset
#   2) Tests are conducted for significance of the correlation of the Schoenfeld residuals against time
#   3) Graphs are created of the Schoenfeld residuals against time
# Input:  cph - Cox PH model to be assessed
#         dat_train - Dataset used to train the Cox PH model
#         var - Name of the variables for which the Schoenfeld residuals are to be computed
#         id - Name of the column uniquely identifying each row of the dat_train
#         time - Name of the column identifying the associated time in dat_train
#         status - Name of the column identifying the target variable in dat_train
#         verbose - Indicator variable for supressing graphs created by the function
#         max_time - The maximum time for which the graph should be plotted
# Output: A list containing the following:
#         data - A dataset containing the id, time, raw variable's value, and the associated Schoenfeld residual
#         CorTest - A list of the results from test(s) for significance of the correlation of the Schoenfeld residuals against time
#         plots - A list of graphs of the Schoenfeld residuals against time
cph_schoen <- function(cph, var=NULL, dat_train, id, time, status, verbose=T, max_time=NULL){
  # UNIT TEST (VARIABLE INITIALISATION)
  # cph <- cph_Default_PH_test; var <- c("Principal_wins"); id <- "PerfSpell_Key"; time <- "TimeInPerfSpell"; status <- "DefaultStatus1"; max_time <- 240; verbose <- F
  
  # Copying the training dataset to ensure no contamination and ensuring that the dataset is of the correct class (this step usually takes a considerable amount of time for large datasets)
  if (any(class(dat_train) %in% "tbl_df")){
    dat_train2 <- as.data.table(dat_train)
  } else {
    dat_train2 <- copy(dat_train)
  }
  if (any(class(dat_train2) %in% "grouped_df")){
    dat_train2 <- ungroup(dat_train2)
  }
  
  cph_sum <- summary(cph) # Getting a summary of the model
  row_names <- rownames(cph_sum$coefficients) # Getting the names of all the variables in the Cox model (includes levels of categorical variables)
  
  # Object for graphs
  gplots <- list()
  
  # Object for correlation test
  cor_test <- list()
  
  # Getting all the variable names in the model if no names are specified
  if (is.null(var)){
    var <- unlist(strsplit(toString(cph_sum$call$formula[[3]]), '[,+ ]+'))
    var <- var[var!=""]
  }
  
  # Initialising the dataset to be returned
  dat_return <- data.table(ID = numeric(),
                           Time = numeric(),
                           Var_Val = numeric(),
                           Sch_Res = numeric(),
                           Var_Name = as.character(),
                           Var_Name_Base = as.character())
  
  # Initialising the temporary dataset used in the loop for all variables
  col_names <- c(id, time, status, var) # Getting the names of the columns from the training dataset to be subsetted
  dat_temp <- dat_train2[, ..col_names] # Subsetting from the main dataset
  colnames(dat_temp)[1:3] <- c("ID", "Time", "Status") # Renaming the columns for conveinience
  dat_temp <- cbind(dat_temp, predict(cph, dat_train2, type="risk")); colnames(dat_temp)[length(colnames(dat_temp))] <- "RiskScore" # Getting the risk score of each observation
  dat_temp <- merge(dat_temp, dat_temp[, list(SumScore=sum(RiskScore)), by=list(Time)], by="Time", all.x=T) # Getting the total risk score at each time point and merging it back into the dataset
  
  # Computing the Schoenfeld residuals for each selected variable and appending the dataset which to return
  k <- 0 # Counting variable for lists in for loop below
  for (i in 1:length(var)){
    # i <- 1
    var_type <- class(dat_temp[[var[i]]]) # The class of the selected variable
    var_name <- var[i]
    
    col_names2 <- c("ID", "Time", "Status", "RiskScore", "SumScore", var_name)
    dat_temp2 <- dat_temp[, ..col_names2] # Subsetting from the temporary dataset to only include the i'th variable's values
    colnames(dat_temp2)[6] <- "Var_Val"
    dat_temp2[, Var_Name := var_name]; dat_temp2[, Var_Name_Base := var_name] # [Var_Name] is the name of the variable as in the training dataset; [Var_Name_Base] is the level of the categorical variable (equal to Var_Name for numeric variables)
    
    # Computing the Schoenfeld residuals for a numeric variable
    if (var_type == "numeric"){
      # Computing the risk-weighted value of the chosen variable
      dat_temp2[, RW_Val := Var_Val*RiskScore/SumScore]
      # Computing the expected value of the numerical variable at each time point and merging it back into the dataset
      dat_temp2 <- merge(dat_temp2, dat_temp2[, list(Exp_Val=sum(RW_Val)), by=list(Time)], by="Time", all.x = TRUE)
      
      # Subsetting to only include defaulted accounts || Change naming accordingly as structure of data changes
      dat_temp3 <- dat_temp2[Status==1,]; rm(dat_temp2)
      # Computing the Schoenfeld residuals
      dat_temp3[, Sch_Res := Var_Val - Exp_Val] # -447886.85507
      
      # Appending the temporary dataset to the main dataset that is to be returned
      dat_return <- rbind(dat_return, dat_temp3[, list(ID, Time, Var_Val, Sch_Res, Var_Name, Var_Name_Base)])
      
      # Updating the counter variable
      k <- k + 1
      
      # Plotting if verbose = True
      if (verbose==F){
        if (is.null(max_time)){
          max_time <- max(dat_return$Time)
        } # if

        g <- ggplot(dat_temp3, aes(x=Time, y=Sch_Res)) +
             theme_minimal() + xlab(expression(Default~Time~italic(tau[d]))) + ylab(expression(italic(hat(Sc)(x[k~i~tau[d]])))) +
             theme(text=element_text(family=chosenFont),legend.position = "bottom",
                   strip.background=element_rect(colour="grey", fill="#D3D3D3")) +
             facet_wrap(~ Var_Name, strip.position = "right") +
             geom_point(shape=1, col="#3F702F") + geom_smooth(aes(linetype="solid"), col="#043927", fill="#043927") +
             scale_x_continuous(limits = c(NA, max_time)) +
             scale_linetype_manual(name = "Loess-smoother", values = "solid", labels = NULL)
        
        g_name <- paste0("SchRes_",var[i])
        
        gplots[[g_name]] <- g; names(gplots[k])  
      } # if

      # Correlation Test (Formal)
      cor_test[[k]] <- cor.test(dat_temp3$Time, dat_temp3$Sch_Res, method = "spearman"); names(cor_test)[k] <- var[i] # Spearman rank correlation is used for robustness
      
      # Clean up
      rm(dat_temp3)
      
    } else if (var_type %in% c("character", "factor")){
      levels_n <- length(grep(var[i], row_names)) # Computing the number of levels (-1) of the categorical variable | k-1 levels of the categorical variable
      levels <-  substring(row_names[grep(var[i], row_names)], nchar(var[i])+1) # Getting the levels of the categorical variable (-1)
      
      # Computing the Schoenfeld residuals for each level of the categorical variable
      for (j in 1:levels_n){
        # j <- 1
        dat_temp3 <- copy(dat_temp2) # Copying the temporary dataset to increase efficiency with multiple levels of a categorical variable
        dat_temp3[, Var_Name := levels[j]] # Setting [Var_Name] to level j of the categorical variable
        
        # Computing the risk-weighted value of the chosen variable
        dat_temp3[, RW_Val := as.numeric(Var_Val==levels[j])*RiskScore/SumScore]
        # Computing the expected value of the chosen variable at each time point and merging it back into the dataset
        dat_temp3 <- merge(dat_temp3, dat_temp3[, list(Exp_Val=sum(RW_Val)), by=list(Time)], by="Time", all.x = TRUE)
        # Subsetting to only include defaulted accounts
        dat_temp4 <- dat_temp3[Status==1,]; rm(dat_temp3)
        # Computing the Schoenfeld residuals
        dat_temp4[, Sch_Res := as.numeric(Var_Val==levels[j]) - Exp_Val]
        
        # Appending the temporary dataset to the main dataset
        dat_return <- rbind(dat_return, dat_temp4[, list(ID, Time, Var_Val, Sch_Res, Var_Name, Var_Name_Base)])
        
        # Updating the counter variable
        k <- k+1
        
        # Plotting if verbose = True
        if (verbose==F){
          if (is.null(max_time)){
            max_time <- max(dat_return$Time)
          } # if
          
          # Modification for faceting to plotting dataset
          dat_temp4[, fac := paste0(Var_Name_Base, " ~ ", Var_Name)]
          
          g <- ggplot(dat_temp4, aes(x=Time, y=Sch_Res)) +
               theme_minimal() + xlab(expression(Default~Time~italic(tau[d]))) + ylab(expression(italic(hat(Sc)(x[k~i~tau[d]])))) +
               theme(text=element_text(family=chosenFont),legend.position = "bottom",
                     strip.background=element_rect(colour="grey", fill="#D3D3D3")) +
               facet_wrap(~ fac, strip.position = "right") +
               geom_point(shape=1, col="#3F702F") + geom_smooth(aes(linetype="solid"), col="#043927", fill="#043927") +
               scale_x_continuous(limits = c(NA, max_time)) +
               scale_linetype_manual(name = "Loess-smoother", values = "solid", labels = NULL)
          g_name <- paste0("SchRes_",var[i], levels[j])  
          gplots[[g_name]] <- g; names(gplots[k])    
        } # if
        
        # Correlation Test
        cor_test[[k]] <- cor.test(dat_temp4$Time, dat_temp4$Sch_Res, method = "spearman"); names(cor_test)[k] <- paste(var[i], levels[j]) # Spearman rank correlation is used for robustness
      }
    } # else if
    
  } # for
  
  dat_return_wider <- dat_return %>% pivot_wider(names_from = c(Var_Name_Base, Var_Name), values_from = c(Var_Val, Sch_Res)) %>% setDT()
  
  # Small correction in naming
  names(dat_return_wider)[grep("Var_Val*", names(dat_return_wider))[names(dat_return_wider)[grep("Var_Val*", names(dat_return_wider))] == paste0("Var_Val", "_", var, "_", var)]] <- paste0("Var_Val_", var[sapply(dat_train2[,..var], is.numeric)])
  names(dat_return_wider)[grep("Sch_Res*", names(dat_return_wider))[names(dat_return_wider)[grep("Sch_Res*", names(dat_return_wider))] == paste0("Sch_Res", "_", var, "_", var)]] <- paste0("Sch_Res_", var[sapply(dat_train2[,..var], is.numeric)])
  
  # Clean up
  suppressWarnings(rm(dat_train2, dat_temp, dat_temp2, dat_temp3, dat_temp4, var, var_type, var_name, col_names, g, g_name, k, levels, levels_n))

  return(list(data = dat_return_wider, CorTest = cor_test, plots = gplots))
  
  # rm(dat_temp, dat_temp2, dat_return, id, var, time, status, max_time, verbose)
} # function

# --- UNIT TEST (Breslow Approximation)
# - Setup
# method <- "breslow"
# cph <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1) ~
#              Principal + Instalment + LN_TPE + slc_pmnt_method,
#              id= PerfSpell_Key, data = dat_train %>% group_by(PerfSpell_Key, TimeInPerfSpell),
#              ties = method)
# cph_bres_scaled_schoenfeld <- residuals(cph, type="schoenfeld") # | Use residuals() for classical Schoenfeld residuals and not Scaled; Schoefeld residuals computed for each level

# - Function Execution
# return <- cph_schoen(cph=cph, id = "PerfSpell_Key", dat_train = dat_train, time = "TimeInPerfSpell", status = "DefaultStatus1", verbose = F, max_time = 240)

# - Comparison of Schoenfeld residuals for a numerical variable
# a <- cph_bres_scaled_schoenfeld[,1]
# b <- return$data[, Sch_Res_Principal]
# ab <- data.table(Sch_Res_Principal_Residuals = a,
#                  Sch_Res_Principal_CustFunc = b)
# all.equal(ab$Sch_Res_Principal_Residuals, ab$Sch_Res_Principal_CustFunc)
### RESULTS:~ TRUE
# - Plot comparison for numeric variable
# par(mfcol=c(1,2)); plot(x=rownames(cph_bres_scaled_schoenfeld), y=cph_bres_scaled_schoenfeld[,1], xlim=c(0,240)); plot(x=return$data$Time, y=return$data$Sch_Res_Principal, xlim=c(0,240))

# - Comparison of Schoenfeld residuals for a categorical variable
# c <- cph_bres_scaled_schoenfeld[,4]
# d <- as.numeric(return$data$`Sch_Res_slc_pmnt_method_Debit Order other bank`)
# cd <- data.table(Sch_Res_slc_pmnt_method_Debit_Order_Other_Bank_Residuals = c,
#                  Sch_Res_slc_pmnt_method_Debit_Order_Other_Bank_CustFunc = d)
# all.equal(cd$Sch_Res_slc_pmnt_method_Debit_Order_Other_Bank_Residuals,cd$Sch_Res_slc_pmnt_method_Debit_Order_Other_Bank_CustFunc)
### RESULTS:~ TRUE
# - Plot comparison for categorical variable
# par(mfcol=c(1,2)); plot(x=rownames(cph_bres_scaled_schoenfeld), y=cph_bres_scaled_schoenfeld[,8], xlim=c(0,240)); plot(x=return$data$Time, y=return$data$Sch_Res_slc_pmnt_method_Suspense, xlim=c(0,240))


# --- function to return the time-dependent ROC curve and time-dependent AUC value (Incidence/ Dynamic) of a Cox survival model
# Input:  dat_explain - The dataset containing the required variables for time-dependent ROC and AUC computation
#         predict.time - Time point for valuation
#         Input_Names - A vector containing the names in dat_explain
#             First element is the name of the column containing the event time
#             Second element is the name of the column containing the event indicator
#          marker - A vector containing the predictions of the Cox Model on the records in dat_explain
#             Third element is the name of the column containing the marker values (predictions); we work specifically with the linear prediction
# Output: A list containing time-dependent marker values, TPs, FPs, and AUCs.
risksetROC_helper <- function(dat_explain, Input_Names = c("TimeInPerfSpell", "DefaultStatus1"), marker = cph_lp, predict.time=12) {
  risksetROC::risksetROC(Stime        = dat_explain[[Input_Names[1]]],  
                         status       = dat_explain[[Input_Names[2]]],                 
                         marker       = marker,                             
                         entry        = dat_explain[[Input_Names[1]]]-1,                              
                         predict.time = predict.time,
                         method       = "Cox",
                         plot = FALSE)
}




### AB [2024-11-28]: Similarly, I am unsure of this part's utility

# --- function for identifying FALSE default spells in a given dataset
# Input: dat_given - The dataset in which false default spells should be identified
#        LoanID - Name of the account ID
#        PerfSpellID - Name of the performance spell ID of an account
#        DefSpellID - Name of the default spell ID of an account
#        Counter - Name of the counter variable (counts the observation of a account)
#        PerfSpell_Counter - Name of the counter variable for an assocaited performance spell
#        DefSpell_Counter - Name of the counter variable for an assocaited default spell
#        DefSpellResol_Type_Hist - Name of variable indicating how the default spell was resolved
#        PerfSpell_Max_Date - Name of variable indicating the last observed date of the assocaited performance spell (optional)
#        DefSpell_Max_Date - Name of variable indicating the last observed date of the assocaited performance spell (optional)
# Output: dat_given - The dataset conatining the variable identifying FALSE default spells
False_Perf_Def <- function(dat_given, LoanID=NA, Date=NA, PerfSpellID=NA, DefSpellID=NA, Counter=NA, PerfSpell_Counter=NA, DefSpell_Counter=NA,
                           PerfSpell_Max_Date=NA, DefSpell_Max_Date=NA){
  # --- Unit test parameters:
  # dat_given <- copy(datCredit_smp); LoanID <- "LoanID"; Date <- "Date"; PerfSpellID <- "PerfSpell_Key"; DefSpellID <- "DefSpell_Key"
  # Counter <- "Counter"; PerfSpell_Counter <- "PerfSpell_Counter"; DefSpell_Counter <- "DefSpell_Counter"; DefSpell_Max_Date <- "DefSpell_Max_Date"; PerfSpell_Max_Date <- "PerfSpell_Max_Date"
  
  # --- Arranging the dataset according to the LoanID and Date
  dat_given <- arrange(dat_given, get(deparse(substitute(LoanID))), get(deparse(substitute(Date)))) %>% setDT(key=c(get(deparse(substitute(LoanID))), get(deparse(substitute(Date)))))
  
  # --- Creating a subset containing only the required column names
  colnames <- c(LoanID, Date, PerfSpellID, DefSpellID, Counter, PerfSpell_Counter, DefSpell_Counter,
                ifelse(!is.na(DefSpell_Max_Date), DefSpell_Max_Date, NA), ifelse(!is.na(PerfSpell_Max_Date), PerfSpell_Max_Date, NA))
  dat_sub <- subset(dat_given, select=colnames[!is.na(colnames)])
  
  # --- Renaming the columns to enable easier coding
  colnames_new <- c("LoanID", "Date", "PerfSpellID", "DefSpellID", "Counter", "PerfSpell_Counter", "DefSpell_Counter",
                    ifelse(!is.na(DefSpell_Max_Date), "DefSpell_Max_Date", NA), ifelse(!is.na(PerfSpell_Max_Date), "PerfSpell_Max_Date", NA))
  colnames(dat_sub) <- colnames_new[!is.na(colnames_new)]
  
  # --- Creating a variable showing the previous/ next counter value
  dat_sub[, Counter_Prev := as.numeric(shift(x=Counter, n=1, type="lag")), by=LoanID] # Creating a temporary variable for checking if the next observation of the account is in the dataset (exists)
  dat_sub[, Counter_Next := as.numeric(shift(x=Counter, n=1, type="lead")), by=LoanID] # Creating a temporary variable for checking if the next observation of the account is in the dataset (exists)
  
  # --- Creating variables for indicating whether the previous/ next record of an account exists
  dat_sub[, Prev_Exist := Counter_Prev==Counter-1]
  dat_sub[is.na(Prev_Exist), Prev_Exist := F] # Checking whether this variable is missing (should then be FALSE)
  
  dat_sub[, Next_Exist := Counter_Next==Counter+1]
  dat_sub[is.na(Next_Exist), Next_Exist := F] # Checking whether this variable is missing (should then be FALSE)
  
  # --- Creating a variable for identifying whether a performance/ default spell should be included in an subsequent analysis
  dat_sub[!is.na(DefSpellID) & !(is.na(PerfSpellID)) & Prev_Exist==F, PerfSpell_F := T] # Identifying all instances where a loan is in default (and performance) and the previous observation doesn't exist
  dat_sub[!is.na(DefSpellID) & !is.na(PerfSpellID) & Next_Exist==F, DefSpell_F := T] # Identifying all instances where a loan is in default (and performance) and the next observation doesn't exist
  
  # - Amending the variables for identifying FALSE performance/ default spells
  dat_sub[is.na(PerfSpell_F), PerfSpell_F := F] # Correcting for all other instances of performance spells (since the FALSE ones have already been identified)
  dat_sub[is.na(DefSpell_F), DefSpell_F := F]
  
  # - Amending the default variable for single observation default spells
  dat_sub[DefSpell_Counter==1 & Date==DefSpell_Max_Date & Prev_Exist==T, DefSpell_F := F] #  Because of the overlap, this one observation will be in both the training and validation dataset. This ensures that only one of the observations is chosen to take into account.
  
  # --- Subsetting dat_sub to only include the required variables (that which is to be returned)
  dat_sub <- subset(dat_sub, select=c("LoanID", "Date", "PerfSpell_F", "DefSpell_F")) %>% setDT(key=c("LoanID", "Date"))
  
  # --- Adding the false performance- and default variables to the given dataset
  dat_given <- cbind(dat_given, dat_sub[, list(PerfSpell_F, DefSpell_F)])
  
  # --- Amending the dataset so that FALSE performance- and default spells do not have an associated spell key, counter, and max date (enables easier subsettin)
  dat_given[PerfSpell_F==T, c((PerfSpellID), (PerfSpell_Counter), (PerfSpell_Max_Date)) := lapply(.SD, function(x) {x=NA}), .SDcols = c((PerfSpellID), (PerfSpell_Counter), (PerfSpell_Max_Date))]
  dat_given[DefSpell_F==T, c((DefSpellID), (DefSpell_Counter), (DefSpell_Max_Date)) := lapply(.SD, function(x) {x=NA}), .SDcols = c((DefSpellID), (DefSpell_Counter), (DefSpell_Max_Date))]

  # --- Returning the dataset
  return(dat_given)

}

# --- Checks
# dat_train1 <- False_Pef_Def(dat_train1, LoanID="LoanID", Date="Date", PerfSpellID="PerfSpell_Key", DefSpellID="DefSpell_Key",
#                             Counter="Counter", PerfSpell_Counter="PerfSpell_Counter", DefSpell_Counter="DefSpell_Counter",
#                             DefSpell_Max_Date="DefSpell_Max_Date", PerfSpell_Max_Date="PerfSpell_Max_Date", DefSpellResol_Type_Hist="DefSpellResol_Type_Hist")
# dat_valid1 <- False_Pef_Def(dat_valid1, LoanID="LoanID", Date="Date", PerfSpellID="PerfSpell_Key", DefSpellID="DefSpell_Key",
#                             Counter="Counter", PerfSpell_Counter="PerfSpell_Counter", DefSpell_Counter="DefSpell_Counter",
#                             DefSpell_Max_Date="DefSpell_Max_Date", PerfSpell_Max_Date="PerfSpell_Max_Date", DefSpellResol_Type_Hist="DefSpellResol_Type_Hist")
# 
# lookup_IDs <- unique(datCredit_real[PerfSpell_Num>=5, LoanID])
# 
# lookup <- datCredit_real[LoanID=="3000002499333", ]
# lookup_t <- dat_train1[LoanID=="3000002499333", ]
# lookup_v <- dat_valid1[LoanID=="3000002499333", ]
# lookup_sub <- dat_sub[LoanID==lookup_IDs[200], ]
# 
# # - Overlaps
# overlaps <- which(dat_train2[DefSpell_Counter==1 & DefSpell_F==F, DefSpell_Key] %in% dat_valid2[DefSpell_Counter==1 & DefSpell_F==F, DefSpell_Key])
# 
# lookup_IDs <- str_sub(dat_train2[DefSpell_Counter==1 & DefSpell_F==F, DefSpell_Key][overlaps], start=1, end=-3)
# 
# lookup <- datCredit_real[LoanID==lookup_IDs[10], ]
# lookup_t <- dat_train2[LoanID==lookup_IDs[10], ]
# lookup_v <- dat_valid2[LoanID==lookup_IDs[10], ]


# --- Checks
# dat_train1 <- False_Pef_Def(dat_train1, LoanID="LoanID", Date="Date", PerfSpellID="PerfSpell_Key", DefSpellID="DefSpell_Key",
#                             Counter="Counter", PerfSpell_Counter="PerfSpell_Counter", DefSpell_Counter="DefSpell_Counter",
#                             DefSpell_Max_Date="DefSpell_Max_Date", PerfSpell_Max_Date="PerfSpell_Max_Date", DefSpellResol_Type_Hist="DefSpellResol_Type_Hist")
# dat_valid1 <- False_Pef_Def(dat_valid1, LoanID="LoanID", Date="Date", PerfSpellID="PerfSpell_Key", DefSpellID="DefSpell_Key",
#                             Counter="Counter", PerfSpell_Counter="PerfSpell_Counter", DefSpell_Counter="DefSpell_Counter",
#                             DefSpell_Max_Date="DefSpell_Max_Date", PerfSpell_Max_Date="PerfSpell_Max_Date", DefSpellResol_Type_Hist="DefSpellResol_Type_Hist")
# 
# lookup_IDs <- unique(datCredit_real[PerfSpell_Num>=5, LoanID])
# 
# lookup <- datCredit_real[LoanID=="3000002499333", ]
# lookup_t <- dat_train1[LoanID=="3000002499333", ]
# lookup_v <- dat_valid1[LoanID=="3000002499333", ]
# lookup_sub <- dat_sub[LoanID==lookup_IDs[200], ]
# 
# # - Performance spells
# # The start of default spells overlaps with the end of performance spells. When sub-setting for performance spells, those overlaps are including (although the entire performance spell is not included, only the last observation).
# # Need to identify the observations with no further performance spell history as they are false.
# 
#
# # - Overlaps
# overlaps <- which(dat_train2[DefSpell_Counter==1 & DefSpell_F==F, DefSpell_Key] %in% dat_valid2[DefSpell_Counter==1 & DefSpell_F==F, DefSpell_Key])
# 
# lookup_IDs <- str_sub(dat_train2[DefSpell_Counter==1 & DefSpell_F==F, DefSpell_Key][overlaps], start=1, end=-3)
# 
# lookup <- datCredit_real[LoanID==lookup_IDs[10], ]
# lookup_t <- dat_train2[LoanID==lookup_IDs[10], ]
# lookup_v <- dat_valid2[LoanID==lookup_IDs[10], ]


