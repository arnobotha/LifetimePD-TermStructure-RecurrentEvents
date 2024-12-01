tROC <- function(datGiven, cox, month_Start=0, month_End, lambda=0.05, method="NN-0/1", numDigits=2, 
                 fld_ID=NA, fld_Event="MainEvent_Ind", eventVal=1, fld_StartTime="Start", fld_EndTime="Stop",
                 Graph=TRUE, graphName="timedROC-Graph", genFigPath=paste0(getwd(),"/")){
  
  # --- Preliminaries 
  # -- Testing Conditions
  datGiven = copy(dat); cox=coxExample; month_End=334; numDigits=2; Graph=TRUE; month_Start=0; method="NNE-0/1";lambda=0.05;
  fld_ID="ID"; fld_Event="Event_Ind"; fld_StartTime="Start"; fld_EndTime="End";
  graphName="coxExample_cgd"; eventVal=1; genFigPath=genFigPath
  
  # -- Error handling
  if (!is.data.table(datGiven)) {
    stop("[datGiven] must be a data table.\n")
  }# Test whether datGiven is a data table
  if (!inherits(cox, "coxph")) {
    stop("[cox] must be a valid 'coxph' model object.\n")
  }# Test whether [cox] is a coxph model
  if (!all(all.vars(formula(cox)) %in% colnames(datGiven))){
    stop("[datGiven] does not contain the variables required within the [cox] object.\n")
  }# Test whether [datGiven] contains the variables on which [cox] was built
  if (!is.numeric(month_Start) || !is.numeric(month_End)) {
    stop("[month_Start] and [month_End] must be numeric.\n")
  }# Test whether [month_Start] and [month_End] are numerical
  if (month_Start < 0 || month_End < 0) {
    stop("[month_Start] and [month_End] must be non-negative.\n")
  }# Test whether [month_Start] and [month_End] are positive
  if (month_Start > month_End) {
    stop("[month_Start] cannot be greater or equal to [month_End].\n")
  }# Test whether [month_Start] is less than [month_End]
  if (anyNA(c(fld_Event, eventVal, fld_StartTime, fld_EndTime))) {
    stop("The arguments [fld_Event], [eventVal], [fld_StartTime], [fld_EndTime], and [lfd_Marker] cannot be missing and must be specified. \n")
  }
  if ((Graph & is.na(graphName)) | (Graph & is.na(genFigPath))) {
    stop("The graphing arguments [graphName] and [genFigPath] cannot be missing and must be specified when desiring an ROC-graph. \n")
  }
  
  # -- Preliminaries
  datGiven[, Marker := round(predict(cox, newdata=datGiven, type="lp"),numDigits)] # Create marker values based on linear predictors
  thresholds <- datGiven$Marker %>% unique() %>% sort() # Let the unique marker values represent different thresholds
  nThresh <- length(thresholds) # number of thresholds for the ROC curve
  
  # -- Reassign given field names to standardised naming conventions, if only within this function
  # NOTE: This will be reverted towards the end of the function
  if (!is.na(fld_ID)){
    setnames(datGiven, old=c(c(fld_StartTime, fld_EndTime, fld_ID, fld_Event)),
             new=c("StartTime", "EndTime", "ID", "Event_Ind"))
    noGroup_Ind <- F # set indicator according to which the rest of the program will function conditionally
  } else {
    setnames(datGiven, old=c(c(fld_StartTime, fld_EndTime, fld_Event)),
             new=c("StartTime", "EndTime", "Event_Ind"))
    noGroup_Ind <- F # set indicator according to which the rest of the program will function conditionally
  }
  
  # - Check if any end points fall within the given range
  if (datGiven[EndTime >= month_Start & EndTime <= month_End, .N] == 0) {
    stop(paste0("The observed endpoints exceed the bounds of the given prediction time interval, [",
                month_Start, ",", month_End, "]"))
  }
  
  
  # -- Obtain various quantities towards implementing the Nearest Neighbour Estimator (NNE) method
  
  # 1. 
  
  # Get the corresponding rank order of the raw end times when sorted ascendingly
  vOrder <- order(datGiven$End)
  # Re-sort the raw end points according to these particular  rank-order indices, whereafter
  # the same is performed to the vectors of event indicators and marker-values
  vEventTimes <- datGiven$EndTime[vOrder]
  vEventInds <- datGiven$Event_Ind[vOrder]
  vMarkers <- datGiven$Marker[vOrder]
  
  
  
  #setorder(datGiven, cols = EndTime) #  Sort data set according to the given end time field
  # Obtain unique end points (or ages) in the [datGiven] object, i.e., the "ordered failure times t_1 < t_2 < ... < t_m"
  #uEnd <- datGiven$EndTime %>% unique %>% sort()
  # Obtain unique end points (or ages) at which the main event of interest occurred
  uDTime <- datGiven$EndTime[datGiven$Event_Ind == eventVal] %>% unique %>% sort() 
  # Obtain unique end points within the range [mont_Start,month_End] during which the main event occurred
  DTimes <- uDTime[uDTime >= month_Start & uDTime <= month_End] 
  
  # - Initialize empty data structures for implementing the NNE-method
  S_t <- numeric(nThresh) # Initialize vector to contain survival estimates
  n <- NROW(datGiven) # Total number of markers (not necessarily subjects)
  weights <- rep(NA, nrow = n) # kernel vector applied as "weights" onto eventual survival estimates
  
  
  # --- Create various indicator matrices towards estimating the survivor functions
  
  # --  Create an indicator matrix for subjects at risk across the spectrum of event times
  # NOTE: Matrix dimensions are all (n (# observed event times) x m (unique event times within range))
  # - Using outer product of two arrays, create an n x m Boolean-valued matrix that indicates whether 
  # observed start points are earlier/equal than/to unique event times
  matStartBefore <- outer(datGiven$StartTime, DTimes, "<=") 
  # - Using outer product of two arrays, create an n x m Boolean-valued matrix that indicates whether 
  # observed end points are later/equal than/to unique event times
  matEndAfter <- outer(datGiven$EndTime, DTimes, ">=") 
  # - Multiply matrices such that when a cell is true, then the observed raw end point t'=1,\dots,n is still at risk at 
  # the unique event time t=1,\dots,m; i.e., If TRUE, then the rawTime-uniqueTime (t',t)-tuple has an at-risk 
  # lifetime as at t, though during which lifetime it experienced the event at some future time t' > t
  matAtRisk_Ind <- (matStartBefore & matEndAfter)
  rownames(matAtRisk_Ind)<- paste0(1:n, ": ", datGiven$EndTime) # relabel rows intuitively
  
  # -- Create an indicator matrix for subjects that experienced the event across the spectrum of event times
  # - Using outer product of two arrays, create an n x m Boolean-valued matrix that indicates whether 
  # observed end points are currently equal to unique event times
  matEndAt <- outer(datGiven$EndTime, DTimes, "==")
  # - Multiply matrices such that when a cell is true, then the observed raw end point t'=1,\dots,n equals
  # the particular unique event time t=1,\dots,m; i.e.,If TRUE, then the rawTime-uniqueTime (t',t)-tuple of the 
  # corresponding subject has experienced the event at exactly t.
  matEventsAt <- (matStartBefore & matEndAt & (datGiven$Event_Ind == eventVal))
  rownames(matEventsAt)<- paste0(1:n, ": ", datGiven$EndTime) # relabel rows intuitively
  
  
  # --- In calculating the ROC-graph, 3 fundamental quantities must be estimate:
  # 1) the classical survivor function S(t) irrespective of Marker values
  # 2) the conditional survivor function S(t| M > c)
  
  # -- 1. Estimating the classical S(t) given each threshold
  # Implement the chosen estimator for S(t) and the choice of kernel (if Nearest Neighbour)
  
  if(method=="NN-0/1"){ # Nearest Neighbour (NN) estimator, using the 0/1 Kernel function from Akritas1994
    
    # -- S(t) is estimated by iterating across unique markers as thresholds, whilst assuming each threshold
    # holds 'universally' across all subjects.
    for (j in 1:69) { # nThresh
      
      # - Create a marker vector where its entries denoting the distance between each marker and the current threshold
      Diff <- datGiven$Marker - thresholds[j] # distance
      sDiff <- Diff[order(Diff)]  # Sort in ascending order
      
      # - Establish a symmetrical neighbourhood around each unique marker value
      # Find an index in the ordered difference vector [sDiff] beyond which point all markers are positively 
      # differenced wrt the current threshold
      Neigh_Mid <- sum(sDiff <= 0) 
      # Find an index for upper bound of neighbourhood bounded by the index of largest marker
      Neigh_UpperB_ind <- min(Neigh_Mid + trunc(n * lambda * 0.5), n)
      # Find an index for lower bound of neighbourhood bounded by the index of smallest marker
      Neigh_LowerB_ind <- max(Neigh_Mid - trunc(n * lambda * 0.5), 1)
      Neigh_UpperB <- sDiff[Neigh_UpperB_ind] # Marker value for upper bound of neighbourhood
      Neigh_LowerB <- sDiff[Neigh_LowerB_ind] # Marker value for lower bound of neighbourhood
      
      ## - Apply chosen kernel function: 1 if the Marker is within the neighbourhood, 0 otherwise
      weights <- ifelse((Diff <= Neigh_UpperB) & (Diff >= Neigh_LowerB), 1,0)
      
      # - Initialize constituent quantities for the eventual Kaplan-Meier estimator of the survivor function S(t)
      # At each unique event time, calculate the number of at-risk raw end points within each neighborhood
      n_values <- colSums(weights * matAtRisk_Ind) 
      # At each unique event time, calculate the number of raw end points within each neighborhood that just experienced 
      # the event
      d_values <- colSums(weights * matEventsAt) 
      
      # - Estimate the survival function S(t) as the nonparametric product-limit (Kaplan-Meier) over unique event times
      # NOTE: These quantities are already filtered such that the TRUE-marked raw end points in 
      # [d_values] equal the current threshold, having used the NNE-method
      # Calculate the KM-based survival estimate at each unique event time, producing a vector
      survival_factors <- 1 - (d_values / n_values)
      survival_factors[is.na(survival_factors)] <- 1  # Set NaN cases to 1
      # Assemble the previous vector into S(t) by taking the product thereof across all survival factors (of raw event times)
      # Implicitly, we assume that the "true" marker is equal to the current threshold at t
      S_t[j] <- prod(survival_factors)
    }
    
  } else if(method=="NN-Exp"){ # NN-estimator using exponential kernel function
    warning("This estimation method [NNE-Exp] is still untested; caution advised.\n")
    weights <- lapply(thresholds,function(thresh) exp(-(datGiven$Marker - thresh)^2 / lambda^2)) # Compute exponential weights for all thresholds values
    
    # Calculate weighted populations and events
    results <- sapply(1:nThresh, function(i) {
      n_value <- colSums(weights[[i]] * matAtRisk_Ind) # Weighted populations at risk for column j
      d_value <- colSums(weights[[i]] * matEventsAt) # Weighted events for column j
      
      # Calculate survival factor for column j
      survival_factor <- 1 - (d_value / n_value)
      
      # Handle division by zero
      if (is.na(survival_factor)) survival_factor <- 1
      
      S_t[i] <- prod(survival_factor)
      
      return(S_t) # Return the survival factor for this column
    })
  } else{ # Fail the execution
    
    # - Revert name changes 
    if (noGroup_Ind==F){
      setnames(datGiven, new=c(c(fld_StartTime, fld_EndTime, fld_ID, fld_Event)),
               old=c("StartTime", "EndTime", "ID", "Event_Ind"))
    } else {
      setnames(datGiven, new=c(c(fld_StartTime, fld_EndTime, fld_Event)),
               old=c("StartTime", "EndTime", "Event_Ind"))
    }
    stop("Unknown estimation method. Exiting ..")
  }
  
  # - Allocate the estimated S(t)-values across subjects over prediction times t (same as unique event times), given
  # a specific marker value under which that particularly S(t) was estimated
  datGiven[,Surv_prob := S_t[match(datGiven$Marker, thresholds)]]
  
  
  # - Estimation of the remaining two quantities depend on whether observations are supposedly independent from one another
  # or whether they are clustered around a common ID-value
  if (noGroup_Ind==F) { # Dependence amongst certain observations
    
    # -- 2. Estimate the conditional survivor function S(t| M > c), given the subset with marker M and cut-off C
    
    # - Calculate the overall survival probability at prediction time t, i.e., given an -\infty marker value
    # Calculate the average [Surv_prob] for each "id" and averaging these "id" specific averages across the portfolio
    # NOTE: This estimator is the grand mean of the average [Surv_prob]-values per ID
    S_Overall <- mean(datGiven[,list(S_Marg = sum(Surv_prob,na.rm=T)/.N), by=list(ID)]$S_Marg) 
    
    # - Initialize an m x 2 results matrix in which true & false positive rates are stored across columns
    # per unique 
    matRates <- matrix(NA, nThresh,2)
    matRates[nThresh, ] <- c(0, 1) # Initialize matrix to start off with 0 sensitivity (TPR) and 1 Specificity (1-FPR)
    
    # Populate matrix with sensitivity and Specificity values
    for (c in 1:(nThresh - 1)) {
      # Empirical distribution of markers being less than the threshold
      cumulMark = mean(datGiven[,list(sum(Marker <= thresholds[c])/.N), by=list(ID)]$V1) # First average the markers being less than the threshold for each "id" before averaging the id-averages, i.e. aggregate to the entire dataset
      
      #cumulMark <- sum(datGiven$Marker <= thresholds[c])/n # Number of observations with a Marker value less < threshold divided by observations
      #S_lam <- sum(datGiven$Surv_prob[datGiven$Marker > thresholds[c]])/n # Sum of survival probabilities for Marker values greater than threshold
      
      # Survival probability at time t provided that the corresponding marker values are greater than the threshold
      S_t <- mean(datGiven[,list(sum(ifelse(Marker > thresholds[c],Surv_prob,0),na.rm=T)/.N), by=list(ID)]$V1,na.rm=T) # First average the Surv_prob values with markers greater than the threshold for each "id" before averaging the id-averages, i.e. aggregate to the entire dataset
      
      # - Populate results-matrix
      matRates[c, 1] <- ((1-cumulMark) - threshSurv)/(1 - S_Overall) # Sensitivity
      matRates[c, 2] <- 1 - threshSurv/S_Overall # Specificity
      
    } else { # Independence amongst observations
      
      # - Calculate the mean survival probability at prediction time t, i.e., given an -\infty marker value
      S_Overall <- mean(datGiven$Surv_prob, na.rm=T)  
      
      ### AB: Busy drafting ..
    }
    
  } else {stop("Implementation halted for ungrouped survival data.\n")}
  
  # Convert Sensitivity and Specificity to TPR and FPR respectively
  sensitivity = matRates[, 1]
  specificity = matRates[, 2]
  x <- 1 - c(0, specificity) # FPR = 1 - Specificity
  y <- c(1, sensitivity) # TPR = Sensitivity
  
  # Trapezoidal rule
  dx <- x[-(nThresh+1)] - x[-1] # Bar width = difference between consecutive points
  mid.y <- (y[-(nThresh+1)] + y[-1])/2 # Bar height = Average between two consecutive points
  area <- sum(dx * mid.y) # AUC
  
  if(Graph){
    # Create a data frame for plotting
    datGraph <- data.frame(x = x, y=y)
    datSegment <- data.frame(x = 0, y = 0, xend = 1, yend = 1)
    conc=percent(as.numeric(concordance(cox)[1]))
    vCol <- brewer.pal(8,"Set1")[c(2)]
    dpi <- 200
    
    # Plot ROC curve
    gg <- ggplot(datGraph,aes(x=x,y=y,group=1)) + theme_minimal() + 
      theme(text = element_text(family="Cambria"), legend.position="inside",
            legend.background = element_rect(fill="snow2", color="black",
                                             linetype="solid")) +
      labs(x = bquote("False Positive Rate "*italic(F^"+")), y = 
             bquote("True Positive Rate "*italic(T^"+"))) + geom_step(color=vCol) +
      geom_segment(data = datSegment,aes(x = x, y = y, xend = xend, yend = yend),
                   color = "grey", linetype = "dashed") +
      annotate("label", x = c(0.75,0.75), y = c(0.375,0.125),label = 
                 c(paste0("AUC(",month_Start,",", month_End,"): ", percent(area)),
                   paste0("Harrell's c-statistic: ", conc)), fill="grey") + 
      scale_y_continuous(label=percent) + scale_x_continuous(label=percent)
    
    # Save graph
    ggsave(gg, file=paste0(genFigPath, graphName,"(",month_Start,",",month_End,").png"), width=1200/dpi, height=1000/dpi, dpi=dpi, bg="white")
    
    retObj <- list(AUC = area, ROC_graph=gg)
  }else{ retObj <- list(AUC = area) }
  
  # - Revert name changes 
  if (noGroup_Ind==F){
    setnames(datGiven, new=c(c(fld_StartTime, fld_EndTime, fld_ID, fld_Event)),
             old=c("StartTime", "EndTime", "ID", "Event_Ind"))
  } else {
    setnames(datGiven, new=c(c(fld_StartTime, fld_EndTime, fld_Event)),
             old=c("StartTime", "EndTime", "Event_Ind"))
  }
  
  # - Conclude program and return results
  return(retObj)
}