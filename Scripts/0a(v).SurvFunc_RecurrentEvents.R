# ============================== SURVIVAL FUNCTIONS ==============================
# Defining custom functions used across various projects that include recurrent
# events.
# --------------------------------------------------------------------------------
# PROJECT TITLE: Default Survival Modelling
# SCRIPT AUTHOR(S): Bernard Scheepers, Dr Arno Botha

# VERSION: 1.0 (November-2024)
# DESCRIPTION: 
# This script defines various functions specific to survival modelling
# that are used elsewhere in this project or, indeed, used across other projects.
# Functions are grouped thematically.
# ================================================================================

# ----------------- 0. Dataset for unit tests ---
# 'The cgd dataset from the survival package contains survival data from a clinical
# trial on patients with chronic granulomatous disease (CGD), a rare immune deficiency.
Test <- FALSE # Toggle for unit tests; Test <- T
if (Test){
  force(data(cgd,package="survival"))
  data(cgd) # Load data set
  dat <- as.data.table(cgd)[, .(id, tstart, tstop, status, sex, age, height, weight, inherit, enum, steroids, treat)] %>% 
            rename(LoanID=id, Start=tstart,End=tstop,Default_Ind=status) #HW:  Generalize variable names.
  #dat <- survSplit(Surv(Start,End,Default_Ind) ~  .,data=cgd,cut=c(1:max(cgd$End)),
  #                start="Start",end="End",event="Default_Ind") %>% as.data.table() # Apply the counting process
  coxExample <- coxph(Surv(Start,End,Default_Ind) ~ sex + age + height + weight + inherit + enum + steroids + treat,data=dat, id=LoanID) # Build a cox model
  summary(coxExample)
}

# ----------------- 1. Cox-Snell residuals analysis ---

# Function to calculate Cox-Snell residuals adjusted for censoring.
# Input: cox - Cox proportional hazard model.
# Output: cs - Cox-Snell residuals

cs_adjusted <- function(cox, dat){
  cs <-  dat[Removed==1,Default_Ind] - 
    residuals(cox,type="martingale",collapse=dat$LoanID) +
    log(2)*(1 - dat[Removed==1,Default_Ind]) # Add log(2) to all observations that have a 0.
  return(cs)
}
# Function to compute the Kolmogorov-Smirnov statistic (1-KS) for Cox-Snell 
# residuals as well as a ggplot graph to display it.
# Input: cox - Cox proportional hazard model.
# Output: KS_stat - 1 - Kolmogorov-Smirnov statistic
#         KS_graph -  Graph of the Cox-Snell empirical cumulative distribution
#                     function and the unit exponential distribution function.

cs_ks_test <- function(cox, dat, GraphInd=T, legPos=c(0.5,0.5)) {
  # Obtain adjusted Cox-Snell residuals
  cs <- cs_adjusted(cox, dat)
  
  # Initialize null distribution
  exp <- rexp(length(cs),1)
  
  # Perform the Kolmogorov-Smirnov test
  KS <- round(suppressWarnings(ks.test(cs,exp))$statistic,4)
  
  # Code to create ggplot object
  if(GraphInd==T){
    # Get the ECDFs of cs
    EmpDist <- ecdf(cs)
    
    # Create a grid of x values for plotting
    x <- sort(unique(c(cs, exp)))
    
    # Calculate CDF values for each sample at each x value
    y1 <- EmpDist(x)
    y2 <- pexp(x,1)
    
    # Find the maximum difference (D statistic)
    D_location <- which.max(abs(y1 - y2))
    
    # Create a data frame for plotting
    datGraph <- data.frame(x = x, cs = y1, exp = y2)
    segment_data <- data.frame(x = x[D_location],xend = x[D_location],
                               y = y1[D_location],yend = y2[D_location],type="Difference")
    # 
    datplot <- rbind( data.table(x=cs,type="1_Cox-Snell"),
                      data.table(x=exp,type="2_Unit_Exponential"))
    vCol <- brewer.pal(8,"Set1")[c(2,3)]
    vLabel <- c("1_Cox-Snell"=bquote("Adjusted Cox-Snell Residual "*italic(r)^(cs)),
                "2_Unit_Exponential"="Unit Exponential")
    dpi <- 
    
    # Plot the ECDFs with ggplot2
    gg <- ggplot(datplot,aes(x=x,group=type)) + theme_minimal() + 
        theme(text = element_text(family="Cambria"), legend.position.inside=legPos,
              legend.position = "inside",
              legend.background = element_rect(fill="snow2", color="black", linetype="solid")) +
        labs(x = "x", y = "Cumulative Distribution Function") +
        stat_ecdf(aes(color=type,linetype=type)) + 
        geom_segment(data=segment_data,aes(x = x, xend = xend, y = y, yend = yend),
                     linetype = "dashed", color = "black") +
        annotate("label", x = x[D_location], y = (y1[D_location] + y2[D_location]) / 2,
                 label = paste("D =", percent(KS)), hjust = -0.1, vjust = -0.1, fill="white", alpha=0.6) +
        scale_color_manual(name = "Distributions", values = vCol, labels=vLabel) +
        scale_linetype_discrete(name = "Distributions",labels=vLabel) +
        scale_y_continuous(label=percent)
    
    # Save figure
    ggsave(gg, file=paste0(genFigPath, "TFD/Kolmogorov-Smirnov/KS",".png"), width=2550/dpi, height=2000/dpi, dpi=dpi, bg="white")
    
    # Prepare return object
    retOb <- list(Stat = as.vector(1 - KS), KS_graph=gg)
  }else{
    retOb <- list(Stat = as.vector(1-KS))
  }
  
  return(retOb)
}

# HW: Unit Test (mine vs standard)

# Function to graphically test the Cox-Snell residuals by plotting them against 
# their respective hazard rate. The line should tend towards the 45 degree line for
# a good fit.
# Input: cox - cox proportional hazard model
# Output: Graph - ggplot object to showcase the relationship

cs_graph <- function(cox){
  # Obtain adjusted Cox-Snell residuals
  cs <- cs_adjusted(cox)
  
  # Create data for graph
  datGraph <- survfit(coxph(Surv(cs, cox$y[, "status"]) ~ 1, method = "breslow"), type = "aalen") %>% tidy() %>%
    mutate(cumu_hazard = -log(.data$estimate)) %>% rename(coxsnell = "time", survival = "estimate") %>%
    subset(select=c("coxsnell","cumu_hazard"))
  
  # Compile ggplot graph of Cox-Snell residuals against their hazard function.
  Graph <-  ggplot(datGraph,aes(x=coxsnell, y=cumu_hazard )) + geom_point() +
    geom_step() + xlab(bquote("Adjusted Cox-Snell Residual "*italic(r)^(cs))) +
    ylab("Cumulative Hazard Function") + 
    geom_abline(slope=1,intercept=0,color='red',linetype="dashed", linewidth=1) + geom_point() + geom_step() + theme_minimal() +
    theme(text = element_text(family="Cambria"))
  
  # Return ggplot object
  return(Graph)
}

# # Unit test
# # cgd dataset: Data from a study on chronic granulomatous disease (CGD), focusing
# # on repeated infections in patients.
# 
# # Load dataset
# data(cgd)
# 
# # Fit a cox model
# coxExample <- coxph(Surv(tstart,tstop,status) ~ sex + age + height + weight,cgd)
# 
# # Test cs_ks_test function
# csResult <- cs_ks_test(coxExample,T)
# csResult$KS_stat;csResult$KS_graph
# ### RESULTS: D=0.1921
# 
# # Test cs_graph function
# cs_graph(coxExample)
# 
# # House keeping
# rm(cgd,cgd0,coxExample)

# # p <- ggplot(datGraph, aes(x = x)) +
# geom_line(aes(y = cs, color = "Residuals")) +
#   geom_line(aes(y = exp, color = "Exponential")) +
#   geom_segment(data=segment_data,aes(x = x, xend = xend, 
#                                      y = y, yend = yend),
#                linetype = "dashed", color = "black") + # Create the maximum distance line
#   annotate("label", x = x[D_location], y = (y1[D_location] + y2[D_location]) / 2,
#            label = paste("D =", 1-KS), hjust = -0.1, vjust = -0.1, fill="white", alpha=0.6) + # Add the distance label
#   labs(x = bquote("Adjusted Cox-Snell Residual "*italic(r)^((cs))), y = "Cumulative Distribution Function") +
#   scale_color_manual(name = "Distributions", values = c("Residuals" = "#4DAF4A", "Exponential" = "#377EB8")) +
#   theme_minimal() + theme(text = element_text(family="Cambria"), legend.position.inside = c(1,0),
#                           legend.justification = c(1,0), legend.background = element_rect(fill = "white", color = "black", linewidth = 0.5))
# HW:
# ----------------- 2. Time-Dependent ROC analysis ---

# Function to graph the time dependent ROC curve and calculate the AUC.
# Input:  dat - validation dataset containing the variables for testing accuracy.
#         cox - cox model used to obtain linear predictor marker values.
#         month_Start - starting month included in time range over which accuracy is tested on.
#         month_End - last month included in time range over which accuracy is tested on.
#         Graph - toggle output to not include the ggplot object.
#         lambda - smoothing parameter representing % neighborhood size symmetrically around a each unique marker.
#         numDigits - rounding scheme applied to markers (specifically if the vectors are memory intensive).
#         method - estimator method to be used to obtain TPR and FPR.
#         name - name for ggplot graph to be saved.
# Output: AUC - Area under the curve
#         ROC_graph - ggplot object for ROC curve

tdROC <- function(dat, cox, month_Start=0, month_End, Graph=TRUE, lambda=0.05, numDigits=2, method="NNE-0/1", Name="tdROC"){
  # Error handling
  if (!is.data.table(dat)) {
    stop("Error: 'dat' must be a data table.")
  }# Test whether dat is a data table
  if (!inherits(cox, "coxph")) {
    stop("Error: 'cox' must be a valid 'coxph' model object.")
  }# Test whether [cox] is a coxph model
  if (!all(all.vars(formula(cox)) %in% colnames(dat))){
    stop("Error: 'dat' does not contain the variables in the 'cox' object")
  }# Test whether [dat] contains the variables on which [cox] was built
  if (!is.numeric(month_Start) || !is.numeric(month_End)) {
    stop("Error: 'month_Start' and 'month_End' must be numeric.")
  }# Test whether [month_Start] and [month_End] are numerical
  if (month_Start < 0 || month_End < 0) {
    stop("Error: 'month_Start' and 'month_End' must be non-negative.")
  }# Test whether [month_Start] and [month_End] are positive
  if (month_Start > month_End) {
    stop("Error: 'month_Start' cannot be greater or equal to 'month_End'.")
  }# Test whether [month_Start] is less than [month_End]
  
  # Unit test
  # dat = dat;cox=coxExample;month_End=12;numDigits=2;Graph=TRUE;month_Start=0;method="NNE-0/1";lambda=0.05;
  
  # Initialize data set
  dat[, Marker := round(predict(cox, newdata=dat, type="lp"),numDigits)] # Create marker values based on linear predictors
  thresholds <- dat$Marker %>% unique() %>% sort() # Let the unique marker values represent different thresholds
  nThresh <- length(thresholds) # number of thresholds for the ROC curve
  
  # Sort data set according to time
  setorder(dat,End) # FOR AB -> make programmable (Stop time = ENd)

  # Initialize Nearest Neighbor Estimation variables
  uEnd <- dat$End %>% unique() %>% sort() # Unique time points in the [dat] data table
  #nTimes <- sum(uEnd <= month) # Number of months before and including the final month
  #uMarker <- dat$Marker %>%  unique() %>% sort() # Obtain unique markers
  uDTime <- dat$End[dat$Default_Ind == 1] %>% unique %>% sort() # Distinct time points in the [dat] data table where the target event takes place
  DTimes <- uDTime[uDTime >= month_Start & uDTime <= month_End] # Distinct time points within range [mont_Start,month_End] in the [dat] data table where the target event takes place
  S_t <- numeric(nThresh) # Initialize vector to contain survival estimates
  n <- NROW(dat) # Total number of markers
  weights <- rep(NA, nrow = n);gc() # (n (# observations in [dat]) x nThresh (# unique markers)) weight matrix containing 1/0's indicating neighbourhoods 
  survivalROC
  # -  Create an indication matrix for subejects at risk at different times
  # NOTE: Matrix dimensions are (n (# observations) x DTimes (unique default times within range))
  start_before_time <- outer(dat$Start, DTimes, "<=")# Indication matrix for observations originating before each DTimes (columns) time point
  end_after_time <- outer(dat$End, DTimes, ">=")# Indication matrix for observations resolving after each DTimes (columns) time point
  at_risk <- (start_before_time & end_after_time) # Indication matrix for observations in each column originating before DTimes (columns) and resolving after column DTimes
  
  end_at_time <- outer(dat$End, DTimes, "==")# Indication matrix for observations currently at each DTimes (columns) time point
  events <- (start_before_time & end_at_time & (dat$Default_Ind == 1)) # Indication matrix for observations resolving at each DTimes (columns) time point
  
  if(method=="NNE-0/1"){
    # Loop through unique markers -> FOR AB: LOGIC FOR 0/1 neighborhood (TRY TO COLLAPSE)
    for (j in 1:nThresh) { # FOR AB: Make Parallel
      # Create vector with entries indicating how close each marker entry is the looping unique marker value
      Diff <- dat$Marker - thresholds[j] # Calculate the difference between unique markers and each other marker
      sDiff <- Diff[order(Diff)]  # Sort the difference vector to ensure that the ordinal vector entries are in ascending order
      
      # Establish a symmetrical neighborhood  around each unique marker value
      Neigh_Mid <- sum(sDiff <= 0) # Index for Marker value equal to looping unique marker value (last marker index for identical markers), since ordinal vector will contain initial differenced marker values, i.e. less than 0
      Neigh_UpperB_ind <- min(Neigh_Mid + trunc(n * lambda * 0.5), n) # Index for upper bound of neighbourhood bounded by the index of largest marker
      Neigh_LowerB_ind <- max(Neigh_Mid - trunc(n * lambda * 0.5), 1) # Index for lower bound of neighbourhood bounded by the index of smallest marker
      Neigh_UpperB <- sDiff[Neigh_UpperB_ind] # Marker value for upper bound of neighbourhood
      Neigh_LowerB <- sDiff[Neigh_LowerB_ind] # Marker value for lower bound of neighbourhood
      
      # Weight vector containing 1 for when the marker is within the neigbourhood and 0 otherwise
      weights <- ifelse((Diff <= Neigh_UpperB) & (Diff >= Neigh_LowerB), 1,0) # Determine whether markers are within a neighbourhood by ascertaining whether it falls within the bounds
      
      # Initialize values for Kaplan-meier estimate of survival function
      n_values <- colSums(weights * at_risk) # number of observations at risk in each neighborhood per unique event time
      d_values <- colSums(weights * events) # number of observations resolving to target event in each neighborhood per unique event time
      
      # Estimate survival function provided the marker is equal to the unique marker
      survival_factors <- 1 - (d_values / n_values) # Kaplan-Meier estimate factor at each DTimes period
      survival_factors[is.na(survival_factors)] <- 1  # Set NaN cases to 1
      S_t[j] <- prod(survival_factors) # Take the product of each survival factor to obtain the survival probability at time t given marker is equal to the unique marker
    }
  } else if(method=="NNE-Exp"){ # Method using exponential kernel function
    warning("Method untested")
    weights <- lapply(thresholds,function(thresh) exp(-(dat$Marker - thresh)^2 / lambda^2)) # Compute exponential weights for all thresholds values
    
    # Calculate weighted populations and events
    results <- sapply(1:nThresh, function(i) {
      n_value <- colSums(weights[[i]] * at_risk) # Weighted populations at risk for column j
      d_value <- colSums(weights[[i]] * events) # Weighted events for column j
      
      # Calculate survival factor for column j
      survival_factor <- 1 - (d_value / n_value)

      # Handle division by zero
      if (is.na(survival_factor)) survival_factor <- 1
      
      S_t[i] <- prod(survival_factor)
      
      return(S_t) # Return the survival factor for this column
    })
  } else{
    stop("Unknown method.")
  }

  # Allocate survival probability at time t give a specific marker value to their corresponding marker value
  dat[,Surv_prob := S_t[match(dat$Marker, thresholds)]]
  
  # Survival probability at time t regardless of the marker value
  S_Overall <- mean(dat[,list(S_Marg = sum(Surv_prob,na.rm=T)/.N), by=list(PerfSpell_Key)]$S_Marg) # Calculate the average [Surv_prob] for each "id" and averaging these "id" specific averages across the portfolio
  
  # Initialize ROC matrix
  roc.matrix <- matrix(NA, nThresh,2)
  roc.matrix[nThresh, ] <- c(0, 1) # Initialize matrix to start off with 0 sensitivity (TPR) and 1 Specificity (1-FPR)
  
  # Populate matrix with sensitivity and Specificity values
  for (c in 1:(nThresh - 1)) {
    # Empirical distribution of markers being less than the threshold
    cumulMark = mean(dat[,list(sum(Marker <= thresholds[c])/.N), by=list(PerfSpell_Key)]$V1) # First average the markers being less than the threshold for each "id" before averaging the id-averages, i.e. aggregate to the entire dataset
    
    #cumulMark <- sum(dat$Marker <= thresholds[c])/n # Number of observations with a Marker value less < threshold divided by observations
    #S_lam <- sum(dat$Surv_prob[dat$Marker > thresholds[c]])/n # Sum of survival probabilities for Marker values greater than threshold
    
    # Survival probability at time t provided that the corresponding marker values are greater than the threshold
    threshSurv <- mean(dat[,list(sum(ifelse(Marker > thresholds[c],Surv_prob,0),na.rm=T)/.N), by=list(PerfSpell_Key)]$V1,na.rm=T) # First average the Surv_prob values with markers greater than the threshold for each "id" before averaging the id-averages, i.e. aggregate to the entire dataset
    
    # Populate roc matrix
    roc.matrix[c, 1] <- ((1-cumulMark) - threshSurv)/(1 - S_Overall) # Sensitivity
    roc.matrix[c, 2] <- 1 - threshSurv/S_Overall # Specificity
    #roc.matrix[c,3] <- cumulMark
    #roc.matrix[c,4] <- S_lam
  }
  # Convert Sensitivity and Specificity to TPR and FPR respectively
  sensitivity = roc.matrix[, 1]
  specificity = roc.matrix[, 2]
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
    conc=percent(as.numeric(concordance(cox1)[1]))
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
                   conc), fill="grey") + 
      scale_y_continuous(label=percent) + scale_x_continuous(label=percent)
    
    # Save graph
    ggsave(gg, file=paste0(genFigPath, "TFD/tdROC/",Name,"(",month_Start,",",month_End,").png"), width=1200/dpi, height=1000/dpi, dpi=dpi, bg="white")
    
    retObj <- list(AUC = area, ROC_graph=gg)
  }else{
    retObj <- list(AUC = area)
  }
  return(retObj)
}

if(Test){
  survivalROC(Stime=dat$Start,status=dat$Default_Ind, entry=dat$Start,
              marker=round(predict(coxExample, type="lp"),2),span=0.05, predict.time=12)
  tdROC(dat,coxExample,12)
  # tdROC(datCredit_valid_TFD, coxDelinq, month_Start=0, month_End=3, Graph=TRUE, 
  #       lambda=0.05, numDigits=0, method="NNE-0/1", Name="coxDelinq",
  #       fld_id="PerfSpell_Key", fld_event="Default_Ind")
}

# ----------------- 2. Schoenfeld residuals ---

# Function to graph the time dependent ROC curve and calculate the AUC.
# Input:  dat - Dataset containing the [Start], [Stop] and [Default_Ind] variables
#         cox - cox model
#         month - desired month to test ROC curve on.
#         span - % size of the neighborhood (will be symmetrical around the same point)
#         numDigits - rounding scheme applied to markers (specifically if the vectors memory consuming)
# Output: AUC - Area under the curve
#         ROC_graph - ggplot object for ROC curve

sfResiduals <- function(cox, dataset, var, legPos = c(50,1), legPosCat = c(0.9,0.1)){
  # Select relevant columns and convert to data.table
  dat <- dataset[, .(LoanID, Start, End, Default_Ind, Var=get(var))]
  
  # Add risk score and total risk score
  dat[, RiskScore := predict(cox, dataset, type = "risk")]
  dat[, TotalScore := sum(RiskScore), by = End]
  
  # Process numeric and categorical variables differently
  if (is.numeric(dat$Var)) {
    # Handling numeric Var
    dat[, RW_Val := Var * RiskScore / TotalScore]
    dat[, Exp_Val := sum(RW_Val), by = End]
    dat[, sfRes := Var - Exp_Val]
    #dat[, RW_V := var(sfRes), by=End]
    #dat[, ssfRes := sfRes/RW_V]
    p <- sfTest(cox)
    
    # Create a data frame for plotting
    datGraph <- dat[Default_Ind == 1,]
    segment_data <- data.frame(x = 1,xend = max(datGraph$End),y = 0,yend = 0)
    
    gg <- ggplot(datGraph,aes(x=End,y=sfRes)) + theme_minimal() + 
      theme(text = element_text(family="Cambria"),
            legend.position = "inside",
            legend.background = element_rect(fill="snow2", color="black", linetype="solid")) +
      geom_point(alpha=0.7, color = "cornflowerblue") + geom_smooth(method="loess", color="navy") +
      annotate("label",x=legPos[1], y=legPos[2], label = paste("p-value for ",var,": ", percent(p)),
               fill="grey", alpha=0.6) +
      geom_segment(data= segment_data, aes(x = x, xend = xend, y = y, yend = yend),
                   linetype = "dashed", color = "black") +
      labs(x = bquote("Default Time "*italic(T)), y = bquote("Schoenfeld Residuals "*italic(r)^(s)))
    
  } else if (is.character(dat$Var) || is.factor(dat$Var)) {
    # # Handling categorical Var
    Levels <- unique(dat$Var)
    # dat[, Var_Name_Level := as.character(Var)]  # Add level information
    # p <- sfTest(cox)
    # 
    # # Calculate RW_Val, Exp_Val, and residuals for all levels at once
    # dat[, RW_Val := (Var == Var_Name_Level) * RiskScore / TotalScore]
    # dat[, Exp_Val := sum(RW_Val), by = .(End, Var_Name_Level)]
    # dat[, sfRes := (Var == Var_Name_Level) - Exp_Val]
    # 
    if(length(Levels) == 2){
      r <- residuals(cox, type = "schoenfeld")
      r <- data.frame(End=as.numeric(names(r)), residuals = r )
      setnames(r, "residuals", paste0(var,Levels[2]))
      vCol <- brewer.pal(ncol(r)-1,"Set1")[1]
    }else{
      r <- residuals(cox,type="schoenfeld") %>% data.table()
      r <- cbind(End=as.numeric(rownames(r)),r) %>% as.data.table()
      vCol <- brewer.pal(ncol(r)-1,"Set1")
    }
    
    # Create a data frame for plotting
    datGraph <- pivot_longer(r,cols=starts_with(var),names_to=var, values_to = "sfRes")
    segment_data <- data.frame(x = 1,xend = max(datGraph$End),y = 0,yend = 0)
    #vLabel <- separate(unique(datGraph[,var]), strings, into=c("before", "after"), sep=var)$after
    

    gg <- ggplot(datGraph,aes(x=End, y=sfRes, color=get(var))) + theme_minimal() +
      theme(text = element_text(family="Cambria"),
            legend.position = "top",legend.background = element_rect(fill="snow2", color="black", linetype="solid")) +
      labs(x = bquote("Default Time "*italic(T)), y = bquote("Schoenfeld Residuals "*italic(r)^(s))) +
      geom_point(alpha=0.7) + geom_smooth(method="loess", se=FALSE) + facet_wrap( ~ get(var), scales="free_y") +
      geom_segment(data= segment_data, aes(x = x, xend = xend, y = y, yend = yend),
                   linetype = "dashed", color = "black") +
      scale_color_manual(name=var, values=vCol)
    
    # Save graph
    ggsave(gg, file=paste0(genFigPath, "TFD/Shoenfeld Residuals(",var,")",".png"), width=2550/dpi, height=2000/dpi, dpi=dpi, bg="white")
  }
      
  # Return final result
  print(gg)
  
  return(list(p_value = round(p,2),sumRes = sum(datGraph$sfRes)))
}

sfTest <- function(cox){
  ans <- cox.zph(cox)$table[1,"p"]
  return(ans)
}


# r <- residuals(cox,type="schoenfeld")
# plot(names(r), r)
# 
# sr <- residuals(cox,type="scaledsch")
# plot(names(sr),sr)
# abline(0,0, col="red")
