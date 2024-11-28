# ============================== SURVIVAL FUNCTIONS ==============================
# Defining various custom functions relating to the estimation, analysis, and
# graphing of ROC-graphs and their summary statistics (AUC), as part of testing
# the discrimination power (prediction accuracy) of a given Cox regression model
# --------------------------------------------------------------------------------
# PROJECT TITLE: Default Survival Modelling
# SCRIPT AUTHOR(S): Bernard Scheepers, Dr Arno Botha
# VERSION: 1.0 (November-2024)
# ================================================================================




# ----------------- 0. Dataset and models for unit tests ---
# The cgd dataset from the survival package contains survival data from a clinical
# trial on patients with chronic granulomatous disease (CGD), a rare immune deficiency.
Test <- FALSE # Toggle for unit tests; Test <- T
if (Test){
  force(data(cgd,package="survival"))
  data(cgd) # Load data set
  # Lightly prepare data into a generic format that can span our eventual credit dataset as well
  dat <- as.data.table(cgd)[, .(id, tstart, tstop, status, sex, age, height, weight, inherit, enum, steroids, treat)] %>% 
    rename(ID=id, Start=tstart,End=tstop,Event_Ind=status)
  #dat <- survSplit(Surv(Start,End,Default_Ind) ~  .,data=cgd,cut=c(1:max(cgd$End)),
  #                start="Start",end="End",event="Default_Ind") %>% as.data.table() # Apply the counting process
  
  # --- Fit Kaplan-Meier (KM) nonparametric (and "empty-of-covariates") model
  # Compute Kaplan-Meier survival estimates (product-limit) for main-event | Spell-level with right-censoring & left-truncation
  # All competing events preclude the main event from happening and are therefore considered as censored
  # ID is set as the spell key, with no stratification
  kmExample <- survfit(Surv(time=Start, time2=End, event=Event_Ind==1,type="counting") ~ 1, 
                       id=ID, data=dat)
  summary(kmExample)$table # overall summary
  ### RESULTS: 76 events, with median survival probability at time 334 \in [280, 373] as a 95% Confidence Interval
  (kmExample_survFitSummary <- surv_summary(kmExample))
  ### RESULTS: Median survival time of 334 has standard error of 5.8%, which is relatively large
  
  
  # --- Graphing survival and related quantities from fitted KM-model | S(t), h(t)
  
  # -- Graphing parameters
  vCol <- brewer.pal(10, "Paired")[c(10)] # for S(t)
  vCol2 <- brewer.pal(10, "Paired")[c(10,9)] # for h(t)
  sSpan <- 0.1; # span for LOESS-smoother in h(t)
  vlabel <- paste0("Loess-smoothed hazard [span: ", sSpan, "]") # for h(t)
  mainEventName <- "CGD"
  chosenFont <- "Cambria"
  
  # -- Survival probability, S(t)=y
  (gsurv1c_a <- ggsurvplot(kmExample, fun="pct", conf.int=T, legend="none", 
                           break.time.by=round(max(kmExample$time)/8), palette=vCol,
                           xlab = bquote(Discrete~time~italic(t)*" (months) in spell: Multi-spell"),
                           ylab = bquote(Survival~probability~"["*.(mainEventName)*"]"*~italic(S(t))*": spell-level (Kaplan-Meier)"), 
                           xlim=c(0, max(kmExample$time)+1), surv.median.line = "hv", censor=F, 
                           ggtheme = theme_bw(base_family=chosenFont), tables.theme = theme_cleantable(),
                           tables.height=0.10, tables.y.text=F, tables.y.text.col=T, risk.table = "abs_pct", risk.table.pos = "out",
                           cumevents=T, cumevents.title="Cumulative number of events", 
                           cumcensor=T, cumcensor.title="Cumulative number of censored observations (incl. competing risks)",
                           risk.table.title = "Number in (% of) sample at risk of main event", font.family=chosenFont, fontsize=2.5))
  
  
  # -- Discrete baseline hazard function: h(t) | Empirical estimation method
  # - create plotting data object
  haz_dat <- data.table(Time=kmExample$time, AtRisk_n=kmExample$n.risk, 
                        Event_n = kmExample$n.event, Censored_n=kmExample$n.censor,
                        hazard=kmExample$n.event/kmExample$n.risk, 
                        CumulHazard = kmExample$cumhaz, #Nelson-Aalen estimator
                        Group="1",Surv_KM = kmExample$surv) %>% 
    filter(Event_n > 0 | Censored_n >0) %>%
    # Discrete-time variants
    mutate(CumulHazard_Disc = -cumsum(log(1-hazard)), Surv_KM_Disc = cumprod(1-hazard)) %>% 
    mutate(Event_KM_Disc = 1-Surv_KM_Disc) %>% as.data.table()
  haz_dat[, Surv_KM_Disc_prev:= shift(Surv_KM_Disc, n=1, type="lag"), by=list(Group)]
  # - create alternative versions for sanity checks
  haz_dat[Time==Time[1], hazard2 := 1 - Surv_KM_Disc]
  haz_dat[Time>Time[1], hazard2 := 1 - Surv_KM_Disc/Surv_KM_Disc_prev]
  # - conduct sanity checks
  all.equal(haz_dat$hazard, haz_dat$hazard2) # Should be TRUE
  all.equal(haz_dat$Surv_KM, haz_dat$Surv_KM_Disc) # Should be TRUE
  all.equal(haz_dat$CumulHazard, haz_dat$CumulHazard_Disc) # usually FALSE
  plot(kmExample$time, haz_dat$CumulHazard - haz_dat$CumulHazard_Disc, type="b")
  ### RESULTS: The discrepancy is very small difference due to estimator method differences
  
  # - Graph object for shorter time, informed by previous graphs
  (gsurv1c_d <- ggplot(haz_dat[Time<=300,], aes(x=Time,y=hazard)) + theme_minimal() +
      geom_line(linetype="solid", colour=vCol2[1]) + geom_point(colour=vCol2[1]) + 
      geom_smooth(aes(colour=Group, fill=Group), se=T, method="loess", span=sSpan, alpha=0.25, linetype="dotted") +
      labs(y=bquote(plain(Estimated~hazard*" function ["*.(mainEventName)*"]"*~italic(h(t))*": spell-level (Kaplan-Meier)")), 
           x=bquote(Discrete~time~italic(t)*" (months) in spell: Multi-spell")) + 
      theme(text=element_text(family=chosenFont),legend.position="bottom") + 
      scale_colour_manual(name="", values=vCol2[2], labels=vlabel) + 
      scale_fill_manual(name="", values=vCol2[2], labels=vlabel) + 
      scale_y_continuous(breaks=breaks_pretty(), label=percent) + 
      scale_x_continuous(breaks=breaks_pretty(n=8), label=comma))
  ### RESULTS: The hazard appears to be near-constant over time, with some notable oscillation over some prediction periods.
  # However, when viewed in tandem with S(t), itself almost a straight downward-sloping line, it makes sense for hazard
  # to be near-flat. The oscillation also seems more pronounced towards later prediction periods than earlier ones.
  
  # -- Save plots
  dpi <- 150 # need to decrease size for risk tables' text
  ggsave(print(gsurv1c_a,newpage=F), file=paste0(genFigPath,"/SurvFig1c_a-", mainEventName,"_Surv-KaplanMeier-SpellLevel-MultiSpell-LatentComp-InclLeftTrunc_Correct.png"),
         width=1200/dpi, height=1000/dpi,dpi=dpi, bg="white")
  dpi <- 180 # reset
  ggsave(gsurv1c_d, file=paste0(genFigPath,"/SurvFig1c_d-", mainEventName,"_Hazard-KaplanMeier-SpellLevel-MultiSpell-LatentComp-InclLeftTrunc_Correct.png"),
         width=1200/dpi, height=1000/dpi,dpi=dpi, bg="white")
  
  
  
  # --- Fit Cox Regression Model correctly, where observations are clustered around a given ID without assuming independence
  coxExample <- coxph(Surv(Start,End,Event_Ind) ~ weight + age + enum + steroids + treat,
                      data=dat, id=ID)
  summary(coxExample)
  
  # --- Fit a Cox regression model incorrectly by ignoring the clustering and heterogeneous variance, assuming that
  # all observations are subject-level (one row per subject) and therefore independent from one another (wrong)
  coxExample_wrong <- coxph(Surv(Start,End,Event_Ind) ~ weight + age + enum + steroids + treat,
                            data=dat)
  summary(coxExample_wrong)
  ### RESULTS: standard errors (as measured using the robust-variant of the log-rank test) are inflated
  #   in this wrongly-fit model, as a result of assuming independence amongst observations within a cluster/spell per subject
}






# ----------------- 1. Functions for conducting time-Dependent ROC-analysis ------

# --- Function to calculate an ROC-graph (True vs false positive rates) for a given prediction time interval,
# as adjusted for right-censoring by estimating both the overall and marker-conditional survivor functions.
# For a given/fitted Cox regression model, this function chiefly implements the Nearest Neighbours estimator from 
# Heagerty2000 (DOI: https://doi.org/10.1111/j.0006-341x.2000.00337.x) in estimating the 
# aforementioned survivor functions. This function then culminates in producing both the associated
# ROC-graph, itself constructed using the trapezoidal rule from Mason2002 (DOI: https://doi.org/10.1256/003590002320603584),
# and the AUC-statistic in summarising the ROC-graph.
# Input:  [datGiven]: A validation dataset containing the variables of interest for testing prediction accuracy.
#         [cox]: A fitted cox model used to obtain marker values (theoretically either the risk scores exp(\beta.X)
#               or simply just the linear combination \beta.X.
#         [month_Start]: The prediction starting period of the time range over which prediction accuracy is tested.
#         [month_End], The last prediction period of the time range over which prediction accuracy is tested.
#         [Graph]: A boolean-valued toggle to produce the ROC-graph as a ggplot-object.
#         [lambda]: A smoothing parameter representing the %-valued neighborhood size,
#                  symmetrically calculated around each unique marker.
#         [numDigits]: The number of digits to which unique marker values are rounded, as an algorithmic efficiency boost
#         [method]: The estimation method by which True Positive Rates (TPR) and False Positive Rates (FPR) are calculated
#         [graphName]: The base name under which the produced ggplot graph will be saved in the given path directory
#         [genFigPath], A given path directory in which the ROC-graph (if produced) will be saved
#         [fld_ID]: An optional field name that designates whether to group certain observations together by subject/spell ID
#         [fld_Event]: A required field name that designates the main event indicator
#         [eventVal]: A required field that denotes the main event value against which [fld_Event] is tested
#         []:
# Output: [AUC]: The time-dependent Area under the curve (AUC) in summarising the corresponding time-dependent ROC-graph
#         [ROC_graph]: The associated ROC-graph as a ggplot-object
timedROC <- function(datGiven, cox, month_Start=0, month_End, lambda=0.05, method="NNE-0/1", numDigits=2, 
                     fld_ID=NA, fld_Event="MainEvent_Ind", eventVal=1, fld_StartTime="Start", fld_EndTime="Stop",
                     Graph=TRUE, graphName="timedROC-Graph", genFigPath=paste0(getwd(),"/")){
  
  # --- Preliminaries 
  # -- Testing Conditions
  datGiven = copy(dat); cox=coxExample; month_End=334; numDigits=2; Graph=TRUE; month_Start=0; method="NNE-0/1";lambda=0.05;
  fld_ID="ID"; fld_Event="Event_Ind"; fld_StartTime="Start"; fld_EndTime="End"; graphName="coxExample_cgd"; 
  eventVal=1; genFigPath=genFigPath
  
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
    stop("The arguments [fld_Event], [eventVal], [fld_StartTime], and [fld_EndTime] cannot be missing and must be specified. \n")
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
  
  # - Obtain various quantities towards implementing the Nearest Neighbour Estimator (NNE) method
  setorder(datGiven, cols = EndTime) #  Sort data set according to the given end time field
  # Get the corresponding rank order of the raw end times when sorted ascendingly
  vOrder <- order(datGiven$End)
  # Obtain unique end points (or ages) in the [datGiven] object, i.e., the "ordered failure times t_1 < t_2 < ... < t_m"
  uEnd <- datGiven$EndTime %>% unique %>% sort()
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
  # 2) he conditional survivor function S(t| M > c)
  
  # -- 1. Estimating the classical S(t) given each threshold
  # Implement the chosen estimator for S(t) and the choice of kernel (if Nearest Neighbour)
  
  if(method=="NNE-0/1"){ # Nearest Neighbour Estimator, using the 0/1 Kernel function from Akritas1994
    
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
    
  } else if(method=="NNE-Exp"){ # Method using exponential kernel function
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


# --- Unit test: timedROC()
# GoF_CoxSnell_graph(coxExample)

if(Test){
  
  # - Calculate AUC at median survival time for correctly-fitted Cox model | survivalROC (NNE)
  # NOTE: Uses the superior Nearest Neighbour Estimator (NNE) method for S(t)
  survivalROC(Stime=dat$End, status=dat$Event_Ind, entry=dat$Start, 
              method = "NNE",  span=0.05, predict.time=334,
              marker=round(predict(coxExample, type="lp"),2))
  ### RESULTS: Survival estimate at median survival time = 62% .. (should be 50%); AUC: 67.06%
  #            This already shows the bias of neglecting the ID-variable in clustering observations 
  #            around each relevant subject. 
  
  # - Calculate AUC at median survival time for wrongly-fitted Cox model | survivalROC (NNE)
  # NOTE: Uses the superior Nearest Neighbour Estimator (NNE) method for S(t)
  survivalROC(Stime=dat$End, status=dat$Event_Ind, entry=dat$Start, 
              method = "NNE", span=0.05, predict.time=334,
              marker=round(predict(coxExample_wrong, type="lp"),2))
  ### RESULTS: Survival estimate at median survival time remains 62% .. (should be 50%); AUC: 67.06%
  #           Therefore, prediction accuracy seems unaffected by model fitting mechanism
  
  # - Calculate AUC at median survival time for correctly-fitted Cox model | survivalROC (KM)
  # NOTE: Uses the classical but flawed Kaplan-Meier (KM) estimator for S(t)
  survivalROC(Stime=dat$End, status=dat$Event_Ind, entry=dat$Start, 
              method = "KM", predict.time=334,
              marker=round(predict(coxExample, type="lp"),2))
  ### RESULTS: Survival estimate at median survival time remains 62% .. (should be 50%); AUC: 69.93%
  
  # - Calculate AUC at median survival time for correctly-fitted Cox model | timeROC
  # NOTE: This function calculates the "Inverse Probability of Censoring Weighting [IPCW]" version
  # of the Cumulative/Dynamic (CD) time-dependence ROc-graph. Results may therefore differ inherently
  # NOTE2: Only implements the the KM-estimator for S(t) instead of the superior
  #   NNE-method from the survivalROC-package
  timeROC::timeROC(T=dat$End, delta=dat$Event_Ind, entry=dat$Start,
  )
  
  
  
  timedROC(datGiven=dat, cox=coxExample, month_End=334, method="NNE-0/11", numDigits=2, 
           fld_ID="ID", fld_Event="Event_Ind", eventVal=1, fld_StartTime="Start", fld_EndTime="End",
           graphName="coxExample_cgd", genFigPath=genFigPath)
  
  
  ### AB: Need to rewire the fields here, though I need an actual dataset to do that.
  timedROC(datGiven=datCredit_valid_TFD, cox=coxDelinq, month_End=12, lambda=0.05, method="NNE-0/1", numDigits=2, 
           fld_ID="PerfSpell_Key", fld_Event="Default_Ind", eventVal=1, fld_StartTime="Start", fld_EndTime="Stop",
           graphName="DefaultSurvModel-Cox", genFigPath=paste0(genFigPath, "TFD/tdROC/"))
  
}