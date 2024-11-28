# ============================== SURVIVAL FUNCTIONS ==============================
# Defining various custom functions relating to the estimation, analysis, and
# graphing of residuals, as part of testing the goodness-of-fit of a given 
# Cox regresion model
# --------------------------------------------------------------------------------
# PROJECT TITLE: Default Survival Modelling
# SCRIPT AUTHOR(S): Bernard Scheepers, Dr Arno Botha
# VERSION: 1.0 (November-2024)
# ================================================================================




# ----------------- 0. Dataset and models for conducting unit tests later --------
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




# ----------------- 1. Functions related to Cox-Snell residuals ------------------

# --- Function to calculate Cox-Snell residuals, as adjusted for censoring.
### AB: Insert article link that explains the rationale for this adjustment, e.g., "... as discussed by John2004 (DOI: ...)"
# Input: [cox], A fitted Cox proportional hazard model.
# Output: [cs], Cox-Snell residuals
CoxSnell_adjusted <- function(cox, dat){
  ### AB: This [Removed]-field seems hard-coded. I could not run GoF_CoxSnell_KS(coxExample) ...
  #       Though I suppose it is because you do not yet know how to program dynamically for given field names
  #       I will therefore move to the timedROC-function, make it dynamically programmable and trust that you
  #       will retro-apply the logic to this function, as well as possibly to GoF_CoxSnell_KS().
  cs <-  dat[Removed==1,Default_Ind] - 
    residuals(cox,type="martingale",collapse=dat$LoanID) +
    log(2)*(1 - dat[Removed==1,Default_Ind]) # Add log(2) to all observations that have a 0.
  return(cs)
}


# --- Function to calculate the Goodness-of-Fit (GoF) of a given Cox regression model, mainly 
# achieved by calculating the degree of similarity between Cox-Snell residuals
# and a random unit exponential distribution. The similarity degree is summarised by
# using the complement of the Kolmogorov-Smirnov test statistic (1-KS), which becomes
# our similarity measure; higher values = greater similarity = better fit.
# Input: [cox]: A fitted Cox proportional hazard model.
# Output: [Stat]: The test statistic value (1 - KS) as a measure of goodness-of-fit
#         [KS_graph]: A graph that combines the Cox-Snell empirical cumulative distribution
#                     with the unit exponential distribution function.
GoF_CoxSnell_KS <- function(cox, dat, GraphInd=T, legPos=c(0.5,0.5)) {
  
  # --- Preliminaries
  # - Obtain adjusted Cox-Snell residuals
  cs <- CoxSnell_adjusted(cox, dat)
  
  # - Initialize a unit exponential distribution
  exp <- rexp(length(cs),1)
  
  # - Perform the two-sample Kolmogorov-Smirnov test of distribution equality
  # H_0: cs and exp originates from the same distribution
  # NOTE: We only desire the KS test statistic in measuring distributional dissimilarity
  KS <- round(suppressWarnings(ks.test(cs,exp))$statistic,4)
  
  # Conditional execution towards creating a graphical output in accompanying main output (dissimilarity degree)
  if (GraphInd==T){
    # Calculate the empirical cumulative distribution of the obtained Cox-Snell residuals (cs)
    EmpDist <- ecdf(cs)
    
    # Create a grid of x-values for plotting purposes
    x <- sort(unique(c(cs, exp)))
    
    # Calculate CDF-values for each observation at each x value
    ### AB: From where is this function EmpDist?? I could not find it within any of the standad packages in script 0 ... 
    ### AB: As such, I'm ceasing my review of this function until we clarify this point.
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

### AB: I provided the following structure for you to complete, after addressing other comments
# --- Unit test: GoF_CoxSnell_KS()
# csResult <- GoF_CoxSnell_KS(coxExample,T)
# csResult$Stat;csResult$KS_graph
# ### RESULTS: D=0.1921

### AB: Rewrite function header according to the previous bits I wrote for you
# Function to graphically test the Cox-Snell residuals by plotting them against 
# their respective hazard rate. The line should tend towards the 45 degree line for
# a good fit.
# Input: cox - cox proportional hazard model
# Output: Graph - ggplot object to showcase the relationship
GoF_CoxSnell_graph <- function(cox){
  # Obtain adjusted Cox-Snell residuals
  cs <- CoxSnell_adjusted(cox)
  
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

### AB: I provided the following structure for you to complete, after addressing other comments
# --- Unit test: GoF_CoxSnell_graph()
# GoF_CoxSnell_graph(coxExample)
# 

### AB: Is the following still truly necessary? If so, then comment what this is all about
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




# ----------------- 2. Functions for estimating Schoenfeld residuals -------------

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



# --- House keeping
if (Test) {
  rm(cgd,cgd0,coxExample) 
}
