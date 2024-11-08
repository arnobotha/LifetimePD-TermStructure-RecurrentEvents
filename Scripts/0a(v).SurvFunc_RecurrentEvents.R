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

# ----------------- 1. Cox-Snell residuals analysis ---

# Function to calculate Cox-Snell residuals adjusted for censoring.
# Input: cox - Cox proportional hazard model.
# Output: cs - Cox-Snell residuals

cs_adjusted <- function(cox){
  cs <- get_csvec(cox) + log(2)*(1 - cox$y[, "status"]) # Add log(2) to all observations that have a 0.
  return(cs)
}

# Function to compute the Kolmogorov-Smirnov statistic (1-KS) for Cox-Snell 
# residuals as well as a ggplot graph to display it.
# Input: cox - Cox proportional hazard model.
# Output: KS_stat - 1 - Kolmogorov-Smirnov statistic
#         KS_graph -  Graph of the Cox-Snell empirical cumulative distribution
#                     function and the unit exponential distribution function.

cs_ks_test <- function(cox,GraphInd=T,legPos=c(0.5,0.5)) {
  # Obtain adjusted Cox-Snell residuals
  cs <- cs_adjusted(cox)
  
  # Initialize null distribution
  exp <- rexp(length(cs),1)
  
  # Perform the Kolmogorov-Smirnov test
  KS <- round(1 - ks.test(cs,exp)$statistic,4)
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
    
    # Plot the ECDFs with ggplot2
    (gg <- ggplot(datplot,aes(x=x,group=type)) + theme_minimal() + 
        theme(text = element_text(family="Cambria"), legend.position.inside=legPos,
              legend.position = "inside",
              legend.background = element_rect(fill="snow2", color="black", linetype="solid")) +
        labs(x = "x", y = "Cumulative Distribution Function") +
        stat_ecdf(aes(color=type,linetype=type)) + 
        geom_segment(data=segment_data,aes(x = x, xend = xend, y = y, yend = yend),
                     linetype = "dashed", color = "black") +
        annotate("label", x = x[D_location], y = (y1[D_location] + y2[D_location]) / 2,
                 label = paste("D =", percent(1-KS)), hjust = -0.1, vjust = -0.1, fill="white", alpha=0.6) +
        scale_color_manual(name = "Distributions", values = vCol, labels=vLabel) +
        scale_linetype_discrete(name = "Distributions",labels=vLabel) +
        scale_y_continuous(label=percent))
    # Prepare return object
    retOb <- list(KS_stat = 1-KS, KS_graph=gg)
  }else{
    retOb <- list(KS_stat = 1-KS)
  }
  
  return(retOb)
}

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

# Unit test
# cgd dataset: Data from a study on chronic granulomatous disease (CGD), focusing
# on repeated infections in patients.

# Load dataset
data(cgd)

# Fit a cox model
coxExample <- coxph(Surv(tstart,tstop,status) ~ sex + age + height + weight,cgd)

# Test cs_ks_test function
cs_ks_test(coxExample,T)
### RESULTS: D=0.1921

# Test cs_graph function
cs_graph(coxExample)

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
