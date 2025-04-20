# ====================== ANALYTICS: HARRELL'S C-STATISTIC ========================
# Compare Harrell's c across variables for each technique
# --------------------------------------------------------------------------------
# PROJECT TITLE: Default Survival Modelling
# SCRIPT AUTHOR(S): Bernard Scheeprs (BS), Dr Arno Botha (AB)
# --------------------------------------------------------------------------------
# -- Script dependencies:
#   - 0.Setup.R
#   - 1.Data_Import.R
#   - 2a.Data_Prepare_Credit_Basic.R
#   - 2b.Data_Prepare_Credit_Advanced.R
#   - 2c.Data_Prepare_Credit_Advanced2.R
#   - 2d.Data_Enrich.R
#   - 2e.Data_Prepare_Macro.R
#   - 2f.Data_Fusion1.R
#   - 3a(i).Data_Transform.R
#   - 3c.Data_Fusion2.R
#   - 5a(i).InputSpace_TFD.R
#   - 5a(ii).InputSpace_AG.R
#   - 5a(iii).InputSpace_PWPST.R

# -- Inputs:
#   - Table_TFD | Prepared from script 3b
#   - datCredit_valid_TFD | Prepared from script 3b
#
# -- Outputs:
#   - <Analytics> | Harrell's c graph
# ================================================================================




# ----------------- 1. Load data

# - Confirm prepared datasets are loaded into memory
if (!exists('Table_TFD')) unpack.ffdf(paste0(genObjPath,"TFD_Univariate_Models"), tempPath);gc()
if (!exists('Table_PWPST')) unpack.ffdf(paste0(genObjPath,"PWPST_Univariate_Models"), tempPath);gc()




# ----------------- 2. Harrell's c-statistic comparative graph across variables and technique
# For each technique, Single factor models are built, whereupon Harrell's c-statistic is calculated 
# for each such a model, and subsequently graphed.

# - Combine metrics across techniques for graphing purposes
datGraph <- rbind(data.table(subset(Table_TFD, select=c("Variable", "Concordance")), Type="a_TFD"),
                  data.table(subset(Table_PWPST, select=c("Variable", "Concordance")), Type="c_PWP"))

# - Rename variables into more readable versions
vVarLabels <- c("g0_Delinq_SD_4", 
  "Arrears", 
  "g0_Delinq_Ave"="g0_Delinq_Avg",      
  "slc_acct_arr_dir_3_Change_Ind"="ArrearsDir_3_Changed", 
  "slc_acct_roll_ever_24_imputed_mean"="RollEver_24",
  "slc_acct_pre_lim_perc_imputed_med"="Prepaid_Pc", 
  "M_DTI_Growth"="M_DebtToIncome",
  "M_DTI_Growth_9"="M_DebtToIncome_9",
  "M_Inflation_Growth_6",
  "M_RealIncome_Growth",
  "M_Repo_Rate_6",
  "pmnt_method_grp"="PayMethod",
  "InterestRate_Nom"="InterestRate_Nominal",
  "Principal",
  "PerfSpell_Num"="Spell_Num",
  "AgeToTerm_Aggr_Mean"="AgeToTerm_Avg",
  "BalanceToPrincipal"
  )

# -- Graphing parameters
vCol <- brewer.pal(8, "Set1")
vlabel <- c("a_TFD"="TFD", "c_PWP"= "PWP")
chosenFont <- "Cambria"

# - Create main graph 
(g_comp_Acc <- ggplot(datGraph, aes(x = Variable, y = Concordance, fill=Type)) + theme_minimal() + 
    labs(x = "Variables",y = "Harrell's c", fill = "Time Definition") +
    theme(text=element_text(family=chosenFont),legend.position="bottom") +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +  # Adjust bar spacing
    scale_fill_manual(values = vCol, labels=vlabel) +  # Define colors for models
    coord_flip() +  # Flip the axes for horizontal bars
    scale_y_continuous(breaks=breaks_pretty(n=8), label=percent, limits=c(0,1)) + 
    scale_x_discrete(labels=vVarLabels)
  )


# - Save plot
dpi <- 220
ggsave(g_comp_Acc, file=paste0(genFigPath,"FULL SET/HarrellC_comp.png"),
       width=1200/dpi, height=1000/dpi,dpi=dpi, bg="white")

