###############################################################################
# Edited by Ariane Jong-Levinger
# Last Modified 8/13/26
###############################################################################

packages <- c("tidyverse", "ggplot2", "here", "ggpubr", "RColorBrewer", "ggtext")

package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)

########### Prep data for plotting #############################################

data <- read.csv(here::here("Data", "all_scatter_data.csv"))

# added 4/29/26
metadata <- read.csv(here::here("Data", "type_by_location_rz.csv"))

rain_zone_mapping <- read.csv(here::here("Data", "rain_zone_mapping.csv"))
N_to_num <- read.csv(here::here("Data", "N_to_num.csv"))

data <- merge(data, rain_zone_mapping, by = intersect(names(data), names(rain_zone_mapping)))
data <- merge(data, N_to_num, by = intersect(names(data), names(N_to_num)))

data$N <- as.factor(data$N)
data$N_numerical <- as.factor(data$N_numerical)
data$Quartile <- as.factor(data$Quartile)

#write.csv(data, here::here("Data", "scatter_data_with_rainzones.csv"))

# old version of data_summary using plyr
# data_summary <- function(data, varname, groupnames){
#   require(plyr)
#   summary_func <- function(x, col){
#     c(mean = mean(x[[col]], na.rm=TRUE),
#       sd = sd(x[[col]], na.rm=TRUE))
#   }
#   data_sum<-plyr::ddply(data, groupnames, .fun=summary_func,
#                   varname)
#   data_sum <- plyr::rename(data_sum, c("mean" = varname))
#   return(data_sum)
# }

# tidyverse version of data_summary
data_summary <- function(data, varname, groupnames) {
  data %>%
    group_by(across(all_of(groupnames))) %>%
    summarise(
      sd     = sd(.data[[varname]], na.rm = TRUE),
      n      = sum(!is.na(.data[[varname]])),
      "{varname}" := mean(.data[[varname]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      se         = sd / sqrt(n),
      moe_95 = qt(0.975, df = n - 1) * se # margin of error for 95% confidence interval
    ) %>%
    select(-se)
}

data$countBins <- cut(data$NumObservations, breaks = c(0, 50, 100, 150, 200, Inf), labels = c("20-50", "51-100", "101-150", "151-200", ">200"))
data$nearDetectBins <- cut(data$NearMedianDetectionLimit, breaks = c(-1, 0, 10, 20, 30, 40, 50, Inf), labels = c("None", "1-10", "11-20", "21-30", "31-40", "41-50", ">50"))
data$belowDetectBins <- cut(data$NumObservationsBelowDetect, breaks = c(-1, 0, 5, 10, 15, Inf), labels = c("None", "1-5", "6-10", "11-15", ">15"))
data$smallerCountBins <- cut(data$NumObservations, breaks = c(0, 30, 60, 90, Inf), labels = c("20-30", "31-60", "61-90", ">90"))
data$variationBins <- cut(data$coefficientOfVariation, breaks = c(0, .50, 1.00, 1.50, 2.00, Inf), labels = c("0-0.5", "0.51-1.00", "1.01-1.50", "1.51-2.00", ">2.00"))
data$numMonthsBins <- cut(data$averageAnnualMonthsSampled, breaks = c(0, 3, 6, 9, Inf), labels = c("1-3", "4-6", "7-9", ">9"))
data$inductionBins <- cut(data$inductionRatio, breaks = c(0, 100, 200, 300, 400, 500, Inf), labels = c("0-100", "101-200", "201-300", "301-400", "401-500", ">500"))
data$smallInductionBins <- cut(data$inductionRatio, breaks = c(0, 30, 60, 90, Inf), labels = c("0-30", "31-60", "61-90", ">90"))

### 6/3/26: Remove data with irrelevant BMP categories
  ## join BMP category to data
  # first, format metadata
  metadata_fmt <- metadata %>%
    rename(Flow = Flow_Type) %>%
    # remove improperly formatted Rain Zone column
    select(-Rain_Zone)
  
  # join properly formatted Rain Zone column
  metadata_rz_fix <- left_join(metadata_fmt, rain_zone_mapping, by = "Location")
  
  # join metadata to data
  data_w_md <- left_join(data, metadata_rz_fix, 
                         by = c("Location", "Analyte", "Flow", "Rain_Zone")) %>%
    # create long description for BMP Category code
    mutate(BMP_Category_long = case_when(
      BMP_Category == "BS" ~ "Grass Swale",
      BMP_Category == "PF" ~ "Permeable Friction",
      BMP_Category == "WC" ~ "Wetland Channel",
      BMP_Category == "WB" ~ "Wetland Basin",
      BMP_Category == "BR" ~ "Bioretention",
      BMP_Category == "CO" ~ "Composite",
      BMP_Category == "DB" ~ "Detention Basin",
      BMP_Category == "RP" ~ "Retention Pond",
      BMP_Category == "MF" ~ "Media Filter",
      BMP_Category == "MP" ~ "Maintenance Practice", # not sure we want to keep since all other BMPs are structural
      BMP_Category == "BI" ~ "Grass Strip",
      BMP_Category == "MD" ~ "Manufactured Device",
      BMP_Category == "PP" ~ "Porous Pavement",
      BMP_Category == "CX" ~ "Control", # documentation says "Control - No BMP/Control Site"
      BMP_Category == "IB" ~ "Infiltration Basin",
      BMP_Category == "GR" ~ "Green Roof",
      BMP_Category == "OT" ~ "Other",
      BMP_Category == "LD" ~ "LID",
      BMP_Category == "Not Applicable" ~ "N/A: NSWQ",
      .default = BMP_Category
    )) 
  
  # remove irrelevant BMP Categories: Composite, Control, LID, Maintenance Practice,   
  # and Other
  exclude_BMPCats <- c("CO", "CX", "LD", "MP", "OT")
  
  data_BMPCat_s <- data_w_md %>%
    filter(!BMP_Category %in% exclude_BMPCats) %>%
    # also exclude data from National Storm Water Database (updated 6/12/26)
    filter(BMP_Category != "Not Applicable") # n = 555 (# subsets * 5 quantiles)
  
  subset_stats <- data_BMPCat_s %>%
    select(N, Rain_Zone, Location, Analyte, Flow, BMP_Category, BMP_Category_long, NumObservations) %>%
    distinct() # 111
    
  # format data for use in statistical analysis of Coefficient of Variation (see
    # CompareCVBetweenGroups.R)
  cv_BMPType <- data_BMPCat_s %>%
    select(Rain_Zone, Location, Analyte, Flow, BMP_Category, BMP_Category_long,
           coefficientOfVariation, NumObservations) %>%
    distinct() 
  
  # save coeff. of variation data w/ BMP Type
  #write.csv(cv_BMPType, "Coefficient of Variation for Parent Datasets with Metadata.csv", row.names = FALSE)

########### Create sample size tables ##########################################
  # create table of sample sizes by rain zone, analyte, flow type
  nPar_RZ_Analyte_Flow <- cv_BMPType %>%
    group_by(Rain_Zone, Analyte, Flow) %>%
    summarize(N = n())
  #write.csv(nPar_RZ_Analyte_Flow, here::here("Data", "Table 1_Number of Parent Datasets by Rain Zone_Analyte_Flow Type.csv"), row.names = FALSE)
  
    # filter by N>=30
    # cv_BMPType_Ngte30 <- data_BMPCat_s %>%
    #   filter(NumObservations >= 30) %>%
    #   select(Rain_Zone, Location, Analyte, Flow, BMP_Category, BMP_Category_long,
    #         coefficientOfVariation, NumObservations) %>%
    #   distinct() 

    # nPar_RZ_Analyte_Flow_Ngte30 <- cv_BMPType_Ngte30 %>%
    #   group_by(Rain_Zone, Analyte, Flow) %>%
    #   summarize(N = n())
    #write.csv(nPar_RZ_Analyte_Flow_Ngte30, here::here("Data", "Table 1_Number of Parent Datasets by Rain Zone_Analyte_Flow Type_N_GTE_30.csv"), row.names = FALSE)

  # create table of sample sizes by BMP Category
  nPar_BMPType <- cv_BMPType %>%
    group_by(BMP_Category_long, BMP_Category) %>%
    summarise(n_parent_datasets = n(), n_observations = sum(NumObservations))
  
  # save as CSV
  #write.csv(nPar_BMPType, here::here("Data", "Table 2_Number of Parent Datasets by BMP Type.csv"), row.names = FALSE)
  
#diff_Ngte20_v_30 <- setdiff(nPar_RZ_Analyte_Flow, nPar_RZ_Analyte_Flow_Ngte30)

  ### remove subset rows from resampled dataset
  parent_only_d <- data_BMPCat_s %>%
    select(Rain_Zone, Location, Analyte, Flow, BMP_Category_long, NumObservations, smallerCountBins) %>%
    distinct() 
# 716 after filtering by CategoryAnalysisScreen_flag == "Y" OR "="

# save as CSV
#write.csv(parent_only_d, here::here("Data", "Parent Dataset Basic Metadata.csv"), row.names = FALSE)
  
  # format parent metadata
  metadata_fmt <- metadata %>%
    rename(Flow = Flow_Type) %>%
    select(-Rain_Zone)
  
  metadata_rz_fix <- left_join(metadata_fmt, rain_zone_mapping, by = "Location")
  
  parent_w_BMPCat <- left_join(parent_only_d, metadata_rz_fix, 
                     by = c("Location", "Analyte", "Flow", "Rain_Zone")) %>%
    # reorder columns
    select(Rain_Zone, BMP_Category, Location:smallerCountBins)
  
  # table: number of parent datasets by Rain Zone, BMP_Category, Analyte, and Flow type
  n_parent_summ_tbl <- parent_w_BMPCat %>% 
    group_by(Rain_Zone, BMP_Category, Analyte, Flow) %>%
    summarise(n_parent_datasets = n())
  
  # check # parent datasets
  n_pd <- sum(n_parent_summ_tbl$n_parent_datasets) # 716
  
  # check # parent datasets with N>20, N>25, N>30
  N_gt20_25 <- parent_w_BMPCat %>%
    mutate(N_cat = case_when(
      NumObservations > 30 ~ "N > 30",
      NumObservations > 25 & NumObservations <= 30 ~ "25 < N <= 30",
      NumObservations > 20 & NumObservations <= 25 ~ "20 < N <= 25",
      .default = "N = 20"
    )) %>%
    mutate(N_cat = factor(N_cat, levels = c("N = 20", "20 < N <= 25", "25 < N <= 30", "N > 30"))) %>%
    group_by(N_cat) %>%
    summarise(N = n()) %>%
    mutate(pct = N/n_pd*100)

  # table: count of parent datasets by categories for Number of Events in Parent Dataset
  n_parent_by_NumEventBin <- parent_only_d %>%
    group_by(Analyte, smallerCountBins) %>%
    summarise(n_parent = n())
  
    # save as CSV
    #write.csv(n_parent_by_NumEventBin, here::here("Data", "Number of Parent Datasets by Analyte and N.csv"), row.names = FALSE)

######## Make Final Figures ####################################################

####### Figure 1: COV by analyte, rain zone, flow type #########################
# add Rain Zone-Flow group column
  cv_rzf <- cv_BMPType %>%
    # make Rain Zone labels
    mutate(Rain_Zone_label = case_when(
      Rain_Zone == 1 ~ "1 (Great Lakes)",
      Rain_Zone == 2 ~ "2 (North East)",
      Rain_Zone == 3 ~ "3 (South East)",
      Rain_Zone == 4 ~ "4 (Lower Mississippi Valley)",
      Rain_Zone == 5 ~ "5 (Texas)",
      Rain_Zone == 6 ~ "6 (South West)",
      Rain_Zone == 7 ~ "7 (Northwest)",
      Rain_Zone == 8 ~ "8 (California)",
      Rain_Zone == 9 ~ "9 (Rocky Mountains)"
    )) %>%
    # convert Rain Zone to factor
    mutate(Rain_Zone_label = factor(Rain_Zone_label, levels = c("1 (Great Lakes)", "2 (North East)",
        "3 (South East)", "4 (Lower Mississippi Valley)", "5 (Texas)", "6 (South West)",
        "7 (Northwest)", "8 (California)", "9 (Rocky Mountains)")))

cdf_rzf_facet <- ggplot(cv_rzf, aes(x = coefficientOfVariation, color = Rain_Zone_label)) +
    stat_ecdf(geom = "step", linewidth = 1.5) +
    facet_grid(Flow ~ Analyte, scales = "free_y") +
    labs(x = "Coefficient of Variation",
         y = "CDF",
         color = "Rain Zone") +
    scale_color_brewer(palette = "Dark2", direction = -1) +
    #scale_color_brewer(palette = "Paired", direction = -1) +
    theme_bw(base_size = 41) +
    # get rid of gray strip
    theme(strip.placement = "outside", panel.spacing = unit(1, "lines"), strip.background = element_blank()) + 
    # bold panel heading, increase legend font size
    theme(strip.text = element_text(face = "bold"),
      legend.text = element_text(size = 40)
  )
  
  cdf_rzf_facet
  
  # save as png
  ggsave(here::here("Plots", "Final", "CDFs_by_Analyte_RainZone_FlowType_Dark2.png"), width = 80, height = 40, units = "cm")

######## Figure 2: COV CDFs by Analyte, Flow Type ##############################

# run Mann-Whitney test to compare Inflow & Outflow for each analyte
    
    # initialize empty output table
    analytes <- unique(cv_BMPType$Analyte)
    n_ana <- length(analytes)  
    flowtypes <- unique(cv_BMPType$Flow)
    n_flow <- length(flowtypes)

    mw_p_tbl <- tibble(
      Analyte = character(n_ana),
      MannWhit_pvalue = numeric(n_ana)
    )
      
    # loop over analytes
    for (i in 1:n_ana) {
      ana_s <- cv_BMPType %>%
        filter(Analyte == analytes[i])
      
      # subset by inflow/outflow
      inf_s <- ana_s %>%
        filter(Flow == "Inflow")
      
      out_s <- ana_s %>%
        filter(Flow == "Outflow")
      
      # run Mann-Whitney test
      mw_result <- wilcox.test(inf_s$coefficientOfVariation,
                               out_s$coefficientOfVariation,
                               paired = FALSE, exact = TRUE)
        
      # save results in output table
      mw_p_tbl$Analyte[i] <- analytes[i]
      mw_p_tbl$MannWhit_pvalue[i] <- mw_result$p.value
        
    }
    # for copper and TSS, p > 0.05, meaning no significant difference
    # for Phosphorus, p = 0.0002, so there IS a significant difference in medians of the distributions
  
# CDFs
    # format dataframe of Mann-Whitney results to use as labels
    mw_p_label_df <- mw_p_tbl %>%
      # format p values for label
      mutate(MW_p_fmt = case_when(
        Analyte == "TSS" ~ signif(MannWhit_pvalue,1),
        Analyte == "Phosphorus" ~ round(MannWhit_pvalue,4),
        Analyte == "Copper" ~ signif(MannWhit_pvalue,1)
      )) %>%
      # format p value text label
      mutate(p_label = paste0("'Mann-Whitney: '~italic(p)~'= ", MW_p_fmt, "'"))
    
    # plot empirical cdfs
    cdf_ana_facet <- ggplot(cv_BMPType, aes(x = coefficientOfVariation, linetype = Flow)) +
      stat_ecdf(geom = "step", linewidth = 1.5) +
      facet_wrap(~ Analyte, scales = "free_y") +
      # plot labels w/ p values (p<0.05 means significant difference)
      geom_text(data = mw_p_label_df,
                aes(x = 2.25, y = 0.3, label = p_label),
                parse = TRUE,
                #hjust = 1.3, vjust = 1.6,
                inherit.aes = FALSE,
                size = 11) +
      labs(x = "Coefficient of Variation",
           y = "CDF",
           linetype = "Flow Type") +
      #scale_color_brewer(palette = "Dark2") +
      theme_bw(base_size = 41) +
      # get rid of gray strip
      theme(strip.placement = "outside", panel.spacing = unit(1, "lines"), strip.background = element_blank()) + 
      # bold panel heading
      theme(strip.text = element_text(face = "bold"),
      legend.text = element_text(size = 40))
    
    cdf_ana_facet
  
  ggsave(here::here("Plots", "Final", "CDFs_by_Analyte_FlowType.png"), width = 80, height = 25, units = "cm", dpi = 350)
   

######## Figure 4: Mean RPD by parent and subset sample sizes ##################
    # and by Analyte, Flow type, and EMC percentile
  
  # first, identify which subsets of RPDs can be grouped together as 1 curve
  
  # create separate facets for Phosphorus:Inflow and Phosphorus:Outflow
  P_flow_s <- data_BMPCat_s %>%
    filter(Analyte == "Phosphorus") %>%
    # create facet group == Analyte + Flow
    mutate(facet_group = paste0(Analyte, ": ", Flow)) %>%
    # create npar_gt90 group == facet_group + whether smallerCountBins == ">90"
    mutate(npar_gt90 = case_when(
      smallerCountBins == ">90" ~ ">90",
      .default = "20-90"
    )) %>%
    # drop Flow column for joining with other data
    select(-Flow) %>%
    # format grouping vars
    mutate(N_numerical = factor(N_numerical, levels = c(5, 10, 15, 20, 25))) %>%
    mutate(Quartile = factor(Quartile, levels = c(0.1, 0.25, 0.5, 0.75, 0.9))) %>%
    select(N_numerical, Quartile, smallerCountBins, Analyte, facet_group, npar_gt90, NumObservationsBelowDetect, NumObservations, rpd, median)
  
  Cu_TSS_s <- data_BMPCat_s %>%
    filter(Analyte != "Phosphorus") %>%
    # create facet group == Analyte
    mutate(facet_group = Analyte) %>%
    # create npar_gt90 group == facet_group + whether smallerCountBins == ">90"
    mutate(npar_gt90 = case_when(
      smallerCountBins == ">90" ~ ">90",
      .default = "20-90"
    )) %>%
    # format grouping vars
    mutate(N_numerical = factor(N_numerical, levels = c(5, 10, 15, 20, 25))) %>%
    mutate(Quartile = factor(Quartile, levels = c(0.1, 0.25, 0.5, 0.75, 0.9))) %>%
    select(N_numerical, Quartile, smallerCountBins, Analyte, facet_group, npar_gt90, NumObservationsBelowDetect, NumObservations, rpd, median)
  
  data_BMPCat_s_fg <- rbind(P_flow_s, Cu_TSS_s) %>%
    # format facet group
    mutate(facet_group = factor(facet_group, levels = c("Copper", "Phosphorus: Inflow",
                                                        "Phosphorus: Outflow", "TSS"))) %>%
    # format npar_gt90 group
    mutate(npar_gt90 = factor(npar_gt90, levels = c(
      "20-90",
      ">90"
    )))

  # histograms of NumObservations by facet_group
  hist_N_fg <- ggplot(data_BMPCat_s_fg, aes(x = NumObservations)) +
    geom_histogram(aes(fill = facet_group)) +
    facet_wrap( ~ facet_group) +
    labs(x = "Parent Dataset Sample Size (N)", y = "Count") +
    theme_bw(base_size = 40)

  hist_N_fg

  # cdfs of NumObservations by facet_group
  cdf_N_fg <- ggplot(data_BMPCat_s_fg, aes(x = NumObservations, color = facet_group)) +
    stat_ecdf(geom = "step", linewidth = 1) +
    facet_wrap( ~ facet_group) +
    scale_color_brewer(palette = "Dark2") +
    labs(x = "Parent Dataset Sample Size (N)", y = "CDF") +
    theme_bw(base_size = 40) +
    theme(legend.title = element_blank())

  cdf_N_fg

  #ggsave(here::here("Plots", "CDFs_Parent_Dataset_Sample_Size_by_Analyte_Flow_Group.png"), width = 55, height = 60, units = "cm")

  # CDFs of parent dataset EMC grouped by smallerCountBins - facet by Quartile
    # format parent data
    parent_dataset_quartiles <- read.csv(here::here("Data", "parent_dataset_quartiles.csv"))
    parent_dataset_quartiles$smallerCountBins <- cut(parent_dataset_quartiles$NumSamples, breaks = c(0, 30, 60, 90, Inf), labels = c("20-30", "31-60", "61-90", ">90"))
    parent_dataset_quartiles <- parent_dataset_quartiles %>% mutate(
      Quartile = case_when(
        Quartile == 0.1 ~ "10th Percentile",
        Quartile == 0.25 ~ "25th Percentile",
        Quartile == 0.5 ~ "50th Percentile",
        Quartile == 0.75 ~ "75th Percentile",
        Quartile == 0.9 ~ "90th Percentile"
      ),
      Analyte_label = case_when(
        Analyte == "Copper" ~ "Copper (μg/L)",
        Analyte == "Phosphorus" ~ "Phosphorus (mg/L)",
        Analyte == "TSS" ~ "TSS (mg/L)"
      )
    )
    
    par_w_md <- left_join(parent_dataset_quartiles, metadata_rz_fix, 
                           by = c("Location", "Analyte", "Flow")) %>%
      # create long description for BMP Category code
      mutate(BMP_Category_long = case_when(
        BMP_Category == "BS" ~ "Grass Swale",
        BMP_Category == "PF" ~ "Permeable Friction",
        BMP_Category == "WC" ~ "Wetland Channel",
        BMP_Category == "WB" ~ "Wetland Basin",
        BMP_Category == "BR" ~ "Bioretention",
        BMP_Category == "CO" ~ "Composite",
        BMP_Category == "DB" ~ "Detention Basin",
        BMP_Category == "RP" ~ "Retention Pond",
        BMP_Category == "MF" ~ "Media Filter",
        BMP_Category == "MP" ~ "Maintenance Practice", # not sure we want to keep since all other BMPs are structural
        BMP_Category == "BI" ~ "Grass Strip",
        BMP_Category == "MD" ~ "Manufactured Device",
        BMP_Category == "PP" ~ "Porous Pavement",
        BMP_Category == "CX" ~ "Control", # documentation says "Control - No BMP/Control Site"
        BMP_Category == "IB" ~ "Infiltration Basin",
        BMP_Category == "GR" ~ "Green Roof",
        BMP_Category == "OT" ~ "Other",
        BMP_Category == "LD" ~ "LID",
        BMP_Category == "Not Applicable" ~ "N/A: NSWQ",
        .default = BMP_Category
      )) 
    
    # remove irrelevant BMP Categories: Composite, Control, LID, Maintenance Practice,   
    # and Other
    exclude_BMPCats <- c("CO", "CX", "LD", "MP", "OT")
    
    par_BMPCat_s <- par_w_md %>%
      filter(!BMP_Category %in% exclude_BMPCats) %>%
      # also exclude data from National Storm Water Database (updated 6/12/26)
      filter(BMP_Category != "Not Applicable") # n = 3615 
    
    # calculate EMC stats by Analyte, percentile, smallerCountBins
    stats_EMC_ana_pctl_N <- par_BMPCat_s %>%
      group_by(Analyte, Quartile, smallerCountBins) %>%
      summarise(mean_EMC = mean(EMC), med_EMC = median(EMC), min_EMC = min(EMC), 
        max_EMC = max(EMC), sd_EMC = sd(EMC), n = n())

    # calculate EMC stats by Analyte, smallerCountBins
    stats_EMC_ana_N <- par_BMPCat_s %>%
      group_by(Analyte, smallerCountBins) %>%
      summarise(mean_EMC = mean(EMC), med_EMC = median(EMC), min_EMC = min(EMC), 
        max_EMC = max(EMC), sd_EMC = sd(EMC), n = n())
  
  # for each analyte, make CDF by percentile and N_class        
  analytes <- unique(par_BMPCat_s$Analyte)
  n_ana <- length(analytes)
  ymax <- c(100, 600, 200)
  fn_prefix <- c("CDFs_EMC_by_Analyte_by_Parent_Dataset_Sample_Size")

  for (i in seq_along(analytes)) {
    # subset by analyte
    ana_s <- par_BMPCat_s %>%
      filter(Analyte == analytes[i])

    cdf_EMC_pctl_ana <- ggplot(ana_s, aes(x = EMC, color = smallerCountBins)) +
      stat_ecdf(geom = "step", linewidth = 1) +
      facet_wrap( ~ Quartile, ncol = 3) +
      scale_color_brewer(name = "N", palette = "Dark2") +
      labs(x = "EMC", y = "CDF", title = analytes[i]) +
      xlim(0, ymax[i]) +
      theme_bw(base_size = 40)
    
    print(cdf_EMC_pctl_ana)

    fn <- paste0(fn_prefix, "_", analytes[i], ".png")
    
    ggsave(here::here("Plots", fn), width = 66, height = 35, units = "cm")
  }

  # CDFs of EMC grouped by smallerCountBins - facet by Quartile
  bp_EMC_pctl_fg <- ggplot(par_BMPCat_s, aes(x = smallerCountBins, y = EMC, fill = smallerCountBins)) +
    geom_boxplot() +
    facet_grid(Quartile ~ Analyte_label) +
    scale_fill_brewer(palette = "Dark2") +
    labs(x = NULL, y = "EMC") +
    ylim(0,600)+
    theme_bw(base_size = 40) +
    theme(legend.title = element_blank())
  
  bp_EMC_pctl_fg #39 outliers removed
  
  #ggsave(here::here("Plots", "Boxplots_EMC_by_Analyte_by_Parent_Dataset_Sample_Size.png"), width = 70, height = 70, units = "cm")
  
  library(moments) # for calculating skew, kurtosis
  stats_EMC_ana_N <- par_BMPCat_s %>%
    rename(N_Class = smallerCountBins) %>%
    group_by(Analyte, N_Class) %>%
    summarise(med_EMC = median(EMC), sd_EMC = sd(EMC), skew_EMC = skewness(EMC), kurtosis_EMC = kurtosis(EMC), n = n())
  
  # save as CSV
 # write.csv(stats_EMC_ana_N, "EMC Stats by Analyte_Parent Sample Size.csv", row.names = FALSE)
   
 ### plot with all analytes
  
  # --- 1. Define which facet_group × Quartile combinations should show the median ---
  # calculate 95% CIs for each facet_group, smallerCountBins group, Quartile, and N_numerical,
    # then determine whether CIs for each SmallerCountBins group overlap across subset
    # sample size (N_numerical)
  fg_Nclass_pctl_CIoverlap <- data_BMPCat_s_fg %>%
    group_by(facet_group, smallerCountBins, Quartile, N_numerical) %>%
    summarise(
      sd  = sd(rpd, na.rm = TRUE),
      n = n(),
      rpd = mean(rpd, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      se         = sd / sqrt(n),
      lowr_95CI = rpd - qt(0.975, df = n - 1) * se, # lower limit of 95% confidence interval
      uppr_95CI = rpd + qt(0.975, df = n - 1) * se # upper limit of 95% confidence interval
    ) %>%
    select(-se) %>%
    # filter out datasets with n < 5
    filter(n >= 5)

  # summary table: for each facet_group x N_numerical x Quartile combo, check
  # whether the 95% CIs (lowr_95CI to uppr_95CI) overlap across the different
  # smallerCountBins groups at that subset sample size -- if they do at every
  # N_numerical value, that supports grouping data across smallerCountBins in
  # subsequent analyses
  fg_Nclass_pctl_CIoverlap_summary <- fg_Nclass_pctl_CIoverlap %>%
    group_by(facet_group, N_numerical, Quartile) %>%
    summarise(
      n_Nclass_levels    = n_distinct(smallerCountBins),
      max_lowr_95CI = max(lowr_95CI),
      min_uppr_95CI = min(uppr_95CI),
      # TRUE only when every smallerCountBins CI in this combo shares a common
      # overlapping region (i.e. the highest lower bound is still below the
      # lowest upper bound)
      CI_overlap_by_N = max_lowr_95CI <= min_uppr_95CI,
      .groups = "drop"
    ) %>%
    # collapse across N_numerical: CI_overlap_all_N counts how many of the n
    # N_numerical values had overlapping smallerCountBins CIs -- when
    # CI_overlap_all_N == n, the CIs overlapped across smallerCountBins at
    # every N_numerical value
    group_by(facet_group, Quartile) %>%
    summarise(CI_overlap_all_N = sum(CI_overlap_by_N), n = n())

  aggregate_pairs <- tribble(
    ~facet_group,          ~Quartile,
    "Copper", 0.1,
    "Copper", 0.25, 
    "Copper", 0.5,
    "Copper", 0.75,
    "Copper", 0.9,
    "Phosphorus: Inflow", 0.1,
    "Phosphorus: Inflow", 0.25,
    "Phosphorus: Inflow", 0.5,
    "Phosphorus: Inflow", 0.75,
    "Phosphorus: Inflow", 0.9,
    "Phosphorus: Outflow", 0.25,
    "Phosphorus: Outflow", 0.5
  )

  aggregate_pairs_gt90 <- tribble(
    ~facet_group,          ~Quartile,
    "Phosphorus: Outflow", 0.75  # added for plot using npar_gt90
  )
  
  # --- 2. Tag rows that belong to an aggregate panel ---
  data_tagged <- data_BMPCat_s %>%
    select(N_numerical, Quartile, smallerCountBins, Analyte, Flow, rpd) %>%
    # create facet group (Analyte + Flow for Phosphorus, Analyte for Copper & TSS)
    mutate(facet_group = case_when(
      Analyte == "Phosphorus" ~ paste0(Analyte, ": ", Flow),
      .default = Analyte
    )) %>%
    mutate(aggregate = interaction(facet_group, Quartile) %in%
             interaction(aggregate_pairs$facet_group, aggregate_pairs$Quartile)) %>%
    # format grouping vars so they plot in the correct order
    mutate(N_numerical = factor(N_numerical, levels = c(5, 10, 15, 20, 25))) %>%
    mutate(Quartile = factor(Quartile, levels = c(0.1, 0.25, 0.5, 0.75, 0.9))) %>%
    mutate(facet_group = factor(facet_group, levels = c("Copper", "Phosphorus: Inflow",
                                                        "Phosphorus: Outflow", "TSS")))
  
  data_tagged_gt90 <- data_BMPCat_s_fg %>%
    select(N_numerical, Quartile, smallerCountBins, Analyte, facet_group, npar_gt90, rpd) %>%
    mutate(aggregate = interaction(facet_group, Quartile) %in%
             interaction(aggregate_pairs_gt90$facet_group, aggregate_pairs_gt90$Quartile)) %>%
    # format grouping vars so they plot in the correct order
    mutate(N_numerical = factor(N_numerical, levels = c(5, 10, 15, 20, 25))) %>%
    mutate(Quartile = factor(Quartile, levels = c(0.1, 0.25, 0.5, 0.75, 0.9))) %>%
    mutate(facet_group = factor(facet_group, levels = c("Copper", "Phosphorus: Inflow",
                                                        "Phosphorus: Outflow", "TSS")))
  
  # --- 3. Subset by aggregate_pairs and calculate mean and 95% margin of error ---
  aggregated_data <- data_tagged %>%
    filter(aggregate) %>%
    group_by(N_numerical, Quartile, facet_group) %>%
    summarise(
      sd  = sd(rpd, na.rm = TRUE),
      n = n(),
      rpd = mean(rpd, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      se         = sd / sqrt(n),
      moe_95 = qt(0.975, df = n - 1) * se # margin of error for 95% confidence interval
    ) %>%
    select(-se) %>%
    mutate(smallerCountBins = "All bins") 
  
  aggregated_data_gt90 <- data_tagged_gt90 %>%
    filter(aggregate) %>%
    group_by(N_numerical, Quartile, facet_group, npar_gt90) %>%
    summarise(
      sd  = sd(rpd, na.rm = TRUE),
      n = n(),
      rpd = mean(rpd, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      se         = sd / sqrt(n),
      moe_95 = qt(0.975, df = n - 1) * se # margin of error for 95% confidence interval
    ) %>%
    select(-se) %>%
    mutate(smallerCountBins = npar_gt90) %>%
    # drop npar_gt90 column for joining
    select(-npar_gt90)
  
  # --- 4. Subset by non-aggregate panels and calculate mean and 95% margin of error ---
  pergroup_data <- data_tagged %>%
    filter(!aggregate) %>%
    # tagged gt90
    filter(!(facet_group == "Phosphorus: Outflow" & Quartile == "0.75")) %>%
    group_by(N_numerical, Quartile, facet_group, smallerCountBins) %>%
    summarise(
      sd  = sd(rpd, na.rm = TRUE),
      n = n(),
      rpd = mean(rpd, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      se         = sd / sqrt(n),
      moe_95 = qt(0.975, df = n - 1) * se # margin of error for 95% confidence interval
    ) %>%
    select(-se)
    
  # --- 5. Combine ---
  data_BMPCat_s_fg_Ngt50 <- data_BMPCat_s %>%
    select(Rain_Zone, Location, N_numerical, Quartile, smallerCountBins, Analyte, Flow, NumObservations, NumObservationsBelowDetect, rpd) %>%
    # create facet group (Analyte + Flow for Phosphorus, Analyte for Copper & TSS)
    mutate(facet_group = case_when(
      Analyte == "Phosphorus" ~ paste0(Analyte, ": ", Flow),
      .default = Analyte
    )) %>%
    mutate(N_numerical = factor(N_numerical, levels = c(5, 10, 15, 20, 25))) %>%
    mutate(Quartile = factor(Quartile, levels = c(0.1, 0.25, 0.5, 0.75, 0.9))) %>%
    mutate(facet_group = factor(facet_group, levels = c("Copper", "Phosphorus: Inflow",
                                                        "Phosphorus: Outflow", "TSS"))) %>%
    mutate(npar_gt50 = case_when(
      NumObservations <= 50 ~ "20-50",
      .default = "51-250"
    )) %>%
    mutate(npar_gt50 = factor(npar_gt50, levels = c("20-50", "51-250")))

  plot_data_by_N_gt50 <- data_BMPCat_s_fg_Ngt50 %>%
    group_by(N_numerical, Quartile, facet_group, npar_gt50) %>%
    summarise(
      sd  = sd(rpd, na.rm = TRUE),
      n = n(),
      rpd = mean(rpd, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      se         = sd / sqrt(n),
      moe_95 = qt(0.975, df = n - 1) * se # margin of error for 95% confidence interval
    ) %>%
    select(-se)

  # export table as CSV for further analysis
  fig4_data_out <- plot_data_by_N_gt50 %>%
    rename(n_subset = N_numerical, n_obs = n)

  #write.csv(fig4_data_out, here::here("Data", "Figure 4_Mean RPD Summary Statistics.csv"), row.names = FALSE)
  
  # --- 6. Plot ---
  
### mean RPD by parent dataset sample size and subset sample size (error bars = 95% CIs)
  #bin_levels <- c("20-30", "31-60", "61-90", ">90", "20-90", "All bins")

# plot_data <- bind_rows(pergroup_data, aggregated_data, aggregated_data_gt90) %>%
#     mutate(smallerCountBins = factor(smallerCountBins, levels = bin_levels))

 #bin_levels_gt90 <- c("20-30", "31-60", "61-90", ">90", "20-90")
  
  # plot_data_gt90 <- bind_rows(pergroup_data, aggregated_data_gt90) %>%
  #   mutate(smallerCountBins = factor(smallerCountBins, levels = bin_levels_gt90))
  
  #bin_levels_byNclass <- c("20-30", "31-60", "61-90", ">90")
 
  bin_levels_gt50 <- c("20-50", "51-250")

  set.seed(36)  # ensures reproducible jitter positions for geom_point

  # panel tags (a, b, c, ...), ordered left-to-right, top-to-bottom across the
  # Quartile (rows) x facet_group (columns) grid
  panel_labels_by_N_gt50 <- expand.grid(
    facet_group = c("Copper", "Phosphorus: Inflow", "Phosphorus: Outflow", "TSS"),
    Quartile    = c("10th Percentile", "25th Percentile", "50th Percentile", "75th Percentile", "90th Percentile")
  ) %>%
    mutate(label = letters[row_number()])

  plot_data_by_N_gt50 %>%
    mutate(Quartile = case_when(
      Quartile == 0.1 ~ "10th Percentile",
      Quartile == 0.25 ~ "25th Percentile",
      Quartile == 0.5 ~ "50th Percentile",
      Quartile == 0.75 ~ "75th Percentile",
      Quartile == 0.9 ~ "90th Percentile"
    )) %>%
    ggplot(aes(x = N_numerical, y = rpd, color = npar_gt50, group = npar_gt50)) +
    geom_line(stat = "smooth", method = "lm", se = FALSE, formula = y ~ I(log(x)), alpha = 0.7, size = 1.5) +
    geom_errorbar(aes(ymin = rpd - moe_95, ymax = rpd + moe_95), width = 0.2, size = 1.5,
                  position = position_jitter(width = 0.15, height = 0, seed = 38)) +
    geom_point(aes(shape = npar_gt50), size = 9, fill = NA, stroke = 2,
               position = position_jitter(width = 0.15, height = 0, seed = 38)) +
    scale_shape_manual(name = "Number of Events in Parent Dataset",
                       labels = bin_levels_gt50,
                       values = c(16, 18)) + # solid shapes: 17, 15, 4, 18,
    scale_size_manual(name = "Number of Events in Parent Dataset",
                       values = c(8, 11)) + #9, 9, 9, 11,
    facet_grid(Quartile~facet_group) +
    # panel letter labels
    geom_text(data = panel_labels_by_N_gt50,
              aes(x = Inf, y = Inf, label = label),
              inherit.aes = FALSE, hjust = 1.6, vjust = 1.5,
              fontface = "bold", size = 14) +
    scale_y_continuous(breaks = seq(0, 60, 20), minor_breaks = seq(0, 60, 10)) +
    coord_cartesian(ylim = c(0,60)) +
    xlab("Number of Monitored Events in Data Subset") +
    ylab("Mean Relative Percent Difference (%)") +
    theme_bw(base_size = 40) +  
    theme(strip.placement = "outside", panel.spacing = unit(1, "lines"), strip.background = element_blank()) + 
    theme(legend.position = "top", legend.box = "vertical") +
    guides(fill = "none") +
    scale_color_manual(name = "Number of Events in Parent Dataset",
                       labels = bin_levels_gt50,
                       values = c("#1D91C0", "#081D58")) + #"#87d069", "#41B6C4", "#225EA8", , "#000000"
    theme(legend.position = "top") +
    theme(panel.grid.minor = element_line(color = "grey75", linewidth = 0.7),
          panel.grid.major = element_line(color = "grey75", linewidth = 0.7),
          axis.text.x = element_text(angle = 45, hjust = 1),
          strip.text = element_text(face = "bold")
    )

  ggsave(here::here("Plots", "Final", "num_samples_analyte_PFlow_summary_aggregated_by_N_gt50.png"), width = 66, height = 70, units = "cm")
  
#display.brewer.all(colorblindFriendly = TRUE)
#brewer.pal(11, "RdBu")
#brewer.pal(9, "Greens")
#"#F7FBFF" "#DEEBF7" "#C6DBEF" "#9ECAE1" "#6BAED6" "#4292C6" "#2171B5" "#08519C" "#08306B"
#brewer.pal(9, "YlGnBu")
# "#FFFFD9" "#EDF8B1" "#C7E9B4" "#7FCDBB" "#41B6C4" "#1D91C0" "#225EA8" "#253494" "#081D58"

######## Figure 5: Mean RPD for 50th percentile for all analytes #####################################################

  # one rpd curve per facet_group, except TSS which is split into
  # "20-90" and ">90" groupings based on smallerCountBins (via npar_gt90)
  # "20-50" and "51-250" groupings based on smallerCountBins (via npar_gt90)
  fig5_data <- data_BMPCat_s_fg_Ngt50 %>%
    filter(Quartile == 0.5) %>%
    mutate(curve_group = case_when(
      facet_group == "TSS" & npar_gt50 == "51-250" ~ "TSS: 51-250 events",
      facet_group == "TSS" ~ "TSS: 20-50 events",
      .default = as.character(facet_group)
    )) %>%
    mutate(curve_group = factor(curve_group, levels = c(
      "Copper", "Phosphorus: Inflow", "Phosphorus: Outflow", "TSS: 20-50 events", "TSS: 51-250 events"
    ))) %>%
    group_by(N_numerical, curve_group) %>%
    summarise(
      sd  = sd(rpd, na.rm = TRUE),
      n = n(),
      rpd = mean(rpd, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      se         = sd / sqrt(n),
      moe_95 = qt(0.975, df = n - 1) * se # margin of error for 95% confidence interval
    ) %>%
    select(-se)

  # export table as CSV for further analysis
  fig5_data_out <- fig5_data %>%
    rename(n_subset = N_numerical, n_obs = n) %>%
    select(n_subset, curve_group, rpd, sd, moe_95, n_obs)

  #write.csv(fig5_data, here::here("Data", "Figure 4_50th Percentile RPD Summary Statistics.csv"), row.names = FALSE)

  set.seed(36)  # ensures reproducible jitter positions for geom_point

  fig5_data %>%
    ggplot(aes(x = N_numerical, y = rpd, color = curve_group, linetype = curve_group, group = curve_group)) + # 
    geom_line(stat = "smooth", method = "lm", se = FALSE, formula = y ~ I(log(x)), alpha = 0.7, size = 2) +
    geom_errorbar(aes(ymin = rpd - moe_95, ymax = rpd + moe_95), width = 0.2, size = 1.5, linetype = "solid") + #,
                  #position = position_jitter(width = 0.2, height = 0, seed = 38)
    geom_point(aes(shape = curve_group, size = curve_group), stroke = 2) + #, #alpha = 0.7, 
               #position = position_jitter(width = 0.2, height = 0, seed = 38)
    scale_shape_manual(name = "",
                       #values = c(16, 17, 15, 18, 4)) + # solid shapes
                       values = c(1, 2, 0, 5, 4)) +  # hollow shapes
    scale_linetype_manual(name = "",
                         values = c("solid", "solid", "solid", "dashed", "dashed")) +
    scale_size_manual(name = "",
                       values = c(11, 11, 11, 13, 11)) + # if using shape=18, use size=13
    guides(size = "none",
           shape = guide_legend(override.aes = list(size = 8), nrow = 2, byrow = TRUE),
           color = guide_legend(nrow = 2, byrow = TRUE),
           linetype = guide_legend(nrow = 2, byrow = TRUE)) +
    scale_color_manual(name = "",
                       values = c("#3F007D", "#41AB5D", "#006D2C", "#225EA8", "#081D58")) +
    scale_y_continuous(breaks = seq(0, 60, 10), minor_breaks = seq(0, 55, 50)) +
    coord_cartesian(ylim = c(0,45)) +
    xlab("Number of Monitored Events in Data Subset") +
    ylab("Mean RPD (%): 50th Percentile EMC") +
    theme_bw(base_size = 40) + theme(legend.position = "top", legend.box = "vertical") +
    guides(fill = "none") +
    theme(panel.grid.minor = element_line(color = "grey75", linewidth = 0.7),
          panel.grid.major = element_line(color = "grey75", linewidth = 0.7),
          #axis.text.x = element_text(angle = 45, hjust = 1)
    )

  ggsave(here::here("Plots", "Final", "rpd_50thpercentile_all_analytes_N_gt_50_hollow.png"), width = 42, height = 33.3, units = "cm", dpi = 350)


######## Figure 6: Mean RPD by % below detection limit ###############################################################
data_BMPCat_s_BDL <- data_BMPCat_s_fg_Ngt50 %>% mutate(
  percentBelowDetect = (NumObservationsBelowDetect/NumObservations) * 100
)

data_BMPCat_s_BDL$belowDetectPercentBins <- cut(data_BMPCat_s_BDL$percentBelowDetect, breaks = c(-Inf, 0, 9.99, 19.99, 39.99, Inf), labels = c("0%", "0.001-10%", "10-20%", "20-40%", ">40%"))

# make table: count of parent datasets by categories for "% of WQ results below detection limit"
n_parent_by_PctBlwDL <- data_BMPCat_s_BDL %>%
  select(Rain_Zone, Location, Analyte, Flow, facet_group, Quartile, belowDetectPercentBins) %>%
  distinct() %>% # 5 percentiles x 716 parent datasets = 3580
  # only select percentiles 10 and 25
  filter(Quartile == 0.1 | Quartile == 0.25) %>%
  group_by(Quartile, Analyte, belowDetectPercentBins) %>%
  summarise(n_parent = n())

# save as CSV
#write.csv(n_parent_by_PctBlwDL, here::here("Data", "Number of Parent Datasets by Below Detection Limit Plot Category.csv"), row.names = FALSE)

# query results for paper
rpd_by_PctBlwDL <- data_BMPCat_s_fg %>%
  mutate(percentBelowDetect = (NumObservationsBelowDetect/NumObservations) * 100) %>%
  mutate(belowDetectPercentBins = cut(percentBelowDetect,
        breaks = c(-Inf, 0, 9.99, 19.99, 39.99, Inf),
        labels = c("0%", "0.001-10%", "10-20%", "20-40%", ">40%")
        )) %>%
  select(N_numerical, facet_group, Quartile, belowDetectPercentBins, rpd) %>%
  #distinct() %>% # 5 percentiles x 723 parent datasets = 3615
  # only select percentiles 10 and 25
  filter(Quartile == 0.1 | Quartile == 0.25) %>%
  group_by(facet_group, Quartile, N_numerical, belowDetectPercentBins) %>%
  summarise(mean_rpd = mean(rpd), n = n())

rpd_diff_fm_0 <- rpd_by_PctBlwDL %>%
  group_by(facet_group, Quartile, N_numerical) %>%
  mutate(diff_from_zero = mean_rpd - mean_rpd[belowDetectPercentBins == "0%"]) %>%
  filter(belowDetectPercentBins == "10-20%")


# check data for issues
gt40_bdl_TP_TSS <- data_BMPCat_s_BDL %>%
  filter(belowDetectPercentBins == ">40%") %>%
  select(Rain_Zone, Location, Analyte, Flow, percentBelowDetect, belowDetectPercentBins, NumObservations) %>%
  distinct() %>%
  arrange(Analyte) # 2 parent datasets w/ Phosphorus w/ BDL >40%, N = 20 and N = 25 -> subsets for n=25 could not be created; only 1 mean RPD for subsets w/ n=20 

# check sample size of rpd stats used in BDL plot
bdl_rpd_check <- data_BMPCat_s_BDL %>%
  group_by(N_numerical, Quartile, belowDetectPercentBins, Analyte) %>%
  summarise(n = n(), mean_rpd = mean(rpd), sd_rpd = sd(rpd))

num_samples_below_detect_summary <- data_summary(data_BMPCat_s_BDL, varname = "rpd", groupnames = c("N_numerical", "Quartile", "belowDetectPercentBins", "Analyte"))
num_samples_below_detect_summary$N_numerical <- factor(num_samples_below_detect_summary$N_numerical, levels = c(5, 10, 15, 20, 25))
num_samples_below_detect_summary$Quartile <- factor(num_samples_below_detect_summary$Quartile, levels = c(0.1, 0.25, 0.5, 0.75, 0.9))
num_samples_below_detect_summary$belowDetectPercentBins <- factor(num_samples_below_detect_summary$belowDetectPercentBins, levels = c("0%", "0.001-10%", "10-20%", "20-40%", ">40%"))

plot_data_BDL <- num_samples_below_detect_summary %>%
  # filter only 10th and 25th percentile (other percentiles had low RPDs)
  filter(Quartile==0.1 | Quartile==0.25) %>%
  # filter by n>2
  filter(n > 2)
  # filter out Phosphorus and TSS where BDL >40%
  # filter(!((Analyte == "Phosphorus" | Analyte == "TSS") & belowDetectPercentBins == ">40%")) %>%
  # filter(!(Analyte == "TSS" & belowDetectPercentBins == "20-40%"))

# panel tags (a, b, c, ...), ordered left-to-right, top-to-bottom across the
  # Quartile (rows) x facet_group (columns) grid
  panel_labels_BDL <- expand.grid(
    Analyte = c("Copper", "Phosphorus", "TSS"),
    Quartile    = c("10th Percentile", "25th Percentile")
  ) %>%
    mutate(label = letters[row_number()])

set.seed(36)  # ensures reproducible jitter positions for geom_point

plot_data_BDL %>%
  mutate(Quartile = case_when(
    Quartile == 0.1 ~ "10th Percentile",
    Quartile == 0.25 ~ "25th Percentile",
    Quartile == 0.5 ~ "50th Percentile",
    Quartile == 0.75 ~ "75th Percentile",
    Quartile == 0.9 ~ "90th Percentile"
  )) %>%
  ggplot(aes(x = N_numerical, y = rpd, color = belowDetectPercentBins, group = belowDetectPercentBins)) + 
  geom_line(stat = "smooth", method = "lm", se = FALSE, formula = y ~ I(log(x)), size = 1.5, alpha = 0.7) + 
  geom_errorbar(aes(ymin = rpd - moe_95, ymax = rpd + moe_95), width = .2, size = 1.5,
                position = position_jitter(width = 0.1, height = 0, seed = 38)) +
  facet_grid(Quartile~Analyte) +
  #geom_hline(yintercept=20, linetype="dashed", color = "black", size=1) + 
  #geom_hline(yintercept=10, linetype="dashed", color = "black", size=1) + 
  geom_point(aes(shape = belowDetectPercentBins), size = 8, stroke = 2, fill = NA,
             position = position_jitter(width = 0.1, height = 0, seed = 38)) + 
  # panel letter labels
    geom_text(data = panel_labels_BDL,
              aes(x = Inf, y = Inf, label = label),
              inherit.aes = FALSE, hjust = 1.6, vjust = 1.5,
              fontface = "bold", size = 14) +
  scale_shape_manual(name = "Percentage of Water Quality \nResults Below Detection Limit",
                     labels = c("0%", "0.001-10%", "10-20%", "20-40%", ">40%"),
                     values = c(16, 17, 15, 4, 7)) +
  xlab("Number of Monitored Events in Data Subset") + ylab("Mean Relative Percent Difference (%)") +
  theme_bw(base_size = 40) +  theme(strip.placement = "outside", panel.spacing = unit(1, "lines"), strip.background = element_blank()) + theme(legend.position = "top", legend.box = "vertical") +
  #labs(color  = "Percentage of Water Quality \nResults Below Detection Limit", shape = "Percentage of Water Quality \nResults Below Detection Limit") +
  guides(fill = "none") +
  scale_color_manual(name = "Percentage of Water Quality \nResults Below Detection Limit",
                     labels = c("0%", "0.001-10%", "10-20%", "20-40%", ">40%"),
                     values = c("#87d069", "#41B6C4", "#2C7FB8","#253494", "#081D58")) +
  scale_y_continuous(breaks = seq(0, 120, 20)) + #, limits = c(-15,116)
  coord_cartesian(ylim = c(0,117)) +
  guides(legend = guide_legend(ncol = 2)) +
  theme(legend.position = "top") +
  theme(panel.grid.minor = element_blank(), 
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold")) # panel.grid.major = element_blank(), 

ggsave(here::here("Plots", "Final", "num_samples_percent_below_detect_summary_low_percentiles.png"), height = 40, width = 50, units = "cm")

