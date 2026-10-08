###############################################################################
# Edited by Ariane Jong-Levinger
# Last Modified 7/17/26
###############################################################################

packages <- c("tidyverse", "ggplot2", "here", "ggpubr", "RColorBrewer")

package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)

# webshot::install_phantomjs()


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

write.csv(data, here::here("Data", "scatter_data_with_rainzones.csv"))

data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-plyr::ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- plyr::rename(data_sum, c("mean" = varname))
  return(data_sum)
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
      .default = BMP_Category
    )) 
  
  # remove irrelevant BMP Categories: Composite, Control, LID, Maintenance Practice,   
  # and Other
  exclude_BMPCats <- c("CO", "CX", "LD", "MP", "OT")
  
  data_BMPCat_s <- data_w_md %>%
    filter(!BMP_Category %in% exclude_BMPCats) %>%
    # also exclude data from National Storm Water Database (updated 6/12/26)
    filter(BMP_Category != "Not Applicable") # n = 16,585 (# subsets * 5 quantiles)
  
  subset_stats <- data_BMPCat_s %>%
    select(N, Rain_Zone, Location, Analyte, Flow, BMP_Category, BMP_Category_long, NumObservations) %>%
    distinct() # 3,317
    
  # format data for use in statistical analysis of Coefficient of Variation (see
    # CompareCVBetweenGroups.R)
  cv_BMPType <- data_BMPCat_s %>%
    select(Rain_Zone, Location, Analyte, Flow, BMP_Category, BMP_Category_long,
           coefficientOfVariation, NumObservations) %>%
    distinct() 
  
  # save coeff. of variation data w/ BMP Type
  #write.csv(cv_BMPType, "Coefficient of Variation for Parent Datasets with Metadata_20260612.csv", row.names = FALSE)

  # create table of sample sizes by BMP Category
  nPar_BMPType <- cv_BMPType %>%
    group_by(BMP_Category_long, BMP_Category) %>%
    summarise(n_parent_datasets = n(), n_observations = sum(NumObservations))
  
  # save as CSV
  #write.csv(nPar_BMPType, here::here("Data", "Number of Parent Datasets by BMP Type_20260616.csv"), row.names = FALSE)
  
    # Analyte & BMP Category
    nPar_BMPType_Analyte <- cv_BMPType %>%
      group_by(BMP_Category_long, BMP_Category, Analyte) %>%
      summarise(n_parent_datasets = n(), n_observations = sum(NumObservations))
  
    # Flow Type & Analyte & BMP Category
    nPar_BMPType_Analyte_Flow <- cv_BMPType %>%
      group_by(BMP_Category_long, BMP_Category, Analyte, Flow) %>%
      summarise(n_parent_datasets = n(), n_observations = sum(NumObservations))
  
# 4/24/26: create inventory of parent datasets (updated 6/3/26 after removing 
  #irrelevant BMP Categories)
  
  ## remove subset rows
  parent_only_d <- data_BMPCat_s %>%
    select(Rain_Zone, Location, Analyte, Flow, NumObservations, smallerCountBins) %>%
    distinct() # originally 834 parent datasets, 723 after removing bad BMP Cats
  
  # format parent metadata
  metadata_fmt <- metadata %>%
    rename(Flow = Flow_Type) %>%
    select(-Rain_Zone)
  
  metadata_rz_fix <- left_join(metadata_fmt, rain_zone_mapping, by = "Location")
  
  parent_w_BMPCat <- left_join(parent_only_d, metadata_rz_fix, 
                     by = c("Location", "Analyte", "Flow", "Rain_Zone")) %>%
    # reorder columns
    select(Rain_Zone, BMP_Category, Location:smallerCountBins)

  
  parent_w_BMPCat_long <- left_join(parent_only_d, metadata_rz_fix, 
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
      .default = BMP_Category
    )) %>%
    # select and reorder columns
    select(Rain_Zone, Location, BMP_Category_long, Analyte:NumObservations)
  
  # save as CSV
  #write.csv(parent_w_BMPCat_long, here::here("Data", "Parent Dataset Metadata for Supp Materials_20260626.csv"), row.names = FALSE)
  
  # table: number of parent datasets by Rain Zone, BMP_Category, Analyte, and Flow type
  n_parent_summ_tbl <- parent_w_BMPCat %>% 
    group_by(Rain_Zone, BMP_Category, Analyte, Flow) %>%
    summarise(n_parent_datasets = n())
  
    # save as CSV
    #write.csv(n_parent_summ_tbl, here::here("Data", "Number of Parent Datasets by RZ_BMPType_Analyte_FlowType_20260612.csv"), row.names = FALSE)
  
  # check # parent datasets 
  sum(n_parent_summ_tbl$n_parent_datasets) 

  # table: number of observations by Location, BMP_Category, Analyte, and Flow type
  parent_w_BMPCat_fmt <- parent_w_BMPCat %>%
    arrange(Rain_Zone, BMP_Category, Analyte, Flow)
  
    # save as CSV
    #write.csv(parent_w_BMPCat_fmt, here::here("Data", "Number of Observations per Parent Dataset_with BMPType_20260612.csv"), row.names = FALSE)
  
  # table: count of parent datasets by categories for Number of Events in Parent Dataset
  n_parent_by_NumEventBin <- parent_only_d %>%
    group_by(Analyte, smallerCountBins) %>%
    summarise(n_parent = n())
  
    # save as CSV
    #write.csv(n_parent_by_NumEventBin, here::here("Data", "Number of Parent Datasets by Plot Category_20260619.csv"), row.names = FALSE)

# updated 6/12/26
  # make a dataframe of metadata for parent datasets for Phosphorus & TSS with
    # 31-60 sampled events
  lg_nevents_P_TSS <- parent_only_d %>%
    filter(Analyte != "Copper") %>%
    filter(smallerCountBins == "31-60") %>%
    arrange(Rain_Zone, Analyte, Flow)
  
  # save as CSV
  #write.csv(lg_nevents_P_TSS, here::here("Data", "Metadata for P and TSS Parent Datasets w 31-60 Events_20260612.csv"), row.names = FALSE)
  
#Range summary without analyte facet
range_summary <- data_summary(data, varname = "rpd", groupnames = c("N_numerical", "Quartile", "inductionBins"))
range_summary$N_numerical <- factor(range_summary$N_numerical, levels = c(5, 10, 15, 20, 25))
range_summary$Quartile <- factor(range_summary$Quartile, levels = c(0.1, 0.25, 0.5, 0.75, 0.9))

range_summary %>%
  ggplot(aes(x = N_numerical, y = rpd, color = Quartile, group = Quartile)) + 
  geom_point() + 
  stat_smooth(method = "lm", se = FALSE, formula = y ~ I(log(x))) + 
  geom_errorbar(aes(ymin = rpd - sd, ymax = rpd + sd), width = .2) +
  facet_grid(.~inductionBins) +
  xlab("Random Sample Size") + ylab("Mean Relative Percent Difference (%)") +
  theme_bw(base_size = 20) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(here::here("Plots", "range_summary.jpg"), width = 30, height = 20, units = "cm")

#Range summary with analyte facet
range_analyte_summary <- data_summary(data, varname = "rpd", groupnames = c("N_numerical", "Quartile", "inductionBins", "Analyte"))
range_analyte_summary$N_numerical <- factor(range_analyte_summary$N_numerical, levels = c(5, 10, 15, 20, 25))
range_analyte_summary$Quartile <- factor(range_analyte_summary$Quartile, levels = c(0.1, 0.25, 0.5, 0.75, 0.9))

range_analyte_summary %>%
  ggplot(aes(x = N_numerical, y = rpd, color = Quartile, group = Quartile)) + 
  geom_point() + 
  stat_smooth(method = "lm", se = FALSE, formula = y ~ I(log(x))) + 
  geom_errorbar(aes(ymin = rpd - sd, ymax = rpd + sd), width = .2) +
  facet_grid(Analyte~inductionBins) +
  xlab("Random Sample Size") + ylab("Mean Relative Percent Difference (%)") +
  theme_bw(base_size = 20) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(here::here("Plots", "range_and_analyte_summary.jpg"), width = 30, height = 20, units = "cm")

#Coefficient of varitation summary without analyte facet
cv_summary <- data_summary(data, varname = "rpd", groupnames = c("N_numerical", "Quartile", "variationBins"))
cv_summary$N_numerical <- factor(cv_summary$N_numerical, levels = c(5, 10, 15, 20, 25))
cv_summary$Quartile <- factor(cv_summary$Quartile, levels = c(0.1, 0.25, 0.5, 0.75, 0.9))

cv_summary %>%
  ggplot(aes(x = N_numerical, y = rpd, color = Quartile, group = Quartile)) + 
  geom_point() + 
  stat_smooth(method = "lm", se = FALSE, formula = y ~ I(log(x))) + 
  geom_errorbar(aes(ymin = rpd - sd, ymax = rpd + sd), width = .2) +
  facet_grid(.~variationBins) +
  xlab("Random Sample Size") + ylab("Mean Relative Percent Difference (%)") +
  theme_bw(base_size = 20) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(here::here("Plots", "cv_summary.jpg"), width = 30, height = 20, units = "cm")

#Coefficient of variation summary with analyte facet
cv_analyte_summary <- data_summary(data, varname = "rpd", groupnames = c("N_numerical", "Quartile", "variationBins", "Analyte"))
cv_analyte_summary$N_numerical <- factor(cv_analyte_summary$N_numerical, levels = c(5, 10, 15, 20, 25))
cv_analyte_summary$Quartile <- factor(cv_analyte_summary$Quartile, levels = c(0.1, 0.25, 0.5, 0.75, 0.9))

cv_analyte_summary %>%
  ggplot(aes(x = N_numerical, y = rpd, color = Quartile, group = Quartile)) + 
  geom_point() + 
  stat_smooth(method = "lm", se = FALSE, formula = y ~ I(log(x))) + 
  geom_errorbar(aes(ymin = rpd - sd, ymax = rpd + sd), width = .2) +
  facet_grid(Analyte~variationBins) +
  xlab("Random Sample Size") + ylab("Mean Relative Percent Difference (%)") +
  theme_bw(base_size = 20) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(here::here("Plots", "cv_and_analyte_summary.jpg"), width = 30, height = 20, units = "cm")

#Number of samples summary without analyte facet
num_samples_summary <- data_summary(data, varname = "rpd", groupnames = c("N_numerical", "Quartile", "smallerCountBins"))
num_samples_summary$N_numerical <- factor(num_samples_summary$N_numerical, levels = c(5, 10, 15, 20, 25))
num_samples_summary$Quartile <- factor(num_samples_summary$Quartile, levels = c(0.1, 0.25, 0.5, 0.75, 0.9))

num_samples_summary %>%
  ggplot(aes(x = N_numerical, y = rpd, color = Quartile, group = Quartile)) + 
  geom_point() + 
  stat_smooth(method = "lm", se = FALSE, formula = y ~ I(log(x))) + 
  geom_errorbar(aes(ymin = rpd - sd, ymax = rpd + sd), width = .2) +
  facet_grid(.~smallerCountBins) +
  xlab("Random Sample Size") + ylab("Mean Relative Percent Difference (%)") +
  theme_bw(base_size = 20) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(here::here("Plots", "num_samples_summary.jpg"), width = 30, height = 20, units = "cm")

#Number of samples summary with analyte facet
num_samples_analyte_summary <- data_summary(data, varname = "rpd", groupnames = c("N_numerical", "Quartile", "smallerCountBins", "Analyte"))
num_samples_analyte_summary$N_numerical <- factor(num_samples_analyte_summary$N_numerical, levels = c(5, 10, 15, 20, 25))
num_samples_analyte_summary$Quartile <- factor(num_samples_analyte_summary$Quartile, levels = c(0.1, 0.25, 0.5, 0.75, 0.9))

num_samples_analyte_summary %>%
  ggplot(aes(x = N_numerical, y = rpd, color = Quartile, group = Quartile)) + 
  geom_point() + 
  stat_smooth(method = "lm", se = FALSE, formula = y ~ I(log(x))) + 
  geom_errorbar(aes(ymin = rpd - sd, ymax = rpd + sd), width = .2) +
  facet_grid(Analyte~smallerCountBins) +
  xlab("Random Sample Size") + ylab("Mean Relative Percent Difference (%)") +
  theme_bw(base_size = 20) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) 
ggsave(here::here("Plots", "num_samples_and_analyte_summary.jpg"), width = 30, height = 20, units = "cm")

#Below detection summary without analyte facet
below_detection_summary <- data_summary(data, varname = "rpd", groupnames = c("N_numerical", "Quartile", "belowDetectBins"))
below_detection_summary$N_numerical <- factor(below_detection_summary$N_numerical, levels = c(5, 10, 15, 20, 25))
below_detection_summary$Quartile <- factor(below_detection_summary$Quartile, levels = c(0.1, 0.25, 0.5, 0.75, 0.9))

below_detection_summary %>%
  ggplot(aes(x = N_numerical, y = rpd, color = Quartile, group = Quartile)) + 
  geom_point() + 
  stat_smooth(method = "lm", se = FALSE, formula = y ~ I(log(x))) + 
  geom_errorbar(aes(ymin = rpd - sd, ymax = rpd + sd), width = .2) +
  facet_grid(.~belowDetectBins) +
  xlab("Random Sample Size") + ylab("Mean Relative Percent Difference (%)") +
  theme_bw(base_size = 20) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(here::here("Plots", "below_detection_summary.jpg"), width = 30, height = 20, units = "cm")

#Below detection summary with analyte facet
below_detection_analyte_summary <- data_summary(data, varname = "rpd", groupnames = c("N_numerical", "Quartile", "belowDetectBins", "Analyte"))
below_detection_analyte_summary$N_numerical <- factor(below_detection_analyte_summary$N_numerical, levels = c(5, 10, 15, 20, 25))
below_detection_analyte_summary$Quartile <- factor(below_detection_analyte_summary$Quartile, levels = c(0.1, 0.25, 0.5, 0.75, 0.9))

below_detection_analyte_summary %>%
  ggplot(aes(x = N_numerical, y = rpd, color = Quartile, group = Quartile)) + 
  geom_point() + 
  stat_smooth(method = "lm", se = FALSE, formula = y ~ I(log(x))) + 
  geom_errorbar(aes(ymin = rpd - sd, ymax = rpd + sd), width = .2) +
  facet_grid(Analyte~belowDetectBins) +
  xlab("Random Sample Size") + ylab("Mean Relative Percent Difference (%)") +
  theme_bw(base_size = 20) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(here::here("Plots", "below_detection_analyte_summary.jpg"), width = 30, height = 20, units = "cm")

#Months summary with analyte facet
months_analyte_summary <- data_summary(data, varname = "rpd", groupnames = c("N_numerical", "Quartile", "numMonthsBins", "Analyte"))
months_analyte_summary$N_numerical <- factor(months_analyte_summary$N_numerical, levels = c(5, 10, 15, 20, 25))
months_analyte_summary$Quartile <- factor(months_analyte_summary$Quartile, levels = c(0.1, 0.25, 0.5, 0.75, 0.9))

months_analyte_summary %>%
  ggplot(aes(x = N_numerical, y = rpd, color = Quartile, group = Quartile)) + 
  geom_point() + 
  stat_smooth(method = "lm", se = FALSE, formula = y ~ I(log(x))) + 
  geom_errorbar(aes(ymin = rpd - sd, ymax = rpd + sd), width = .2) +
  facet_grid(Analyte~numMonthsBins) +
  xlab("Random Sample Size") + ylab("Mean Relative Percent Difference (%)") +
  theme_bw(base_size = 20) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(here::here("Plots", "months_analyte_summary.jpg"), width = 30, height = 20, units = "cm")


#Rain zone summary with analyte facet
rainzone_analyte_summary <- data_summary(data, varname = "rpd", groupnames = c("N_numerical", "Quartile", "Rain_Zone", "Analyte"))
rainzone_analyte_summary$N_numerical <- factor(rainzone_analyte_summary$N_numerical, levels = c(5, 10, 15, 20, 25))
rainzone_analyte_summary$Quartile <- factor(rainzone_analyte_summary$Quartile, levels = c(0.1, 0.25, 0.5, 0.75, 0.9))

rainzone_analyte_summary %>%
  ggplot(aes(x = N_numerical, y = rpd, color = Quartile, group = Quartile)) + 
  geom_point() + 
  stat_smooth(method = "lm", se = FALSE, formula = y ~ I(log(x))) + 
  geom_errorbar(aes(ymin = rpd - sd, ymax = rpd + sd), width = .2) +
  facet_grid(Analyte~Rain_Zone) +
  xlab("Random Sample Size") + ylab("Mean Relative Percent Difference (%)") +
  theme_bw(base_size = 20) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(here::here("Plots", "rainzone_analyte_summary.jpg"), width = 30, height = 20, units = "cm")


#Bonus scatter plot that shows an improvement in RPD with increasing sample size for the 25, 50, 75 quartiles when looking at num samples:scatter
data %>% 
  ggplot(aes(x = NumObservations, y = inductionRatio)) +
  geom_point(aes(color = rpd)) +
  facet_grid(Quartile~N_numerical) + 
  theme_bw(base_size = 20) +
  scale_color_gradient(low = "green", high = "red") +
  xlab("Number of Observations") + 
  ylab("Spread of Values") +
  labs(color = "Relative\nPercent\nDifference") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_smooth(method = "lm")
  


rainzone_analyte_summary %>%
  mutate(Quartile = case_when(
    Quartile == 0.1 ~ "10th Percentile",
    Quartile == 0.25 ~ "25th Percentile",
    Quartile == 0.5 ~ "50th Percentile",
    Quartile == 0.75 ~ "75th Percentile",
    Quartile == 0.9 ~ "90th Percentile"
  )) %>%
  ggplot(aes(x = N_numerical, y = rpd, color = as.character(Rain_Zone), group = as.character(Rain_Zone))) + 
  geom_point(aes(shape = as.character(Rain_Zone)), size = 3, fill = NA) + 
  #scale_shape_manual(values = c(0, 1, 2, 5, 6)) +
  stat_smooth(method = "lm", se = FALSE, formula = y ~ I(log(x))) + 
  geom_errorbar(aes(ymin = rpd - sd, ymax = rpd + sd), width = .2) +
  facet_grid(Analyte~Quartile) +
  geom_hline(yintercept=20, linetype="dashed", color = "black", size=1) + 
  geom_hline(yintercept=10, linetype="dashed", color = "black", size=1) + 
  xlab("Number of Monitored Events in Data Subset") + ylab("Mean Relative Percent Difference (%)") +
  theme_bw(base_size = 24) +
  labs(color  = "EPA Rain Zone", shape = "EPA Rain Zone") +
  guides(fill = "none") +
  theme(legend.position = "top") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(here::here("Plots", "rainzone_analyte_summary.jpg"), width = 40, height = 30, units = "cm")



#Number of samples summary with analyte facet
num_samples_analyte_flow_summary <- data_summary(data, varname = "rpd", groupnames = c("N_numerical", "Quartile", "smallerCountBins", "Analyte", "Flow"))
num_samples_analyte_flow_summary$N_numerical <- factor(num_samples_analyte_flow_summary$N_numerical, levels = c(5, 10, 15, 20, 25))
num_samples_analyte_flow_summary$Quartile <- factor(num_samples_analyte_flow_summary$Quartile, levels = c(0.1, 0.25, 0.5, 0.75, 0.9))

num_samples_analyte_flow_summary %>%
  mutate(Quartile = case_when(
    Quartile == 0.1 ~ "10th Percentile",
    Quartile == 0.25 ~ "25th Percentile",
    Quartile == 0.5 ~ "50th Percentile",
    Quartile == 0.75 ~ "75th Percentile",
    Quartile == 0.9 ~ "90th Percentile"
  )) %>% filter(Analyte == "TSS") %>%
  ggplot(aes(x = N_numerical, y = rpd, color = smallerCountBins, group = smallerCountBins)) + 
  geom_point(aes(shape = smallerCountBins), size = 3, fill = NA) + 
  #scale_shape_manual(values = c(0, 1, 2, 5, 6)) +
  stat_smooth(method = "lm", se = FALSE, formula = y ~ I(log(x))) + 
  geom_errorbar(aes(ymin = rpd - sd, ymax = rpd + sd), width = .2) +
  facet_grid(Flow~Quartile) +
  geom_hline(yintercept=20, linetype="dashed", color = "black", size=1) + 
  geom_hline(yintercept=10, linetype="dashed", color = "black", size=1) + 
  xlab("Number of Monitored Events in Data Subset") + ylab("Mean Relative Percent Difference (%)") +
  theme_bw(base_size = 24) +
  labs(color  = "Number of Events in Parent Dataset", shape = "Number of Events in Parent Dataset") +
  guides(fill = "none") +
  theme(legend.position = "top") + ggtitle("TSS") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(here::here("Plots", "num_samples_analyte_flow_summary_TSS.jpg"), width = 40, height = 30, units = "cm")

num_samples_analyte_flow_summary %>%
  mutate(Quartile = case_when(
    Quartile == 0.1 ~ "10th Percentile",
    Quartile == 0.25 ~ "25th Percentile",
    Quartile == 0.5 ~ "50th Percentile",
    Quartile == 0.75 ~ "75th Percentile",
    Quartile == 0.9 ~ "90th Percentile"
  )) %>% filter(Analyte == "Copper") %>%
  ggplot(aes(x = N_numerical, y = rpd, color = smallerCountBins, group = smallerCountBins)) + 
  geom_point(aes(shape = smallerCountBins), size = 3, fill = NA) + 
  #scale_shape_manual(values = c(0, 1, 2, 5, 6)) +
  stat_smooth(method = "lm", se = FALSE, formula = y ~ I(log(x))) + 
  geom_errorbar(aes(ymin = rpd - sd, ymax = rpd + sd), width = .2) +
  facet_grid(Flow~Quartile) +
  geom_hline(yintercept=20, linetype="dashed", color = "black", size=1) + 
  geom_hline(yintercept=10, linetype="dashed", color = "black", size=1) + 
  xlab("Number of Monitored Events in Data Subset") + ylab("Mean Relative Percent Difference (%)") +
  theme_bw(base_size = 24) +
  labs(color  = "Number of Events in Parent Dataset", shape = "Number of Events in Parent Dataset") +
  guides(fill = "none") +
  theme(legend.position = "top") + ggtitle("Copper") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(here::here("Plots", "num_samples_analyte_flow_summary_Copper.jpg"), width = 40, height = 30, units = "cm")

num_samples_analyte_flow_summary %>%
  mutate(Quartile = case_when(
    Quartile == 0.1 ~ "10th Percentile",
    Quartile == 0.25 ~ "25th Percentile",
    Quartile == 0.5 ~ "50th Percentile",
    Quartile == 0.75 ~ "75th Percentile",
    Quartile == 0.9 ~ "90th Percentile"
  )) %>% filter(Analyte == "Phosphorus") %>%
  ggplot(aes(x = N_numerical, y = rpd, color = smallerCountBins, group = smallerCountBins)) + 
  geom_point(aes(shape = smallerCountBins), size = 3, fill = NA) + 
  #scale_shape_manual(values = c(0, 1, 2, 5, 6)) +
  stat_smooth(method = "lm", se = FALSE, formula = y ~ I(log(x))) + 
  geom_errorbar(aes(ymin = rpd - sd, ymax = rpd + sd), width = .2) +
  facet_grid(Flow~Quartile) +
  geom_hline(yintercept=20, linetype="dashed", color = "black", size=1) + 
  geom_hline(yintercept=10, linetype="dashed", color = "black", size=1) + 
  xlab("Number of Monitored Events in Data Subset") + ylab("Mean Relative Percent Difference (%)") +
  theme_bw(base_size = 24) +
  labs(color  = "Number of Events in Parent Dataset", shape = "Number of Events in Parent Dataset") +
  guides(fill = "none") +
  theme(legend.position = "top") + ggtitle("Phosphorus") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) 
ggsave(here::here("Plots", "num_samples_analyte_flow_summary_Phosphorus.jpg"), width = 40, height = 30, units = "cm")



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
  Analyte = case_when(
    Analyte == "Copper" ~ "Copper (μg/L)",
    Analyte == "Phosphorus" ~ "Phosphorus (mg/L)",
    Analyte == "TSS" ~ "TSS (mg/L)"
  )
)

parent_dataset_quartiles <- merge(parent_dataset_quartiles, rain_zone_mapping, by = "Location")

parent_dataset_quartiles$Quartile <- factor(parent_dataset_quartiles$Quartile, levels = c("10th Percentile",  "25th Percentile","50th Percentile", "75th Percentile", "90th Percentile"))

parent_dataset_quartiles %>% ggplot() + 
  geom_boxplot(aes(x = Flow, y = EMC), outlier.shape = NA) +
  geom_jitter(aes(x = Flow, y = EMC, color = smallerCountBins, shape = smallerCountBins), size = 2, position=position_jitter(width=0.1, height=0.1), alpha = 0.9, stroke=2) +
  #scale_shape_manual(values = c(0, 1, 2, 5, 6)) +
  facet_grid(Analyte~Quartile, scales = "free", switch = "y") + 
  scale_y_continuous(trans = "log10") + 
  stat_compare_means(aes(x = Flow, y = EMC), label = "p.signif", size = 6, vjust = 0.7, label.x.npc = "center") + 
  theme_bw(base_size = 24) +
  theme(legend.position = "top") + 
  labs(color  = "Number of Events in Parent Dataset", shape = "Number of Events in Parent Dataset") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), strip.placement = "outside", strip.background = element_blank())
ggsave(here::here("Plots", "parent_dataset_by_analyte_quartile.jpg"), width = 40, height = 50, units = "cm")

parent_dataset_quartiles %>% ggplot() + 
  geom_boxplot(aes(x = Flow, y = EMC), outlier.shape = NA) +
  geom_jitter(aes(x = Flow, y = EMC, color = smallerCountBins, shape = smallerCountBins), size = 2, position=position_jitter(width=0.1, height=0.1), alpha = 0.9, stroke = 2) +
  #scale_shape_manual(values = c(0, 1, 2, 5, 6)) +
  facet_grid(Analyte~., scales = "free", switch = "y") + 
  scale_y_continuous(trans = "log10") + 
  stat_compare_means(aes(x = Flow, y = EMC), label = "p.signif", size = 6, vjust = 0.7, label.x.npc = "center") + 
  theme_bw(base_size = 24) +
  theme(legend.position = "top") + 
  labs(color  = "Number of Events in Parent Dataset", shape = "Number of Events in Parent Dataset") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), strip.placement = "outside", strip.background = element_blank())
ggsave(here::here("Plots", "parent_dataset_by_analyte.jpg"), width = 40, height = 50, units = "cm")


parent_dataset_quartiles %>% ggplot() + 
  geom_boxplot(aes(x = as.character(Rain_Zone), y = EMC), outlier.shape = NA) +
  geom_jitter(aes(x = as.character(Rain_Zone), y = EMC, color = smallerCountBins, shape = smallerCountBins), size = 2, position=position_jitter(width=0.1, height=0.1), alpha = 0.9, stroke = 2) +
  #scale_shape_manual(values = c(0, 1, 2, 5, 6)) +
  facet_grid(Analyte~Quartile, scales = "free", switch = "y") + 
  scale_y_continuous(trans = "sqrt") + 
  stat_compare_means(aes(x = as.character(Rain_Zone), y = EMC), label = "p.format", size = 6, vjust = 0.6) + 
  theme_bw(base_size = 24) +
  theme(legend.position = "top") + xlab("Rain Zone") +
  labs(color  = "Number of Events in Parent Dataset", shape = "Number of Events in Parent Dataset") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), strip.placement = "outside", strip.background = element_blank())
ggsave(here::here("Plots", "parent_dataset_by_rainzone_analyte_quartile.jpg"), width = 40, height = 50, units = "cm")

parent_dataset_quartiles %>% ggplot() + 
  geom_boxplot(aes(x = as.character(Rain_Zone), y = EMC), outlier.shape = NA) +
  geom_jitter(aes(x = as.character(Rain_Zone), y = EMC, color = smallerCountBins, shape = smallerCountBins), size = 2, position=position_jitter(width=0.1, height=0.1), alpha = 0.9, stroke = 2) +
  #scale_shape_manual(values = c(0, 1, 2, 5, 6)) +
  facet_grid(Analyte~., scales = "free", switch = "y") + 
  scale_y_continuous(trans = "sqrt") + 
  stat_compare_means(aes(x = as.character(Rain_Zone), y = EMC), label = "p.format", size = 6, vjust = 0.6) + 
  theme_bw(base_size = 24) +
  theme(legend.position = "top") + xlab("Rain Zone") +
  labs(color  = "Number of Events in Parent Dataset", shape = "Number of Events in Parent Dataset") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), strip.placement = "outside", strip.background = element_blank())
ggsave(here::here("Plots", "parent_dataset_by_rainzone_analyte.jpg"), width = 40, height = 50, units = "cm")


parent_dataset_quartiles %>% filter(Flow == "Inflow") %>% ggplot() + 
  geom_boxplot(aes(x = as.character(Rain_Zone), y = EMC), outlier.shape = NA) +
  geom_jitter(aes(x = as.character(Rain_Zone), y = EMC, color = smallerCountBins, shape = smallerCountBins), size = 2, position=position_jitter(width=0.1, height=0.1), alpha = 0.9, stroke = 2) +
  #scale_shape_manual(values = c(0, 1, 2, 5, 6)) +
  facet_grid(Analyte~Quartile, scales = "free", switch = "y") + 
  scale_y_continuous(trans = "sqrt") + 
  stat_compare_means(aes(x = as.character(Rain_Zone), y = EMC), label = "p.format", size = 6, vjust = 0.6) + 
  theme_bw(base_size = 24) +
  theme(legend.position = "top") + xlab("Rain Zone") +
  labs(color  = "Number of Events in Parent Dataset", shape = "Number of Events in Parent Dataset") + ggtitle("Inflow") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), strip.placement = "outside", strip.background = element_blank())
ggsave(here::here("Plots", "parent_dataset_by_rainzone_analyte_quartile_inflow.jpg"), width = 40, height = 50, units = "cm")

parent_dataset_quartiles %>% filter(Flow == "Inflow") %>% ggplot() + 
  geom_boxplot(aes(x = as.character(Rain_Zone), y = EMC), outlier.shape = NA) +
  geom_jitter(aes(x = as.character(Rain_Zone), y = EMC, color = smallerCountBins, shape = smallerCountBins, size = Quartile), position=position_jitter(width=0.1, height=0.1), alpha = 0.9, stroke = 2) +
  #scale_shape_manual(values = c(0, 1, 2, 5, 6)) + scale_size_manual(values = c(1, 2.5, 5, 7.5, 10, 12.5, 15)) +
  facet_grid(Analyte~., scales = "free", switch = "y") + 
  scale_y_continuous(trans = "sqrt") + 
  stat_compare_means(aes(x = as.character(Rain_Zone), y = EMC), label = "p.format", size = 6, vjust = 0.6) + 
  theme_bw(base_size = 24) +
  theme(legend.position = "top") + xlab("Rain Zone") +
  labs(color  = "Number of Events in Parent Dataset", shape = "Number of Events in Parent Dataset") + ggtitle("Inflow") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), strip.placement = "outside", strip.background = element_blank(), legend.box = "vertical")
ggsave(here::here("Plots", "parent_dataset_by_rainzone_analyte_inflow.jpg"), width = 40, height = 50, units = "cm")

parent_dataset_quartiles %>% filter(Flow == "Outflow") %>% ggplot() + 
  geom_boxplot(aes(x = as.character(Rain_Zone), y = EMC), outlier.shape = NA) +
  geom_jitter(aes(x = as.character(Rain_Zone), y = EMC, color = smallerCountBins, shape = smallerCountBins), size = 2, position=position_jitter(width=0.1, height=0.1), alpha = 0.9, stroke = 2) +
  #scale_shape_manual(values = c(0, 1, 2, 5, 6)) +
  facet_grid(Analyte~Quartile, scales = "free", switch = "y") + 
  scale_y_continuous(trans = "sqrt") + 
  stat_compare_means(aes(x = as.character(Rain_Zone), y = EMC), label = "p.format", size = 6, vjust = 0.6) + 
  theme_bw(base_size = 24) +
  theme(legend.position = "top") + xlab("Rain Zone") +
  labs(color  = "Number of Events in Parent Dataset", shape = "Number of Events in Parent Dataset") + ggtitle("Outflow") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), strip.placement = "outside", strip.background = element_blank())
ggsave(here::here("Plots", "parent_dataset_by_rainzone_analyte_quartile_outflow.jpg"), width = 40, height = 50, units = "cm")

parent_dataset_quartiles %>% filter(Flow == "Outflow") %>% ggplot() + 
  geom_boxplot(aes(x = as.character(Rain_Zone), y = EMC), outlier.shape = NA) +
  geom_jitter(aes(x = as.character(Rain_Zone), y = EMC, color = smallerCountBins, shape = smallerCountBins, size = Quartile), position=position_jitter(width=0.1, height=0.1), alpha = 0.9, stroke = 2) +
  #scale_shape_manual(values = c(0, 1, 2, 5, 6)) +  scale_size_manual(values = c(1, 2.5, 5, 7.5, 10, 12.5, 15)) +
  facet_grid(Analyte~., scales = "free", switch = "y") + 
  scale_y_continuous(trans = "sqrt") + 
  stat_compare_means(aes(x = as.character(Rain_Zone), y = EMC), label = "p.format", size = 6, vjust = 0.6) + 
  theme_bw(base_size = 24) +
  theme(legend.position = "top") + xlab("Rain Zone") +
  labs(color  = "Number of Events in Parent Dataset", shape = "Number of Events in Parent Dataset") + ggtitle("Outflow") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), strip.placement = "outside", strip.background = element_blank(), legend.box = "vertical")
ggsave(here::here("Plots", "parent_dataset_by_rainzone_analyte_outflow.jpg"), width = 40, height = 50, units = "cm")


for_table <- data %>% select(Location, Flow, Rain_Zone, Analyte, smallerCountBins) %>% unique() %>% mutate(
  Rain_Zone = case_when(
    Rain_Zone == 1 ~ "1 (Great Lakes)",
    Rain_Zone == 2 ~ "2 (North East)",
    Rain_Zone == 3 ~ "3 (South East)",
    Rain_Zone == 4 ~ "4 (Lower Mississippi Valley)",
    Rain_Zone == 5 ~ "5 (Texas)",
    Rain_Zone == 6 ~ "6 (South West)",
    Rain_Zone == 7 ~ "7 (Northwest)",
    Rain_Zone == 8 ~ "8 (California)",
    Rain_Zone == 9 ~ "9 (Rocky Mountains)"
  ),
  Analyte = case_when(
    Analyte == "Copper" ~ "Dissolved Copper",
    Analyte == "Phosphorus" ~ "Total Phosphorus",
    TRUE ~ Analyte
  )
) %>% select(-Location) %>% group_by(Rain_Zone, Flow, Analyte, smallerCountBins) %>% dplyr::summarize(
  Count = n()
) %>% mutate(
  smallerCountBins = factor(smallerCountBins, c("20-30", "31-60", "61-90", ">90")),
  Flow = factor(Flow, c("Inflow", "Outflow"))
) %>% pivot_wider(names_from = c(Analyte, smallerCountBins), values_from = Count)

for_table <- for_table[, c("Rain_Zone", "Flow", "Dissolved Copper_20-30", "Total Phosphorus_20-30", "TSS_20-30", "Dissolved Copper_31-60", "Total Phosphorus_31-60", "TSS_31-60", "Dissolved Copper_61-90", "Total Phosphorus_61-90", "TSS_61-90", "Dissolved Copper_>90", "Total Phosphorus_>90", "TSS_>90")]

## Final Figures ##
# 6/3/26: changed "data" to "data_BMPCat_s" (subset after removing irrelevant BMP Categories)
num_samples_analyte_summary <- data_summary(data_BMPCat_s, varname = "rpd", groupnames = c("N_numerical", "Quartile", "smallerCountBins", "Analyte"))
num_samples_analyte_summary$N_numerical <- factor(num_samples_analyte_summary$N_numerical, levels = c(5, 10, 15, 20, 25))
num_samples_analyte_summary$Quartile <- factor(num_samples_analyte_summary$Quartile, levels = c(0.1, 0.25, 0.5, 0.75, 0.9))

# sanity check: mean RPD b/w subsets with 25 monitored events and parent datasets with 20-30 events should be close to 0
s25_p20_30 <- num_samples_analyte_summary %>%
  filter(N_numerical == 25) %>%
  filter(smallerCountBins == "20-30")

  # create boxplots for each analyte to visualize RPD
  bp_s25 <- ggplot(s25_p20_30, aes(x = Analyte, y = rpd, fill = Analyte)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.5) +  # Hide default outliers
    geom_jitter(width = 0.2, size = 2, alpha = 0.7) + # Add jittered points
    labs(x = NULL, y = "RPD (%)", title = "RPD for subsets w/ 25 events &\nparent datasets w/ 20-30 events") +
    theme_bw(base_size = 20)
  
  bp_s25
  
  ggsave(here::here("Plots", "Final", "rpd_boxplots_subsets_n25_parent_n20-30_byAnalyte_20260612.jpg"), width = 20, height = 16, units = "cm")
  

# plot with all analytes
  # create separate facets for Phosphorus:Inflow and Phosphorus:Outflow
  num_samples_P_flow_summary <- data_summary(data_BMPCat_s, varname = "rpd", groupnames = c("N_numerical", "Quartile", "smallerCountBins", "Analyte", "Flow")) %>%
    filter(Analyte == "Phosphorus") %>%
    # create facet group == Analyte + Flow
    mutate(facet_group = paste0(Analyte, ": ", Flow)) %>%
    # drop Flow column for joining with other data
    select(-Flow) %>%
    # format grouping vars
    mutate(N_numerical = factor(N_numerical, levels = c(5, 10, 15, 20, 25))) %>%
    mutate(Quartile = factor(Quartile, levels = c(0.1, 0.25, 0.5, 0.75, 0.9)))
  
  num_samples_Cu_TSS_summary <- num_samples_analyte_summary %>%
    filter(Analyte != "Phosphorus") %>%
    # create facet group == Analyte
    mutate(facet_group = Analyte)
    
  num_samples_analyte_summary <- rbind(num_samples_P_flow_summary, num_samples_Cu_TSS_summary) %>%
    # format facet group
    mutate(facet_group = factor(facet_group, levels = c("Copper", "Phosphorus: Inflow",
                                                        "Phosphorus: Outflow", "TSS")))
  
set.seed(36)  # ensures reproducible jitter positions for geom_point

num_samples_analyte_summary %>%
  mutate(Quartile = case_when(
    Quartile == 0.1 ~ "10th Percentile",
    Quartile == 0.25 ~ "25th Percentile",
    Quartile == 0.5 ~ "50th Percentile",
    Quartile == 0.75 ~ "75th Percentile",
    Quartile == 0.9 ~ "90th Percentile"
  )) %>%
  ggplot(aes(x = N_numerical, y = rpd, color = smallerCountBins, group = smallerCountBins)) + 
  geom_line(stat = "smooth", method = "lm", se = FALSE, formula = y ~ I(log(x)), alpha = 0.7, size = 1.5) + 
  geom_errorbar(aes(ymin = rpd - sd, ymax = rpd + sd), width = 0.2, size = 1.5, 
                position = position_jitter(width = 0.1, height = 0, seed = 38)) +
  geom_point(aes(shape = smallerCountBins), size = 8, fill = NA, stroke = 2, 
             position = position_jitter(width = 0.1, height = 0, seed = 38)) + 
  scale_shape_manual(name = "Number of Events in Parent Dataset",
                     labels = c("20-30", "31-60", "61-90", ">90"),
                     values = c(16, 17, 15, 4)) +
  facet_grid(Quartile~facet_group) +
  #geom_hline(yintercept=20, linetype="dashed", color = "black", size=1) + 
  #geom_hline(yintercept=10, linetype="dashed", color = "black", size=1) + 
  scale_y_continuous(breaks = seq(0, 70, 10), limits = c(-3,71)) +
  #ylim(c(-5, 100)) +
  xlab("Number of Monitored Events in Data Subset") + ylab("Mean Relative Percent Difference (%)") +
  theme_bw(base_size = 40) +  theme(strip.placement = "outside", panel.spacing = unit(1, "lines"), strip.background = element_blank()) + theme(legend.position = "top", legend.box = "vertical") +
  #labs(color  = "Number of Events in Parent Dataset", shape = "Number of Events in Parent Dataset") +
  guides(fill = "none") +
  scale_color_manual(name = "Number of Events in Parent Dataset",
                     labels = c("20-30", "31-60", "61-90", ">90"),
                     values = c("#C7E9B4", "#41B6C4", "#225EA8","#081D58")) +
  theme(legend.position = "top") +
  theme(panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) # panel.grid.major = element_blank(), 
ggsave(here::here("Plots", "Final", "num_samples_analyte_PFlow_summary_ymax71_jitter_Brewer_20260617.jpg"), width = 66, height = 70, units = "cm")

# plot with TSS - inflow and outflow separately
  # group by flow type
num_samples_analyte_flow_summary <- data_summary(data_BMPCat_s, varname = "rpd", groupnames = c("N_numerical", "Quartile", "smallerCountBins", "Analyte", "Flow"))
num_samples_analyte_flow_summary$N_numerical <- factor(num_samples_analyte_flow_summary$N_numerical, levels = c(5, 10, 15, 20, 25))
num_samples_analyte_flow_summary$Quartile <- factor(num_samples_analyte_flow_summary$Quartile, levels = c(0.1, 0.25, 0.5, 0.75, 0.9))

set.seed(36)  # ensures reproducible jitter positions for geom_point

num_samples_analyte_flow_summary %>%
  mutate(Quartile = case_when(
    Quartile == 0.1 ~ "10th Percentile",
    Quartile == 0.25 ~ "25th Percentile",
    Quartile == 0.5 ~ "50th Percentile",
    Quartile == 0.75 ~ "75th Percentile",
    Quartile == 0.9 ~ "90th Percentile"
  )) %>%
  # filter for TSS only
  filter(Analyte == "TSS") %>%
  # change Flow categories to include analyte type
  mutate(Flow = case_when(
    Flow == "Inflow" ~ "TSS: Inflow",
    Flow == "Outflow" ~ "TSS: Outflow",
    .default = Flow
  )) %>%
  ggplot(aes(x = N_numerical, y = rpd, color = smallerCountBins, group = smallerCountBins)) + 
  geom_line(stat = "smooth", method = "lm", se = FALSE, formula = y ~ I(log(x)), alpha = 0.7, size = 1.5) + 
  geom_errorbar(aes(ymin = rpd - sd, ymax = rpd + sd), width = 0.2, size = 1.5, 
                position = position_jitter(width = 0.1, height = 0, seed = 38)) +
  geom_point(aes(shape = smallerCountBins), size = 8, fill = NA, stroke = 2, 
             position = position_jitter(width = 0.1, height = 0, seed = 38)) + 
  scale_shape_manual(name = "Number of Events in Parent Dataset",
                     labels = c("20-30", "31-60", "61-90", ">90"),
                     values = c(16, 17, 15, 4)) +
  facet_grid(Quartile~Flow) +
  #geom_hline(yintercept=20, linetype="dashed", color = "black", size=1) + 
  #geom_hline(yintercept=10, linetype="dashed", color = "black", size=1) + 
  scale_y_continuous(breaks = seq(0, 70, 10), limits = c(-3,71)) +
  #ylim(c(-5, 100)) +
  xlab("Number of Monitored Events in Data Subset") + ylab("Mean Relative Percent Difference (%)") +
  #ggtitle("Phosphorus") +
  theme_bw(base_size = 40) +  theme(strip.placement = "outside", panel.spacing = unit(1, "lines"), strip.background = element_blank()) + theme(legend.position = "top", legend.box = "vertical") +
  #labs(color  = "Number of Events in Parent Dataset", shape = "Number of Events in Parent Dataset") +
  guides(fill = "none") +
  scale_color_manual(name = "Number of Events in Parent Dataset",
                     labels = c("20-30", "31-60", "61-90", ">90"),
                     values = c("#C7E9B4", "#41B6C4", "#225EA8","#081D58")) +
  theme(legend.position = "top") +
  theme(panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) # panel.grid.major = element_blank(), 
ggsave(here::here("Plots", "Final", "num_samples_TSS_flow_summary_jitter_Brewer_20260618.jpg"), width = 50, height = 60, units = "cm")

#display.brewer.all(colorblindFriendly = TRUE)
#brewer.pal(9, "Blues")
#"#F7FBFF" "#DEEBF7" "#C6DBEF" "#9ECAE1" "#6BAED6" "#4292C6" "#2171B5" "#08519C" "#08306B"
#brewer.pal(9, "YlGnBu")
# "#FFFFD9" "#EDF8B1" "#C7E9B4" "#7FCDBB" "#41B6C4" "#1D91C0" "#225EA8" "#253494" "#081D58"

num_samples_analyte_summary %>%
  mutate(Quartile = case_when(
    Quartile == 0.1 ~ "10th Percentile",
    Quartile == 0.25 ~ "25th Percentile",
    Quartile == 0.5 ~ "50th Percentile (Median)",
    Quartile == 0.75 ~ "75th Percentile",
    Quartile == 0.9 ~ "90th Percentile"
  )) %>% filter(Analyte == "TSS", Quartile == "50th Percentile (Median)") %>%
  ggplot(aes(x = N_numerical, y = rpd, color = smallerCountBins, group = smallerCountBins)) + 
  geom_point(aes(shape = smallerCountBins), size = 3, fill = NA) + 
  #scale_shape_manual(values = c(0, 1, 2, 5, 6)) +
  stat_smooth(method = "lm", se = FALSE, formula = y ~ I(log(x))) + 
  geom_errorbar(aes(ymin = rpd - sd, ymax = rpd + sd), width = .2) +
  facet_grid(Analyte~Quartile) +
  ylim(c(-5, 100)) +
  geom_hline(yintercept=20, linetype="dashed", color = "black", size=1) + 
  geom_hline(yintercept=10, linetype="dashed", color = "black", size=1) + 
  geom_rect(aes(xmin = 1.5, xmax = 2.5,  ymin = 5, ymax = 48),
            fill = "transparent", color = "red", size = 0.6) + 
  geom_text(
    aes(x = 2, y = 55, label = "b"),
    size = 4, vjust = 0, hjust = 0.5, color = "black", check_overlap = TRUE
  ) +
  geom_rect(aes(xmin = 4.5, xmax = 5.5,  ymin = -2, ymax = 32),
            fill = "transparent", color = "red", size = 0.6) + 
  geom_text(
    aes(x = 5, y = 39, label = "a"),
    size = 4, vjust = 0, hjust = 0.5, color = "black", check_overlap = TRUE
  ) +
  xlab("Number of Monitored Events in Data Subset") + ylab("Mean Relative Percent Difference (%)") +
  theme_bw(base_size = 12) +  theme(strip.placement = "outside", panel.spacing = unit(1, "lines"), strip.background = element_blank()) + 
  theme(legend.position = "right", legend.box = "vertical") +
  labs(color  = "Number of Events\nin Parent Dataset", shape = "Number of Events\nin Parent Dataset") +
  guides(fill = "none") + guides(shape = guide_legend(ncol = 2), color = guide_legend(ncol = 2)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(here::here("Plots", "Final", "num_samples_analyte_summary_tss_50th.jpg"), width = 18, height = 10, units = "cm")

# 6/3/26: changed "data" to "data_BMPCat_s" (subset after removing irrelevant BMP Categories)
data_BMPCat_s <- data_BMPCat_s %>% mutate(
  percentBelowDetect = (NumObservationsBelowDetect/NumObservations) * 100
)

data_BMPCat_s$belowDetectPercentBins <- cut(data_BMPCat_s$percentBelowDetect, breaks = c(-Inf, 0, 9.99, 19.99, 39.99, Inf), labels = c("0%", "0.001-10%", "10-20%", "20-40%", ">40%"))

num_samples_below_detect_summary <- data_summary(data_BMPCat_s, varname = "rpd", groupnames = c("N_numerical", "Quartile", "belowDetectPercentBins", "Analyte"))
num_samples_below_detect_summary$N_numerical <- factor(num_samples_below_detect_summary$N_numerical, levels = c(5, 10, 15, 20, 25))
num_samples_below_detect_summary$Quartile <- factor(num_samples_below_detect_summary$Quartile, levels = c(0.1, 0.25, 0.5, 0.75, 0.9))
num_samples_below_detect_summary$belowDetectPercentBins <- factor(num_samples_below_detect_summary$belowDetectPercentBins, levels = c("0%", "0.001-10%", "10-20%", "20-40%", ">40%"))

# make table: count of parent datasets by categories for "% of WQ results below detection limit"
  n_parent_by_PctBlwDL <- num_samples_below_detect_summary %>%
    group_by(Analyte, Quartile, belowDetectPercentBins) %>%
    summarise(n_parent = n())
  
  # save as CSV
  #write.csv(n_parent_by_PctBlwDL, here::here("Data", "Number of Parent Datasets by Below Detection Limit Plot Category_20260619.csv"), row.names = FALSE)
  
  # # check data for issues
  # nbd_Phos_10pctl <- num_samples_below_detect_summary %>%
  #   filter(Analyte=="Phosphorus" & Quartile == 0.1)

set.seed(36)  # ensures reproducible jitter positions for geom_point

num_samples_below_detect_summary %>%
  mutate(Quartile = case_when(
    Quartile == 0.1 ~ "10th Percentile",
    Quartile == 0.25 ~ "25th Percentile",
    Quartile == 0.5 ~ "50th Percentile",
    Quartile == 0.75 ~ "75th Percentile",
    Quartile == 0.9 ~ "90th Percentile"
  )) %>%
  # filter only 10th and 25th percentile (other percentiles had low RPDs)
  filter(Quartile=="10th Percentile" | Quartile=="25th Percentile") %>%
  ggplot(aes(x = N_numerical, y = rpd, color = belowDetectPercentBins, group = belowDetectPercentBins)) + 
  geom_line(stat = "smooth", method = "lm", se = FALSE, formula = y ~ I(log(x)), size = 1.5, alpha = 0.7) + 
  geom_errorbar(aes(ymin = rpd - sd, ymax = rpd + sd), width = .2, size = 1.5,
                position = position_jitter(width = 0.1, height = 0, seed = 38)) +
  facet_grid(Quartile~Analyte) +
  #geom_hline(yintercept=20, linetype="dashed", color = "black", size=1) + 
  #geom_hline(yintercept=10, linetype="dashed", color = "black", size=1) + 
  geom_point(aes(shape = belowDetectPercentBins), size = 8, stroke = 2, fill = NA,
             position = position_jitter(width = 0.1, height = 0, seed = 38)) + 
  scale_shape_manual(name = "Percentage of Water Quality \nResults Below Detection Limit",
                     labels = c("0%", "0.001-10%", "10-20%", "20-40%", ">40%"),
                     values = c(16, 17, 15, 4, 7)) +
  xlab("Number of Monitored Events in Data Subset") + ylab("Mean Relative Percent Difference (%)") +
  theme_bw(base_size = 40) +  theme(strip.placement = "outside", panel.spacing = unit(1, "lines"), strip.background = element_blank()) + theme(legend.position = "top", legend.box = "vertical") +
  #labs(color  = "Percentage of Water Quality \nResults Below Detection Limit", shape = "Percentage of Water Quality \nResults Below Detection Limit") +
  guides(fill = "none") +
  scale_color_manual(name = "Percentage of Water Quality \nResults Below Detection Limit",
                     labels = c("0%", "0.001-10%", "10-20%", "20-40%", ">40%"),
                     values = c("#C7E9B4", "#41B6C4", "#2C7FB8","#253494", "#081D58")) +
  scale_y_continuous(breaks = seq(0, 120, 20)) + #, limits = c(-15,116)
  #ylim(c(-10, 100)) +
  guides(legend = guide_legend(ncol = 2)) +
  theme(legend.position = "top") +
  theme(panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) # panel.grid.major = element_blank(), 
ggsave(here::here("Plots", "Final", "num_samples_percent_below_detect_summary_jitter_Brewer_low_percentiles_20260617.jpg"), height = 40, width = 50, units = "cm")

# num_samples_below_detect_summary %>%
#   mutate(Quartile = case_when(
#     Quartile == 0.1 ~ "10th Percentile",
#     Quartile == 0.25 ~ "25th Percentile",
#     Quartile == 0.5 ~ "50th Percentile",
#     Quartile == 0.75 ~ "75th Percentile",
#     Quartile == 0.9 ~ "90th Percentile"
#   )) %>%
#   ggplot(aes(x = N_numerical, y = rpd, color = belowDetectPercentBins, group = belowDetectPercentBins)) + 
#   geom_line(stat = "smooth", method = "lm", se = FALSE, formula = y ~ I(log(x)), size = 1.5, alpha = 0.7) + 
#   geom_errorbar(aes(ymin = rpd - sd, ymax = rpd + sd), width = .2, size = 1.5) +
#   facet_grid(Quartile~Analyte) +
#   geom_hline(yintercept=20, linetype="dashed", color = "black", size=1) + 
#   geom_hline(yintercept=10, linetype="dashed", color = "black", size=1) + 
#   geom_point(aes(shape = belowDetectPercentBins), size = 8, stroke = 2, fill = NA) + 
#   #scale_shape_manual(values = c(16, 17, 15, 4, 7)) +
#   xlab("Number of Monitored Events in Data Subset") + ylab("Mean Relative Percent Difference (%)") +
#   theme_bw(base_size = 40) +  theme(strip.placement = "outside", panel.spacing = unit(1, "lines"), strip.background = element_blank()) + theme(legend.position = "top", legend.box = "vertical") +
#   labs(color  = "Proportion of Water Quality \nResults Below Detection Limit", shape = "Proportion of Water Quality \nResults Below Detection Limit") +
#   guides(fill = "none") +
#   ylim(c(-10, 100)) +
#   guides(legend = guide_legend(ncol = 2)) +
#   theme(legend.position = "top") +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) 
# ggsave(here::here("Plots", "Final", "num_samples_percent_below_detect_summary.jpg"), height = 60, width = 50, units = "cm")

coefficient_of_variation <- data %>% select(Location, Flow, Rain_Zone, Analyte, coefficientOfVariation, NumObservations) %>% unique() %>% mutate(
  Rain_Zone = case_when(
    Rain_Zone == 1 ~ "1 (Great Lakes)",
    Rain_Zone == 2 ~ "2 (North East)",
    Rain_Zone == 3 ~ "3 (South East)",
    Rain_Zone == 4 ~ "4 (Lower Mississippi Valley)",
    Rain_Zone == 5 ~ "5 (Texas)",
    Rain_Zone == 6 ~ "6 (South West)",
    Rain_Zone == 7 ~ "7 (Northwest)",
    Rain_Zone == 8 ~ "8 (California)",
    Rain_Zone == 9 ~ "9 (Rocky Mountains)"
  )
) %>% select(-Location)

coefficient_of_variation %>% ggplot() + 
  geom_point(aes(x = NumObservations, y = coefficientOfVariation, shape = Flow, color = Flow), stroke = 1.75, size = 9) + 
  #scale_shape_manual(values = c(9, 13)) +
  facet_grid(Analyte~str_wrap(Rain_Zone, 15)) +
  xlab("Number of Monitored Events in Parent Dataset") + ylab("Coefficient of Variation") +
  theme_bw(base_size = 30) +  theme(strip.placement = "outside", panel.spacing = unit(1, "lines"), strip.background = element_blank()) + theme(legend.position = "top", legend.box = "vertical") +
  labs(color  = "Flow Type", shape = "Flow Type") +
  geom_vline(linetype = "dashed", alpha = 0.4, xintercept = c(50, 100, 150, 200)) + 
  guides(fill = "none") + coord_trans(x = "log10") + scale_x_continuous(breaks = c(0, 20, 50, 100, 200, 600)) +
  theme(legend.position = "top") + ylim(c(0, 3)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(here::here("Plots", "Final", "coefficient_of_variation_summary.jpg"), width = 64, height = 40, units = "cm")

# stats to see if CV is significantly different for subsets by Rain Zone and
  # subsets by Inflow vs. Outflow

  ### Rain Zone
  # plot histograms

  # calculate binwidth using Freedman–Diaconis Rule (good for skewed data)
  fd_binwidth <- function(x) {
    2 * IQR(x, na.rm = TRUE) / length(x)^(1/3)
  }

  hist_by_rz <- ggplot(coefficient_of_variation, aes(x = coefficientOfVariation, fill = Rain_Zone)) + 
    geom_histogram(color = "black", alpha = 0.5, position = "identity", binwidth = fd_binwidth) +
    facet_wrap(~ Rain_Zone, scales = "free_y") +
    labs(x = "Coefficient of Variation",
         y = "Count",
         fill = "Rain Zone") +
    theme_bw(base_size = 20)
  
  hist_by_rz
  
  # save as png
  ggsave(here::here("Plots", "Rain_zone_histograms_FD_binwidth.png"), width = 1600, height = 1100, units = "px", dpi = 110)
  
  # check for normality
  rzs <- unique(coefficient_of_variation$Rain_Zone)
  n_rz <- length(rzs)
  
  st_result_rz <- tibble(
    Rain_Zone = rzs,
    p_value = numeric(n_rz)
  )
  
  for (i in 1:n_rz) {
    rz <- rzs[i]
    
    # subset by Rain Zone
    rz_s <- coefficient_of_variation %>%
      filter(Rain_Zone == rz)
    
    # run Shapiro-Wilk test for normality
    st_rz <- shapiro.test(rz_s$coefficientOfVariation)
    
    # save p-value in output table
    st_result_rz$p_value[i] <- st_rz$p.value
  }
  # all p-values < 0.05 --> data are not normally distributed
  
  # Kruskal-Wallis test for significant differences b/w groups (non-parametric)
  kw_rz <- kruskal.test(coefficientOfVariation ~ Rain_Zone, data = coefficient_of_variation)
  kw_rz # p = 3.6e-05 << 0.05 --> significant difference b/w groups
  
  # Posthoc: Dunn's Test w/ Holm correction for p-values
  library(FSA)
  
  dt_rz <- dunnTest(coefficientOfVariation ~ Rain_Zone, 
                    data = coefficient_of_variation,
                    method = "holm")
  
  dunn_result <- dt_rz$res
  
  # save as CSV
  #write.csv(dunn_result, "Dunn's Test Results_Differences in COV by Rain Zone.csv", row.names = FALSE)

  ### Inflow vs. Outflow
  # plot histograms
  hist_by_flow <- ggplot(coefficient_of_variation, aes(x = coefficientOfVariation, fill = Flow)) + 
    geom_histogram(color = "black", alpha = 0.5, position = "identity") +
    labs(x = "Coefficient of Variation",
         y = "Count",
         fill = "Flow Type") +
    theme_bw(base_size = 20)
  
  hist_by_flow
  
  # save as png
  ggsave(here::here("Plots", "Inflow_Outflow_histograms.png"), width = 900, height = 700, units = "px", dpi = 110)

  # run Shapiro-Wilk test for normality
  CV_inf_s <- coefficient_of_variation %>%
    filter(Flow == "Inflow")
  
  st_inf <- shapiro.test(CV_inf_s$coefficientOfVariation)
  # p-value < 2.2e-16 --> data are not normally distributed
  
  CV_out_s <- coefficient_of_variation %>%
    filter(Flow == "Outflow")
  
  st_out <- shapiro.test(CV_out_s$coefficientOfVariation)
  # p-value < 2.2e-16 --> data are not normally distributed
  
  # since non-normal, test difference with Mann-Whitney Test
    # Mann-Whitney interpretation
      # Null hypothesis (H₀): The two groups come from the same distribution (no shift in medians).
      # Alternative hypothesis (H₁): The distributions differ (median shift ≠ 0 by default).
  
  MW_flow <- wilcox.test(CV_inf_s$coefficientOfVariation, CV_out_s$coefficientOfVariation,
                         paired = FALSE, exact = TRUE)
  
  MW_flow # p = 0.08 > 0.05 --> Medians and distributions do not significantly differ
  

ratios <- read.csv(here::here("Data", "parent_dataset_quartiles.csv")) %>% filter(Quartile %in% c(0.1, 0.25, 0.75, 0.9)) %>%
  pivot_wider(names_from = Quartile, values_from = EMC) %>% mutate(
    `P[90]/P[10]` = `0.9`/`0.1`,
    `P[75]/P[25]` = `0.75`/`0.25`
  ) %>% select(-`0.1`, -`0.25`, -`0.75`, -`0.9`)

ratios_w_rz <- merge(ratios, rain_zone_mapping, by = "Location")

ratios_w_rz <- ratios_w_rz %>% mutate(
  Rain_Zone = case_when(
    Rain_Zone == 1 ~ "1 (Great Lakes)",
    Rain_Zone == 2 ~ "2 (North East)",
    Rain_Zone == 3 ~ "3 (South East)",
    Rain_Zone == 4 ~ "4 (Lower Mississippi Valley)",
    Rain_Zone == 5 ~ "5 (Texas)",
    Rain_Zone == 6 ~ "6 (South West)",
    Rain_Zone == 7 ~ "7 (Northwest)",
    Rain_Zone == 8 ~ "8 (California)",
    Rain_Zone == 9 ~ "9 (Rocky Mountains)"
  )
) %>% select(-Location) %>% pivot_longer(cols = c(`P[90]/P[10]`, `P[75]/P[25]`), names_to = "Ratio", values_to = "P[75]/P[25] and P[90]/P[10] Ratio")

ratios_w_rz %>% ggplot() +
  geom_point(aes(x = NumSamples, y = `P[75]/P[25] and P[90]/P[10] Ratio`, color = Ratio, shape = Flow), stroke = 1.75, size = 9) + 
  geom_line(aes(x = NumSamples, y = `P[75]/P[25] and P[90]/P[10] Ratio`, group = NumSamples)) + 
  #scale_shape_manual(values = c(9, 13)) +
  scale_color_manual(values = c("#F8766D", "#00BFC4"), labels = c(expression(P[75]/P[25]), expression(P[90]/P[10]))) +
  facet_grid(str_wrap(Rain_Zone, 15)~Analyte, scales = "free_y") +
  xlab("Number of Monitored Events in Parent Dataset") + ylab(expression(P[75]/P[25]~~and~~P[90]/P[10]~~Ratio)) +
  theme_bw(base_size = 30) +  theme(strip.placement = "outside", panel.spacing = unit(1, "lines"), strip.background = element_blank()) + theme(legend.position = "top", legend.box = "vertical") +
  labs(color  = "Ratio", shape = "Flow") +
  geom_vline(linetype = "dashed", alpha = 0.4, xintercept = c(50, 100, 150, 200)) + 
  guides(fill = "none") +
  theme(legend.position = "top") + coord_trans(x = "log10", y = "log10") + scale_x_continuous(breaks = c(0, 20, 50, 100, 200, 600)) +
  scale_y_continuous(breaks = c(0, 10, 50, 100, 200)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(here::here("Plots", "Final", "ratios_summary.jpg"), width = 50, height = 65, units = "cm")

  ### BMP Category
# 4/29/26: make coefficient of variation dataframe with BMP_Category
cv_BMPType <- data %>%
  select(Rain_Zone, Location, Analyte, Flow, coefficientOfVariation, NumObservations) %>%
  distinct() %>%
  left_join(metadata_rz_fix, 
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
    .default = BMP_Category
  )) 

cv_BMPType_plt <- cv_BMPType %>%
  # filter out NSWQ record (n=1)
  filter(BMP_Category!= "Not Applicable")

  # save coeff. of variation data w/ BMP Type
  #write.csv(cv_BMPType, "Coefficient of Variation for Parent Datasets with Metadata.csv", row.names = FALSE)

  # histogram plot
  # calculate binwidth using Freedman–Diaconis Rule (good for skewed data)
  fd_binwidth <- function(x) {
    2 * IQR(x, na.rm = TRUE) / length(x)^(1/3)
  }
  
  # calculate sample sizes for labeling
  counts_df <- cv_BMPType_plt %>%
    group_by(BMP_Category_long) %>%
    summarise(n = n(), .groups = "drop") %>%
    mutate(label = paste0("n = ", n))
  
  hist_by_BMPcat <- ggplot(cv_BMPType_plt, aes(x = coefficientOfVariation, fill = BMP_Category_long)) + 
    geom_histogram(color = "black", alpha = 0.5, position = "identity", binwidth = fd_binwidth) +
    facet_wrap(~ BMP_Category_long, scales = "free_y") +
    geom_text(data = counts_df,
              aes(x = Inf, y = Inf, label = label),
              hjust = 1.3, vjust = 1.6,
              inherit.aes = FALSE,
              size = 5) +
    labs(x = "Coefficient of Variation",
         y = "Count",
         fill = "BMP Category") +
    theme_bw(base_size = 20)
  
  hist_by_BMPcat
  
  # save as png
  ggsave(here::here("Plots", "BMP_Category_histograms_FD_binwidth.png"), width = 2600, height = 1000, units = "px", dpi = 150)

  # Kruskal-Wallis test for significant differences b/w groups (non-parametric)
  kw_BMPcat <- kruskal.test(coefficientOfVariation ~ BMP_Category, data = cv_BMPType_plt)
  kw_BMPcat # p = 1.273e-06 << 0.05 --> significant difference b/w groups
  
  # Posthoc: Dunn's Test w/ Holm correction for p-values
  library(FSA)
  
  dt_BMPcat <- dunnTest(coefficientOfVariation ~ BMP_Category, 
                    data = cv_BMPType_plt,
                    method = "holm")
  
  dunn_result <- dt_BMPcat$res
  
  # save as CSV
  #write.csv(dunn_result, "Dunn's Test Results_Differences in COV by BMP Category.csv", row.names = FALSE)
  
  
  