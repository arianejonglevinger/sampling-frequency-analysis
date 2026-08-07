##################################################
## Author: Matt McGauley
## Credits: Emily Darin for ggplot code
## Email: mmcgau01@villanova.edu
## Edited by: Ariane Jong-Levinger
##### Modified: 8/7/2026
##################################################
# For each site/analyte/flow combination, combines the original (parent) dataset with its
# Monte Carlo-resampled datasets to compute percentile estimates and relative percent difference
# (RPD) summary statistics, then exports the results as "all_scatter_data.csv" and
# "merged_original_and_resampled_data.csv".
import glob
import os
import random
import shutil
import warnings

import matplotlib.pyplot as plt
import pandas as pd
import plotnine
from matplotlib import gridspec
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import ListedColormap
from matplotlib.lines import Line2D
from plotnine import *

random.seed(999)

warnings.filterwarnings("ignore")

# --- Output directories ---
plotPath = os.path.join(os.getcwd(), "Plots")
if not os.path.exists(plotPath):
    os.makedirs(plotPath)

dataPath = os.path.join(os.getcwd(), "Data")
if not os.path.exists(plotPath):
    os.makedirs(plotPath)

# Holds intermediate per-location RPD/percentile summary CSVs that feed the scatter-data
# export below, not plots.
siteScatterDataPath = os.path.join(dataPath, "SiteScatterData")
if not os.path.exists(siteScatterDataPath):
    os.makedirs(siteScatterDataPath)

# --- Input directories, populated by upstream data-prep/resampling scripts ---
outputPath = os.path.join(dataPath, "OutputData")
if not os.path.exists(outputPath):
    os.makedirs(outputPath)

resampledPath = os.path.join(dataPath, "ResampledData")
if not os.path.exists(resampledPath):
    os.makedirs(resampledPath)

metadataPath = os.path.join(dataPath, "MetaData")
if not os.path.exists(metadataPath):
    os.makedirs(metadataPath)

# --- More output directories ---
quartilePlotPath = os.path.join(plotPath, "QuartilePlots")
if not os.path.exists(quartilePlotPath):
    os.makedirs(quartilePlotPath)

scatterDataPath = os.path.join(dataPath, "ScatterData")
if not os.path.exists(scatterDataPath):
    os.makedirs(scatterDataPath)

scatterInflow = os.path.join(scatterDataPath, "Inflow")
if not os.path.exists(scatterInflow):
    os.makedirs(scatterInflow)

scatterOutflow = os.path.join(scatterDataPath, "Outflow")
if not os.path.exists(scatterOutflow):
    os.makedirs(scatterOutflow)

exportPath = os.path.join(os.getcwd(), "Export")
if not os.path.exists(exportPath):
    os.makedirs(exportPath)

# Convenience copies of the final CSV outputs for downstream summary analysis scripts.
summaryAnalysisDataPath = os.path.join(os.getcwd(), "Summary_Analysis", "Data")
if not os.path.exists(summaryAnalysisDataPath):
    os.makedirs(summaryAnalysisDataPath)

# Accumulates original-vs-resampled EMC pairs across every site/analyte/flow processed by
# create_output(); populated as a side effect of rpd() and exported at the end as
# "merged_original_and_resampled_data.csv".
all_merged_resample = []

def rpd(original, resampled):
    """
    Compute the mean relative percent difference (RPD) between each resampled percentile
    estimate and the corresponding original (parent-dataset) percentile value, grouped by
    sample size (N) and percentile.

    As a side effect, also appends the merged original/resampled rows to the module-level
    all_merged_resample list for later export.
    """
    resampled = resampled.rename(columns = {"EMC" : "resampledEMC"})
    for_export = resampled.merge(original, how = "left")
    all_merged_resample.append(for_export[["Run", "N", "Location", "Analyte", "Quartile", "resampledEMC", "EMC", "Flow"]])
    resampled = resampled.merge(original, how = "left", on = "Quartile").drop(columns = ["Location_x", "Analyte_x", "Flow_x", "Location_y", "Analyte_y", "Flow_y"])
    resampled["rpd"] = resampled.apply(lambda x: abs(x["resampledEMC"] - x["EMC"])/((x["resampledEMC"] + x["EMC"])/2), axis = 1)
    return resampled.groupby(["N", "Quartile"])["rpd"].agg("mean")

def create_output():
    """
    Main pipeline. For every original (parent) dataset file, find its matching resampled
    Monte Carlo draws and metadata, compute per-N/percentile RPD summary statistics, and
    accumulate them into per-analyte scatter dataframes. Exports the combined result as
    "all_scatter_data.csv" and the raw original/resampled pairs as
    "merged_original_and_resampled_data.csv".
    """
    # Collect every original/resampled/metadata CSV across all site subdirectories.
    all_original_files = []
    subDirs = glob.glob(os.path.join(outputPath, "*"))
    for dir in subDirs:
        csv_files = glob.glob(os.path.join(dir, "*.csv"))
        all_original_files.extend(csv_files)

    all_resampled_files = []
    resampledSubDirs = glob.glob(os.path.join(resampledPath, "*"))
    for dir in resampledSubDirs:
        csv_files = glob.glob(os.path.join(dir, "*.csv"))
        all_resampled_files.extend(csv_files)

    all_metadata_files = []
    metadataSubDirs = glob.glob(os.path.join(metadataPath, "*"))
    for dir in metadataSubDirs:
        csv_files = glob.glob(os.path.join(dir, "*.csv"))
        all_metadata_files.extend(csv_files)

    # Accumulators for the per-analyte scatter data, populated inside the loop below.
    Copper_ScatterData = pd.DataFrame(columns = ["N", "Quartile", "median", "std", "rpd", "Location", "Flow", "Analyte", "NumObservations", "NearMedianDetectionLimit", "NumObservationsBelowDetect","changeRPD", "percentOriginal", "skewnessCategory", "skewnessValue", "coefficientOfVariation", "averageAnnualMonthsSampled", "inductionRatio", "distributionCategory", "RainZone"])
    TSS_ScatterData = pd.DataFrame(columns = ["N", "Quartile", "median", "std", "rpd", "Location", "Flow", "Analyte", "NumObservations", "NearMedianDetectionLimit", "NumObservationsBelowDetect","changeRPD", "percentOriginal", "skewnessCategory", "skewnessValue", "coefficientOfVariation", "averageAnnualMonthsSampled", "inductionRatio", "distributionCategory", "RainZone"])
    Phosphorus_ScatterData = pd.DataFrame(columns = ["N", "Quartile", "median", "std", "rpd", "Location", "Flow", "Analyte", "NumObservations", "NearMedianDetectionLimit", "NumObservationsBelowDetect", "changeRPD", "percentOriginal", "skewnessCategory", "skewnessValue", "coefficientOfVariation", "averageAnnualMonthsSampled", "inductionRatio",  "distributionCategory", "RainZone"])

    for dataFile in all_original_files:
        # Match this original file to its resampled draws and metadata by shared filename.
        fileName = os.path.basename(dataFile).split(".")[0]
        resampledDataFiles = [resampledDataFile for resampledDataFile in all_resampled_files if fileName in resampledDataFile]
        metadataFiles = [metadataDataFile for metadataDataFile in all_metadata_files if fileName in metadataDataFile]

        originalDataframe = pd.read_csv(dataFile)
        resampledDataframes = [pd.read_csv(file) for file in resampledDataFiles]
        metadataDataframe = [pd.read_csv(file) for file in metadataFiles][0]

        # Each resampled file contains 500 Monte Carlo runs at one fixed sample size; label
        # every row with its sample-size category ("n=5", "n=10", etc.) before combining them.
        for frame in resampledDataframes:
            frame["N"] =  "n=" + str(int(len(frame)/500))
            frame["N"] = frame["N"].astype("category")

        fullResampledDataframe = pd.concat(resampledDataframes, ignore_index = True)
        fullResampledDataframe["N"] = fullResampledDataframe["N"].astype("category")
        values = fullResampledDataframe["N"].unique()

        # Not every site was resampled at every sample size, so the category order (and the
        # plot color palette) must be built dynamically based on which "n=" values are present.
        if "n=20" not in values and "n=25" not in values and "n=15" not in values:
            fullResampledDataframe["N"] = fullResampledDataframe["N"].cat.reorder_categories(["n=5", "n=10"])
            color_list = ["#2c925f", "#2e4057", "#000000"]
        if "n=20" not in values and "n=25" not in values and "n=15" in values:
            fullResampledDataframe["N"] = fullResampledDataframe["N"].cat.reorder_categories(["n=5", "n=10", "n=15"])
            color_list = ["#2c925f", "#2e4057", "#810f7c", "#000000"]
        if "n=20" in values and "n=25" not in values:
            fullResampledDataframe["N"] = fullResampledDataframe["N"].cat.reorder_categories(["n=5", "n=10", "n=15", "n=20"])
            color_list = ["#2c925f", "#2e4057", "#810f7c", "#de2d26", "#000000"]
        if "n=20" in values and "n=25" in values:
            fullResampledDataframe["N"] = fullResampledDataframe["N"].cat.reorder_categories(["n=5", "n=10", "n=15", "n=20", "n=25"])
            color_list = ["#2c925f", "#2e4057", "#810f7c", "#de2d26", "#41b6c4", "#000000"]

        Quartiles = [0.1, 0.25, 0.5, 0.75, 0.9]

        # Quartile EMC estimates per Monte Carlo run (resampled) vs. the single "true" percentile
        # from the full parent dataset (original) — the latter is what RPD is measured against.
        resampledQuartiles = fullResampledDataframe.groupby(["Run", "N", "Location", "Analyte", "Flow"])["EMC"].quantile(Quartiles).reset_index().rename(columns = {"level_5" : "Quartile"})
        originalQuartiles = originalDataframe.groupby(["Location", "Analyte", "Flow"])["EMC"].quantile(Quartiles).reset_index().rename(columns = {"level_3" : "Quartile"})
        originalQuartiles["resampled_vline_label"] = "Quartile\nEMC"
        originalQuartiles["original_vline_label"] = "Parent\ndataset\nmedian\ndetection\nlimit"

        # Per-N/percentile summary: median & std across Monte Carlo runs, plus mean RPD (as a
        # percentage) against the parent-dataset percentile.
        siteScatterData = resampledQuartiles.groupby(["N", "Quartile"])["EMC"].agg(["median", "std"])
        siteScatterData["rpd"] = rpd(originalQuartiles, resampledQuartiles) * 100
        siteScatterData["Location"] = originalQuartiles["Location"].iloc[0]
        siteScatterData["Flow"] = originalQuartiles["Flow"].iloc[0]
        siteScatterData["Analyte"] = originalQuartiles["Analyte"].iloc[0]
        siteScatterData["NumObservations"] = metadataDataframe["Number_of_final_Observations"].iloc[0]
        siteScatterData["NumObservationsBelowDetect"] = metadataDataframe["Number_of_below_detect_Observations"].iloc[0]
        siteScatterData["NearMedianDetectionLimit"] = metadataDataframe["Number_of_Obs_near_dl"].iloc[0]
        siteScatterData["skewnessCategory"] = metadataDataframe["skewnessCategory"].iloc[0]
        siteScatterData["distributionCategory"] = metadataDataframe["distributionCategory"].iloc[0]
        siteScatterData["skewnessValue"] = metadataDataframe["skewnessValue"].iloc[0]
        siteScatterData["coefficientOfVariation"] = metadataDataframe["coefficientOfVariation"].iloc[0]
        siteScatterData["averageAnnualMonthsSampled"] = metadataDataframe["averageAnnualMonthsSampled"].iloc[0]
        siteScatterData["inductionRatio"] = metadataDataframe["inductionRatio"].iloc[0]
        siteScatterData["RainZone"] = metadataDataframe["RainZone"].iloc[0]
        siteScatterData = siteScatterData.reset_index()
        # Change in RPD between consecutive sample sizes (within the same percentile), and the
        # sample size as a percentage of the total observations available for this site.
        siteScatterData["changeRPD"] = siteScatterData.groupby(["Quartile"])["rpd"].transform(lambda row: abs(row - row.shift())).fillna(0)
        siteScatterData["percentOriginal"] = siteScatterData.apply(lambda row: (int(row["N"].split("=")[1])/row["NumObservations"]) * 100, axis = 1)

        originalQuartiles["mdl"] = metadataDataframe["mdl"].iloc[0]

        # Write this site's summary as its own intermediate CSV (not read back elsewhere in
        # this script — kept for ad hoc inspection/debugging).
        siteScatterDataDestination = os.path.join(siteScatterDataPath, originalQuartiles["Location"].iloc[0])
        if not os.path.exists(siteScatterDataDestination):
            os.makedirs(siteScatterDataDestination)

        siteScatterData.to_csv(os.path.join(siteScatterDataDestination, originalQuartiles["Location"].iloc[0] + "_" + originalQuartiles["Flow"].iloc[0] + "_" + originalQuartiles["Analyte"].iloc[0] + ".csv"))

        # Route this site's summary into the scatter dataframe for its analyte.
        if originalQuartiles["Analyte"].iloc[0] == "Copper":
            Copper_ScatterData = pd.concat([Copper_ScatterData, siteScatterData])
        if originalQuartiles["Analyte"].iloc[0] == "Phosphorus":
            Phosphorus_ScatterData = pd.concat([Phosphorus_ScatterData, siteScatterData])
        if originalQuartiles["Analyte"].iloc[0] == "TSS":
            TSS_ScatterData = pd.concat([TSS_ScatterData, siteScatterData])

        LocAnalyte = os.path.join(quartilePlotPath, originalQuartiles["Location"].iloc[0] + "_" + originalQuartiles["Analyte"].iloc[0])
        if not os.path.exists(LocAnalyte):
            os.makedirs(LocAnalyte)

        # Legacy per-site percentile density/scatter PDF plots — disabled, kept for reference.
        # with PdfPages(os.path.join(LocAnalyte, originalQuartiles["Location"].iloc[0] + "_" + originalQuartiles["Flow"].iloc[0] + "_" + originalQuartiles["Analyte"].iloc[0] + ".pdf")) as pdf:
        #     for sample in resampledQuartiles["Quartile"].unique():
        #         if originalQuartiles["Analyte"].iloc[0] == "Copper":
        #             xlabel = "EMC (µg/L)"
        #         else:
        #             xlabel = "EMC (mg/L)"

        #         originalQuartiles["Analyte"].replace({
        #                 "Copper" : "Dissolved Copper",
        #                 "Phosphorus" : "Total Phosphorus",
        #                 "TSS" : "Total Suspended Solids"
        #                 }, inplace = True)

        #         finalPlot, ax1 = plt.subplots()
        #         left, bottom, width, height = [0.25, 0.6, 0.2, 0.2]
        #         # ax2 = finalPlot.add_axes([left, bottom, width, height])

        #         thisPlot = (ggplot(resampledQuartiles[resampledQuartiles["Quartile"] == sample], plotnine.aes(x = "EMC", label = "Location", color = "N")) +
        #                         geom_line(stat = "density", size = 2) +
        #                         geom_vline(originalQuartiles[originalQuartiles["Quartile"] == sample], plotnine.aes(xintercept = "EMC", color = "resampled_vline_label"), linetype = "dashed", size = 2, show_legend = {"color" : True}) +
        #                         scale_color_manual(values = color_list) +
        #                         guides(color = guide_legend(title = "Number of\nMonitoring\nEvents\n\n\n")) +
        #                         ggtitle(originalQuartiles["Location"].iloc[0] + "\n" + originalQuartiles["Flow"].iloc[0] + "\n" + "Quartile " + str(int(sample * 100)) + "\n" + originalQuartiles["Analyte"].iloc[0])+
        #                         xlab(xlabel) +
        #                         expand_limits(x = 0) +
        #                         theme_light() +
        #                         theme(legend_position = (0.95, .75), axis_title_y = element_blank(), axis_text_y = element_blank(), text = element_text(size = 15, family = "serif"))
        #                     )


        #         scatter = (ggplot(originalDataframe, plotnine.aes(x = "EMC", y = originalDataframe.index)) +
        #                         geom_point() +
        #                         geom_vline(originalQuartiles[originalQuartiles["Quartile"] == sample], plotnine.aes(xintercept = "mdl", color = "original_vline_label"), linetype = "dashed", size = 2, show_legend = {"color" : True}) +
        #                         scale_color_manual(values = ["#FF0000"]) +
        #                         expand_limits(x = 0) +
        #                         guides(color = guide_legend(title = "")) +
        #                         theme_light() +
        #                         theme(axis_title_y = element_blank(), axis_text_y = element_blank(), text = element_text(size = 15, family = "serif"))
        #                     )

        #         # fig = (ggplot() + theme_void() + theme(figure_size = (60, 40))).draw()

        #         # gs = gridspec.GridSpec(40, 40)
        #         # ax1 = fig.add_subplot(gs[:,:])

        #         # _ = thisPlot._draw_using_figure(fig, [ax1])
                
        #         # 4/23/26
        #         fig = thisPlot.draw()

        #         # ax1.text(0.3 * (ax1.get_xlim()[1]), 1.01 * (ax1.get_ylim()[1]),
        #         # "Number of raw observations: " + str(metadataDataframe["Number_of_raw_Observations"].iloc[0]) + "\n" +
        #         # "Number of below-detection observations: " + str(metadataDataframe["Number_of_below_detect_Observations"].iloc[0]) + "\n" +
        #         # "Number of excluded values: " + str(metadataDataframe["Number_of_excluded_Observations"].iloc[0]) + "\n" +
        #         # "Number of final observations: " + str(metadataDataframe["Number_of_final_Observations"].iloc[0]),
        #         # size = 15, fontfamily = "serif")

        #         # ax1.set_title(originalQuartiles["Location"].iloc[0] + "\n" + originalQuartiles["Flow"].iloc[0] + "\n" + "Quartile " + str(int(sample * 100)) + "\n" + originalQuartiles["Analyte"].iloc[0], fontsize = 30, fontfamily = "serif", loc = "left")
        #         # ax1.set_xlabel(xlabel, size = 51, fontfamily = "serif")

        #         # if originalQuartiles[originalQuartiles["Quartile"] == sample]["EMC"].iloc[0]/ax1.get_xlim()[1] > 0.5:
        #         #     ax2 = fig.add_subplot(gs[2:10:, 2:10])
        #         # elif originalQuartiles[originalQuartiles["Quartile"] == sample]["EMC"].iloc[0]/ax1.get_xlim()[1] < 0.5:
        #         #     ax2 = fig.add_subplot(gs[2:10, -10:-2])
        #         # else:
        #         #     ax2 = fig.add_subplot(gs[2:10:, 2:10])

        #         # 4/23/26
        #         #_ = scatter._draw_using_figure(fig, [ax2])
        #         scatter.draw()

        #         # ax2.set_title("Original Data Distribution", fontsize = 5, fontfamily = "serif")
        #         # ax2.set_xlabel(xlabel, size = 15, fontfamily = "serif")

        #         pdf.savefig(fig)
        #     plt.close("all")

    # Per-analyte scatter data, split by flow direction.
    Copper_ScatterData[Copper_ScatterData["Flow"] == "Inflow"].to_csv(os.path.join(scatterInflow, "Copper_inflow.csv"), index = False)
    TSS_ScatterData[TSS_ScatterData["Flow"] == "Inflow"].to_csv(os.path.join(scatterInflow, "TSS_inflow.csv"), index = False)
    Phosphorus_ScatterData[Phosphorus_ScatterData["Flow"] == "Inflow"].to_csv(os.path.join(scatterInflow, "Phosphorus_inflow.csv"), index = False)

    Copper_ScatterData[Copper_ScatterData["Flow"] == "Outflow"].to_csv(os.path.join(scatterOutflow, "Copper_outflow.csv"), index = False)
    TSS_ScatterData[TSS_ScatterData["Flow"] == "Outflow"].to_csv(os.path.join(scatterOutflow, "TSS_outflow.csv"), index = False)
    Phosphorus_ScatterData[Phosphorus_ScatterData["Flow"] == "Outflow"].to_csv(os.path.join(scatterOutflow, "Phosphorus_outflow.csv"), index = False)

    # Primary output #1: all analytes combined into a single scatter data CSV.
    combined_ScatterData = pd.concat([Copper_ScatterData, TSS_ScatterData, Phosphorus_ScatterData])
    combined_ScatterData.to_csv(os.path.join(scatterDataPath, "all_scatter_data.csv"), index = False)
    combined_ScatterData.to_csv(os.path.join(summaryAnalysisDataPath, "all_scatter_data.csv"), index = False)

    # Primary output #2: raw original-vs-resampled EMC pairs, populated by rpd() above.
    # Only written if rpd() actually ran (i.e. at least one site was processed).
    if all_merged_resample:
        merged_original_and_resampled_data = pd.concat(all_merged_resample)
        merged_original_and_resampled_data.to_csv(
            os.path.join(scatterDataPath, "merged_original_and_resampled_data.csv"),
            index=False
        )
        merged_original_and_resampled_data.to_csv(
            os.path.join(summaryAnalysisDataPath, "merged_original_and_resampled_data.csv"),
            index=False
        )
    else:
        print("No merged original/resampled data available — skipping merged export.")


create_output()

# Zip the plot/data output folders for easy download/handoff.
shutil.make_archive(os.path.join(os.getcwd(), exportPath, "QuartilepPlots"), "zip", quartilePlotPath)
shutil.make_archive(os.path.join(os.getcwd(), exportPath, "ScatterData"), "zip", scatterDataPath)
