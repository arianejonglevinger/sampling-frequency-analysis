##################################################
## Author: Matt McGauley
## Email: mmcgau01@villanova.edu
## Edited by: Ariane Jong-Levinger
##### Modified: 8/7/2026
##################################################

# Reads raw water-quality Excel workbooks (per rain zone / site) from Data/ZoneData,
# splits each into inflow/outflow EMC (event mean concentration) datasets per analyte,
# writes per-location metadata + bootstrap-resampled datasets, and rolls everything up
# into summary CSVs used by the downstream R analysis scripts in Summary_Analysis/.

import glob
import os
import random
import shutil

import numpy as np
import pandas as pd
from scipy.stats import variation

random.seed(999)  # fixed seed so bootstrap resampling is reproducible across runs

# Set up (and create, if missing) the folder structure the rest of the script reads from / writes to
dataPath = os.path.join(os.getcwd(), "Data")
if not os.path.exists(dataPath):
    os.makedirs(dataPath)

locationsDataPath = os.path.join(dataPath, "LocationsData")
if not os.path.exists(locationsDataPath):
    os.makedirs(locationsDataPath)

zoneDataPath = os.path.join(dataPath, "ZoneData")
if not os.path.exists(zoneDataPath):
    os.makedirs(zoneDataPath)

outputPath = os.path.join(dataPath, "OutputData")
if not os.path.exists(outputPath):
    os.makedirs(outputPath)

resampledPath = os.path.join(dataPath, "ResampledData")
if not os.path.exists(resampledPath):
    os.makedirs(resampledPath)

metadataPath = os.path.join(dataPath, "MetaData")
if not os.path.exists(metadataPath):
    os.makedirs(metadataPath)

exportPath = os.path.join(os.getcwd(), "Export")
if not os.path.exists(exportPath):
    os.makedirs(exportPath)


# Module-level accumulator dataframes, populated across all locations/analytes as
# createDataframes() runs, then written out once at the end of the script
median_detection_limits = pd.DataFrame(columns = ["Location", "Analyte", "Flow", "median_detection_limit"])
parent_quartiles = pd.DataFrame(columns = ["Location", "Analyte", "Flow", "NumSamples", "Quartile", "EMC"])

# Categorizes an EMC distribution's skew using standard skewness magnitude thresholds
# (|skew| > 1 = heavy, 0.5-1 = moderate, < 0.5 = negligible)
def getSkew(output_df):
    skewValue = output_df["EMC"].skew()
    if abs(skewValue) > 1:
        if skewValue < 0:
            return "heavy_left_skew"
        if skewValue > 0:
            return "heavy_right_skew"
    if abs(skewValue) < 1 and abs(skewValue) >= 0.5:
        if skewValue < 0:
            return "moderate_left_skew"
        if skewValue > 0:
            return "moderate_right_skew"
    if abs(skewValue) < 0.5:
        return "no_skew"

# Average number of distinct calendar months sampled per year, across all years in the raw data
def getNumMonths(original):
    calc_average = original.copy()
    calc_average["datesample"] = pd.to_datetime(calc_average["datesample"])
    calc_average["month"] = calc_average["datesample"].dt.month
    calc_average["year"] = calc_average["datesample"].dt.year

    return(int(calc_average.groupby(calc_average["year"])["month"].nunique().mean()))

# Ratio of max to min EMC value; a simple measure of how much concentrations vary within a dataset
def getInduction(output_df):
    induction = max(output_df["EMC"])/min(output_df["EMC"])
    return induction

# Computes summary metadata (sample counts, skew, coefficient of variation, detection-limit info,
# etc.) for one location/analyte/flow combination and writes it to its own CSV under Data/MetaData.
# Also appends this combination's median detection limit to the module-level accumulator.
def getMetaData(output_df, original, fileName, dirName, analyte, flow, RainZone):
    global median_detection_limits

    is_below_detect = original["wqqualifier"] == "U"
    is_screening_ok = original["initial analysis screening"] != "No"
    is_category_flag_ok = original["categoryanalysisscreen_flag"].isin(["Y", "="])

    below_detect = original[is_below_detect]["value_subhalfdl"] * 2
    mdl = below_detect.median()

    median_detection_limits = pd.concat([median_detection_limits, pd.DataFrame([[fileName, analyte, flow, mdl]], columns = ["Location", "Analyte", "Flow", "median_detection_limit"])])

    # Number_of_Obs_near_dl can only be computed when a median detection limit exists
    if pd.isna(mdl):
        metaData_df = pd.DataFrame({"Number_of_raw_Observations": len(original),
                                    "Number_of_below_detect_Observations": len(original[is_below_detect & is_screening_ok & is_category_flag_ok]),
                                    "Number_of_excluded_Observations": len(original[~is_screening_ok & is_category_flag_ok]),
                                    "Number_of_final_Observations": len(output_df),
                                    "Number_of_Obs_near_dl": np.nan,
                                    "mdl": mdl,
                                    "skewnessCategory": getSkew(output_df),
                                    "skewnessValue": output_df["EMC"].skew(),
                                    "coefficientOfVariation": variation(output_df["EMC"]),
                                    "averageAnnualMonthsSampled": getNumMonths(original),
                                    "inductionRatio": getInduction(output_df),
                                    "distributionCategory": None,
                                    "RainZone": RainZone },
                                   index=[0])
    else:
        metaData_df = pd.DataFrame({"Number_of_raw_Observations": len(original),
                                    "Number_of_below_detect_Observations": len(original[is_below_detect & is_screening_ok & is_category_flag_ok]),
                                    "Number_of_excluded_Observations": len(original[~is_screening_ok & is_category_flag_ok]),
                                    "Number_of_final_Observations": len(output_df),
                                    "Number_of_Obs_near_dl": len(output_df[(output_df["EMC"] >= mdl - (0.5 * mdl)) & (output_df["EMC"] <= mdl + (0.5 * mdl))]),
                                    "mdl": mdl,
                                    "skewnessCategory": getSkew(output_df),
                                    "skewnessValue": output_df["EMC"].skew(),
                                    "coefficientOfVariation": variation(output_df["EMC"]),
                                    "averageAnnualMonthsSampled": getNumMonths(original),
                                    "inductionRatio": getInduction(output_df),
                                    "distributionCategory": None,
                                    "RainZone": RainZone },
                                   index=[0])

    LocAnalyteMetaData = os.path.join(metadataPath, fileName + "_" + analyte)
    if not os.path.exists(LocAnalyteMetaData):
        os.makedirs(LocAnalyteMetaData)

    metaData_df.to_csv(os.path.join(LocAnalyteMetaData, dirName + "_"
                                    + fileName + "_" + flow + "_" + analyte.lower() + ".csv"), index=False)


# Bootstrap-resamples a parent EMC dataset at increasing sample sizes (5 up to the largest
# tier the dataset qualifies for) and writes 500 draws per sample size to its own CSV, for
# use in the downstream sampling-frequency uncertainty analysis.
def resampleData(df, type, dirName, fileName, analyte):
    numObs = len(df)

    # sample-size tiers scale with how many observations are available
    if numObs == 15:
        sampleSizes = [5, 10]
    if numObs <= 20 and numObs > 15:
        sampleSizes = [5, 10, 15]
    if numObs <= 25 and numObs > 20:
        sampleSizes = [5, 10, 15, 20]
    if numObs > 25:
        sampleSizes = [5, 10, 15, 20, 25]

    ResampledLocAnalyte = os.path.join(resampledPath, fileName + "_" + analyte)
    if not os.path.exists(ResampledLocAnalyte):
        os.makedirs(ResampledLocAnalyte)

    for sampleSize in sampleSizes:
        resampled = pd.DataFrame()
        for i in range(1, 501):  # 500 independent random draws per sample size
            sample = df.sample(n=sampleSize)
            sample["Run"] = i
            resampled = pd.concat([resampled, sample])

        resampled.to_csv(os.path.join(ResampledLocAnalyte, dirName + "_" + fileName + "_"
                                      + type + "_" + analyte.lower() + "_resampled_" + str(sampleSize) + ".csv"), index=False)


# Builds the inflow EMC dataset for one location/analyte (filtering to the given inflow
# monitoring-station names and excluding records flagged out by initial analysis screening 
# and category analysis screening). If the result qualifies as a parent dataset (>= 20 events),
#  writes metadata + quantiles + the resampled datasets; datasets that are empty or too small are skipped.
def getInflow(df, dirName, fileName, analyte, names, RainZone):
    global parent_quartiles

    is_screening_ok = df["initial analysis screening"] != "No"
    is_category_flag_ok = df["categoryanalysisscreen_flag"].isin(["Y", "="])

    original_inflow = df[df["msname"].isin(names)]
    inflow = df[(df["msname"].isin(names)) & is_screening_ok & is_category_flag_ok]
    
    if inflow.empty:
        return
    else:
        output_df = pd.DataFrame()
        output_df["EMC"] = inflow["value_subhalfdl"]
        output_df["Location"] = fileName
        output_df["Analyte"] = analyte
        output_df["Flow"] = "Inflow"

        output_df = output_df.dropna()
        output_df = output_df[output_df["EMC"] > 0]
        
        # if data set has >= 20 events, consider it a parent data set and calculate quantiles
        if len(output_df) < 20:
            return
        else:
            getMetaData(output_df, original_inflow,
                        fileName, dirName, analyte, "inflow", RainZone)
            LocAnalyte = os.path.join(outputPath, fileName + "_" + analyte)
            if not os.path.exists(LocAnalyte):
                os.makedirs(LocAnalyte)

            Quartiles = [0.1, 0.25, 0.5, 0.75, 0.9]
            quartile_values = output_df["EMC"].quantile(Quartiles).reset_index().rename(columns = {"index" : "Quartile"})
            quartile_values["NumSamples"] = len(output_df)
            quartile_values["Flow"] = "Inflow"
            quartile_values["Analyte"] = analyte
            quartile_values["Location"] = fileName

            parent_quartiles = pd.concat([parent_quartiles, quartile_values])

            output_df.to_csv(os.path.join(LocAnalyte, dirName + "_"
                                          + fileName + "_inflow_" + analyte.lower() + ".csv"), index=False)
            resampleData(output_df, "inflow", dirName, fileName, analyte)


# Mirrors getInflow() above, but for outflow monitoring stations
def getOutflow(df, dirName, fileName, analyte, names, RainZone):
    global parent_quartiles
    
    is_screening_ok = df["initial analysis screening"] != "No"
    is_category_flag_ok = df["categoryanalysisscreen_flag"].isin(["Y", "="])

    original_outflow = df[df["msname"].isin(names)]
    outflow = df[(df["msname"].isin(names)) & is_screening_ok & is_category_flag_ok]

    if outflow.empty:
        return
    else:
        output_df = pd.DataFrame()
        output_df["EMC"] = outflow["value_subhalfdl"]
        output_df["Location"] = fileName
        output_df["Analyte"] = analyte
        output_df["Flow"] = "Outflow"

        output_df = output_df.dropna()
        output_df = output_df[output_df["EMC"] > 0]

        # if data set has >= 20 events, consider it a parent data set and calculate quantiles
        if len(output_df) < 20:
            return
        else:
            getMetaData(output_df, original_outflow,
                        fileName, dirName, analyte, "outflow", RainZone)
            LocAnalyte = os.path.join(outputPath, fileName + "_" + analyte)
            if not os.path.exists(LocAnalyte):
                os.makedirs(LocAnalyte)

            Quartiles = [0.1, 0.25, 0.5, 0.75, 0.9]
            quartile_values = output_df["EMC"].quantile(Quartiles).reset_index().rename(columns = {"index" : "Quartile"})
            quartile_values["NumSamples"] = len(output_df)
            quartile_values["Flow"] = "Outflow"
            quartile_values["Analyte"] = analyte
            quartile_values["Location"] = fileName

            parent_quartiles = pd.concat([parent_quartiles, quartile_values])

            output_df.to_csv(os.path.join(LocAnalyte, dirName + "_"
                                          + fileName + "_outflow_" + analyte.lower() + ".csv"), index=False)
            resampleData(output_df, "outflow", dirName, fileName, analyte)


# Main driver: walks every rain-zone sub-directory under Data/ZoneData, reads each site's
# Excel workbook, extracts the Copper/Phosphorus/TSS sheets (if present) and runs them
# through getInflow()/getOutflow(), then rolls up per-location inflow counts and BMP
# category/flow-type combinations into summary CSVs.
def createDataframes():

    print("createDataframes() started")

    # all raw data stored in Data > ZoneData in the following sub-directories
    desiredDirs = ["1", "2", "3", "4", "5", "6", "7", "8", "9"]

    data = os.path.join(dataPath, "ZoneData")
    subdirs = glob.glob(os.path.join(data, "*"))

    # monitoring-station names that identify inflow vs outflow records within each workbook
    in_names = pd.read_csv(os.path.join(locationsDataPath, "in_msnames.csv"), names=[
                           "InflowName"])["InflowName"].to_list()
    out_names = pd.read_csv(os.path.join(locationsDataPath, "out_msnames.csv"), names=[
                            "OutflowName"])["OutflowName"].to_list()

    # summary accumulators, rolled up across all rain zones/locations/analytes
    numInflows = pd.DataFrame(columns = ["Location", "Analyte", "Rain_Zone", "numInflow"])
    allCats = pd.DataFrame(columns = ["Location", "Analyte", "Rain_Zone", "Flow_Type", "BMP_Category"])

    for dir in subdirs:
        # 4/23/26: change dirName to work across different OS's
            # original written for Unix-style paths
        #dirName = dir.split("/")[-1]
        dirName = os.path.basename(dir)

        if dirName in desiredDirs:

            files = glob.glob(os.path.join(dir, "*.xlsx"))

            for file in files:

                fileName = os.path.basename(file).split(".")[0]

                copperSheet = None
                tssSheet = None
                phosphorusSheet = None
                keySheet = "key"

                # the "key" sheet maps each analyte name to the sheet index containing its data
                key = pd.read_excel(file, sheet_name=keySheet).rename(columns=str.lower)

                key["value"] = key["value"].str.lower()

                # look up which sheet (if any) holds each analyte's data in this workbook
                copperSheet = np.where(key.loc[key["value"] == "copper, dissolved"]["index"].empty,
                                       copperSheet, key.loc[key["value"] == "copper, dissolved"]["index"])
                phosphorusSheet = np.where(key.loc[key["value"] == "phosphorus as p, total"]["index"].empty,
                                           phosphorusSheet, key.loc[key["value"] == "phosphorus as p, total"]["index"])
                tssSheet = np.where(key.loc[key["value"] == "total suspended solids"]["index"].empty,
                                    tssSheet, key.loc[key["value"] == "total suspended solids"]["index"])

                RainZone = dirName

                # For each analyte present in this workbook: process inflow/outflow into
                # per-location datasets, then (if the location has >= 20 screened-in events)
                # record its inflow station count and its BMP category for each flow direction
                # it qualifies on. The TSS and Copper blocks below repeat this same pattern.
                if phosphorusSheet.size != 0:
                    phosphorus = pd.read_excel(file, sheet_name=str(phosphorusSheet[0])).rename(columns=str.lower)
                    getInflow(phosphorus, dirName, fileName,
                              "Phosphorus", in_names, RainZone)
                    getOutflow(phosphorus, dirName, fileName,
                               "Phosphorus", out_names, RainZone)

                    phosphorus_screening_ok = phosphorus["initial analysis screening"] != "No"
                    phosphorus_category_flag_ok = phosphorus["categoryanalysisscreen_flag"].isin(["Y", "="])
                    phosphorus_in = phosphorus[phosphorus["msname"].isin(in_names) & phosphorus_screening_ok & phosphorus_category_flag_ok]
                    phosphorus_out = phosphorus[phosphorus["msname"].isin(out_names) & phosphorus_screening_ok & phosphorus_category_flag_ok]

                    if len(phosphorus_in) >= 20:
                        num_inflow = len(phosphorus_in["msname"].unique())
                        numInflows = pd.concat([numInflows, pd.DataFrame([[fileName, "Phosphorus", dirName, num_inflow]], columns = ["Location", "Analyte", "Rain_Zone", "numInflow"])])

                    if len(phosphorus_in) >= 20:
                        cat = phosphorus["bmpcategory_code"].iloc[0]
                        allCats = pd.concat([allCats, pd.DataFrame([[fileName, "Phosphorus", "Inflow", dirName, cat]], columns = ["Location", "Analyte", "Rain_Zone", "Flow_Type", "BMP_Category"])])

                    if len(phosphorus_out) >= 20:
                        cat = phosphorus["bmpcategory_code"].iloc[0]
                        allCats = pd.concat([allCats, pd.DataFrame([[fileName, "Phosphorus", "Outflow", dirName, cat]], columns = ["Location", "Analyte", "Rain_Zone", "Flow_Type", "BMP_Category"])])

                # TSS: same processing pattern as Phosphorus above
                if tssSheet.size != 0:
                    tss = pd.read_excel(file, sheet_name=str(tssSheet[0])).rename(columns=str.lower)
                    getInflow(tss, dirName, fileName, "TSS", in_names, RainZone)
                    getOutflow(tss, dirName, fileName, "TSS", out_names, RainZone)

                    tss_screening_ok = tss["initial analysis screening"] != "No"
                    tss_category_flag_ok = tss["categoryanalysisscreen_flag"].isin(["Y", "="])
                    tss_in = tss[tss["msname"].isin(in_names) & tss_screening_ok & tss_category_flag_ok]
                    tss_out = tss[tss["msname"].isin(out_names) & tss_screening_ok & tss_category_flag_ok]

                    if len(tss_in) >= 20:
                        num_inflow = len(tss_in["msname"].unique())
                        numInflows = pd.concat([numInflows, pd.DataFrame([[fileName, "TSS", dirName, num_inflow]], columns = ["Location", "Analyte", "Rain_Zone", "numInflow"])])

                    if len(tss_in) >= 20:
                        cat = tss["bmpcategory_code"].iloc[0]
                        allCats = pd.concat([allCats, pd.DataFrame([[fileName, "TSS", "Inflow", dirName, cat]], columns = ["Location", "Analyte", "Rain_Zone", "Flow_Type", "BMP_Category"])])

                    if len(tss_out) >= 20:
                        cat = tss["bmpcategory_code"].iloc[0]
                        allCats = pd.concat([allCats, pd.DataFrame([[fileName, "TSS", "Outflow", dirName, cat]], columns = ["Location", "Analyte", "Rain_Zone", "Flow_Type", "BMP_Category"])])


                # Copper: same processing pattern as Phosphorus above
                if copperSheet.size != 0:
                    copper = pd.read_excel(file, sheet_name=str(copperSheet[0])).rename(columns=str.lower)
                    getInflow(copper, dirName, fileName, "Copper", in_names, RainZone)
                    getOutflow(copper, dirName, fileName, "Copper", out_names, RainZone)

                    copper_screening_ok = copper["initial analysis screening"] != "No"
                    copper_category_flag_ok = copper["categoryanalysisscreen_flag"].isin(["Y", "="])
                    copper_in = copper[copper["msname"].isin(in_names) & copper_screening_ok & copper_category_flag_ok]
                    copper_out = copper[copper["msname"].isin(out_names) & copper_screening_ok & copper_category_flag_ok]

                    if len(copper_in) >= 20:
                        num_inflow = len(copper_in["msname"].unique())
                        numInflows = pd.concat([numInflows, pd.DataFrame([[fileName, "Copper", dirName, num_inflow]], columns = ["Location", "Analyte", "Rain_Zone", "numInflow"])])

                    if len(copper_in) >= 20:
                        cat = copper["bmpcategory_code"].iloc[0]
                        allCats = pd.concat([allCats, pd.DataFrame([[fileName, "Copper", "Inflow", dirName, cat]], columns = ["Location", "Analyte", "Rain_Zone", "Flow_Type", "BMP_Category"])])

                    if len(copper_out) >= 20:
                        cat = copper["bmpcategory_code"].iloc[0]
                        allCats = pd.concat([allCats, pd.DataFrame([[fileName, "Copper", "Outflow", dirName, cat]], columns = ["Location", "Analyte", "Rain_Zone", "Flow_Type", "BMP_Category"])])
        else:
            continue

    # write the rolled-up summary dataframes once all rain zones/locations have been processed
    parent_quartiles.to_csv(os.path.join(metadataPath, "parent_dataset_quartiles.csv"), index = False)
    numInflows.to_csv(os.path.join(metadataPath, "num_inflow_locations.csv"), index = False)
    allCats.to_csv(os.path.join(metadataPath, "type_by_location_rz.csv"), index = False)

    print("createDataframes() done")


# Run the full extraction pipeline, then write out the median-detection-limit accumulator
# (populated as a side effect of getMetaData() calls inside createDataframes())
createDataframes()
median_detection_limits.to_csv(os.path.join(metadataPath, "median_detection_limits.csv"), index = False)

# Zip each output folder for easy download/sharing via the Export/ directory
shutil.make_archive(os.path.join(os.getcwd(), exportPath, "LocationsData"), "zip", locationsDataPath)
shutil.make_archive(os.path.join(os.getcwd(), exportPath,
                                 "ZoneData"), "zip", zoneDataPath)
shutil.make_archive(os.path.join(os.getcwd(), exportPath,
                                 "OutputData"), "zip", outputPath)
shutil.make_archive(os.path.join(os.getcwd(), exportPath,
                                 "ResampledData"), "zip", resampledPath)
shutil.make_archive(os.path.join(os.getcwd(), exportPath,
                                 "MetaDataData"), "zip", metadataPath)
