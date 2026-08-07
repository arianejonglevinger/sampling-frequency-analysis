##################################################
## Author: Matt McGauley
## Email: mmcgau01@villanova.edu
##################################################
import random
random.seed(999)

import os
import glob
import csv
import re
import shutil
import numpy as np
import pandas as pd

exportPath = os.path.join(os.getcwd(), "Export")
if not os.path.exists(exportPath):
    os.makedirs(exportPath)

dataPath = os.path.join(os.getcwd(), "Data")
if not os.path.exists(dataPath):
    os.makedirs(dataPath)

locationsDataPath = os.path.join(dataPath, "LocationsData")
if not os.path.exists(locationsDataPath):
    os.makedirs(locationsDataPath)

zoneDataPath = os.path.join(dataPath, "ZoneData")
if not os.path.exists(zoneDataPath):
    os.makedirs(zoneDataPath)

def get_msnames():
    desiredDirs = ["1", "2", "3", "4", "5", "6", "7", "8", "9"]

    subdirs = glob.glob(os.path.join(zoneDataPath, "*"))

    in_msnames = []
    out_msnames = []

    for dir in subdirs:
        dirName = os.path.basename(dir)

        if dirName in desiredDirs:
            files = glob.glob(os.path.join(dir, "*.xlsx"))
            for file in files:
                print(file)
                fileName = os.path.basename(file).split(".")[0]

                copperSheet = None
                tssSheet = None
                phosphorusSheet = None
                keySheet = "key"

                key = pd.read_excel(file, sheet_name = keySheet, engine="openpyxl").rename(columns=str.lower)

                key["value"] = key["value"].str.lower()

                copperSheet = np.where(key.loc[key["value"] == "copper, dissolved"]["index"].empty, copperSheet, key.loc[key["value"] == "copper, dissolved"]["index"])
                phosphorusSheet = np.where(key.loc[key["value"] == "phosphorus as p, total"]["index"].empty, phosphorusSheet, key.loc[key["value"] == "phosphorus as p, total"]["index"])
                tssSheet = np.where(key.loc[key["value"] == "total suspended solids"]["index"].empty, tssSheet, key.loc[key["value"] == "total suspended solids"]["index"])

                if phosphorusSheet.size != 0:
                    phosphorus = pd.read_excel(file, sheet_name = str(phosphorusSheet[0])).rename(columns=str.lower)
                    in_msnames.extend(phosphorus[phosphorus["mstype"].str.contains("in", flags = re.IGNORECASE)]["msname"].unique())
                    out_msnames.extend(phosphorus[phosphorus["mstype"].str.contains("out", flags = re.IGNORECASE)]["msname"].unique())

                if tssSheet.size != 0:
                    tss = pd.read_excel(file, sheet_name = str(tssSheet[0])).rename(columns=str.lower)
                    in_msnames.extend(tss[tss["mstype"].str.contains("in", flags = re.IGNORECASE)]["msname"].unique())
                    out_msnames.extend(tss[tss["mstype"].str.contains("out", flags = re.IGNORECASE)]["msname"].unique())

                if copperSheet.size != 0:
                    copper = pd.read_excel(file, sheet_name = str(copperSheet[0])).rename(columns=str.lower)
                    in_msnames.extend(copper[copper["mstype"].str.contains("in", flags = re.IGNORECASE)]["msname"].unique())
                    out_msnames.extend(copper[copper["mstype"].str.contains("out", flags = re.IGNORECASE)]["msname"].unique())

        else:
            continue

    in_df = pd.DataFrame(list(set(in_msnames)), columns = ["InflowName"])
    in_df.to_csv(os.path.join(locationsDataPath, "in_msnames.csv"), index = False)

    out_df = pd.DataFrame(list(set(out_msnames)), columns = ["OutflowName"])
    out_df.to_csv(os.path.join(locationsDataPath, "out_msnames.csv"), index = False)

get_msnames()

shutil.make_archive(os.path.join(os.getcwd(), exportPath, "LocationsData"), "zip", locationsDataPath)
shutil.make_archive(os.path.join(os.getcwd(), exportPath, "ZoneData"), "zip", zoneDataPath)




main_folder_path = zoneDataPath

# Create an empty list to store the data
data = []

# Walk through the main folder
for subfolder in os.listdir(main_folder_path):
    subfolder_path = os.path.join(main_folder_path, subfolder)
    if os.path.isdir(subfolder_path):  # Check if it's a directory
        for file_name in os.listdir(subfolder_path):
            file_path = os.path.join(subfolder_path, file_name)
            if os.path.isfile(file_path):  # Check if it's a file
                data.append([subfolder, file_name.strip(".xlsx")])

# Create a DataFrame from the collected data
df = pd.DataFrame(data, columns=["Rain_Zone", "Location"])

# Display the DataFrame
print(df)

# Save to a CSV if needed
df.to_csv("rain_zone_mapping.csv", index=False)