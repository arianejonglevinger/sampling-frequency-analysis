# sampling-frequency-analysis

The purpose of this repository is to answer the question: how many storm events 
need to be sampled to achieve a given level of uncertainty in event-mean 
concentration (EMC) for a given pollutant? It does this by resampling large water 
quality datasets and determining the relative uncertainty in EMCs between the full
dataset and smaller subsamples of that dataset.

## Description

The scripts in this repository filter EMC data from the 2024 version of the  
International Stormwater BMP Database (https://bmpdatabase.org/urban) for 
datasets with sample size >= 20 ("parent datasets"), resample the parent datasets 
500 times to create hypothetical smaller monitoring programs ("subsamples") with 
sample sizes of 5, 10, 15, 20, and 25 events, and calculate the Relative Percent
Difference (RPD) between EMCs in the parent datasets and their subsamples as a 
measure of uncertainty. Program managers can compare uncertainty in EMC 
percentiles for different pollutants across hypothetical program sample sizes 
to select a sample size that appropriately balances program resources with 
statistical confidence in the resulting EMC monitoring data.

### A Note About the Data

The data in this repository were extracted from the International Stormwater BMP
Database (2024 version), which is a Microsoft Access Database. An existing query
in the database called "qry2WQAnalysisFlatFile" was modified to filter for the 
relevant data for this analysis. The custom query (written in SQL) is shown
below:
```
SELECT *
FROM qry2WQAnalysisFlatFile
WHERE (((qry2WQAnalysisFlatFile.[ParameterName]) In ('Copper, Dissolved','Total suspended solids','Phosphorus AS P, Total','Phosphorus AS PO4, Total')) AND ((qry2WQAnalysisFlatFile.[BMPID]) In (SELECT BMPID 
    FROM qry2WQAnalysisFlatFile
    WHERE MSType IN ('Inflow', 'Outflow') 
      AND ParameterName IN ('Copper, Dissolved', 'Total suspended solids', 'Phosphorus AS P, Total', 'Phosphorus AS PO4, Total')
      AND SampleType = 'EMC-Flow Weighted'
    GROUP BY BMPID, ParameterName
    HAVING COUNT(IIF(MSType = 'Inflow', 1, NULL)) >= 20
       OR COUNT(IIF(MSType = 'Outflow', 1, NULL)) >= 20
  )));
```
The data are saved in Excel spreadsheets in subfolders corresponding to their 
Rain Zone in the directory ~/Data/ZoneData.

Data from the following Location IDs were excluded from the analysis following
discussion with the American Society of Civil Engineers/Environmental and Water 
Resources Institute BMP Database Committee because they did not meet criteria
for inclusion:
702665323_TX
-1762711317_WA
-43160592_WA
488336797_FL
58597072_TX
1841653012_CA
-2055029627_CA
-357056449_CA
-1166835566_CO

Data are filtered for acceptable BMP Categories in the script **`FinalFigures.R`**

## Getting Started

### Dependencies

#### Python 
Requires **Python 3.11 or later** (this floor comes from `pandas` 3.x's own minimum
supported version, not from anything in these scripts' own code).

**`getData.py`**
- `pandas` (tested: 3.0.0)
- `numpy` (tested: 2.4.2)
- `scipy` (tested: 1.17.0) — specifically `scipy.stats.variation`
- `openpyxl` (tested: 3.1.5) — not imported directly, but required by `pandas.read_excel()`
  as the engine for reading `.xlsx` files; without it, Excel reads fail even though
  `openpyxl` never appears in the script

**`getPlots.py`**
- `pandas` (tested: 3.0.0)
- `matplotlib` (tested: 3.10.8) and `plotnine` (tested: 0.15.3) — still imported and
  required for the script to run without an `ImportError`, but not currently exercised:
  the only code that calls into them (per-site quartile density/scatter plots) is commented
  out. If that plotting code is removed for good rather than re-enabled, these two
  dependencies can be dropped from this list.

#### R
**`FinalFigures.R`** requires **R 4.5.1+**, run via a standalone system R installation
(separate from any R that might be installed elsewhere, e.g. via conda for other tools).
Required packages — `tidyverse`, `ggplot2`, `here`, `ggpubr`, `RColorBrewer`, `ggtext`,
`moments` — are all auto-installed by the script itself on first run if missing.

### Installing

* All program files are available for download at: https://github.com/arianejonglevinger/sampling-frequency-analysis

### Executing program

* Step 1: Run getData.py
* Step 2: Run getPlots.py
* Step 3: Run FinalFigures.R

## Troubleshooting

Please do not alter the folder structure or the structure of the data, as this 
may result in errors when running the scripts.

## Authors

* Original Author: Matthew McGauley
* Editor: Ariane Jong-Levinger (arijl@sccwrp.org)

## Version History

* 1.0
    * Initial Release

## License

This project is licensed under the GNU Affero General Public License (AGPL) v3 License - 
see the COPYING.txt and the COPYRIGHT.txt files for details.

## Acknowledgments

This work was partially made possible with support from the Villanova Center for 
Resilient Water Systems and the Southern California Stormwater Monitoring Coalition.