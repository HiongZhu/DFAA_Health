# DFAA_Health

## 1. System Requirements
To run this code, you will need the following environment:
*   **Operating System:** [e.g., Windows 10/11, macOS 12+, or Ubuntu 20.04]
*   **Software:** R version [e.g., 4.5.0] or higher
*   **R Package Dependencies:** 
    *   `dlnm` (Version [e.g., 2.4.7])
    *   `mgcv` (Version [e.g., 1.8-40])
    *   `splines`
    *   `dplyr`, `data.table`

## 2. Installation Guide
1. Install R from [CRAN](https://cran.r-project.org/).
2. Open R or RStudio and install the required packages by running the following command in the console:
   ```R
   install.packages(c("dlnm", "mgcv", "dplyr", "ggplot2", "data.table"))

## 3. Files
The .csv file demonstrates the data format used in the two-stage time-series analysis (hospitalization counts are simulated demo values).

two-stage-DLNM.R: two-stage time-series analysis.

DFAA_definiton.R: definition of drought-to-flood events.
