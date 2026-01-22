# Modelling-Extreme-Precipitation-in-Aotearoa-New-Zealand-with-Bootstrapping-Regression-and-BART

This repository contains R code corresponding to the work 'Modelling Extreme Precipitation in Aotearoa (New Zealand) from 1951–2025 with Bootstrapping, Negative Binomial Regression, and Bayesian Additive Regression Trees' by Floe Foxon. 

## Analyses

All code required to run the analyses in this work are contained within the 'Analysis-Code.R' file located in this repository. No generative AI was used for this code.

## Data

All data required to run these analyses are publicly available from the following sources.

Precipitation data were sourced from the European Centre for Medium-Range Weather Forecasts (ECMWF) Reanalysis v5 (ERA5) product 'ERA5 post-processed daily-statistics on single levels from 1940 to present' from the Copernicus Climate Data Store at https://doi.org/10.24381/cds.4991cf48 using the following request parameters: 
Variable: Total precipitation; Year: 1951, 1952, ..., 2025; Month: January, February, ..., December; Day: 01, 02, ..., 31; Daily statistic: Daily mean; Time zone: UTC+00:00; Frequency: 6-hourly; Geographical area: North: -34°, West: 166°, South: -48°, East: 179°. Data were downloaded for each year as individual netCDF files (.nc) named after each year, e.g. '1951.nc'.

Average annual Aotearoa New Zealand Temperature data for the period 1951-2018 were sourced from the National Institute of Water and Atmospheric Research (NIWA) product "'Seven-station' series temperature data" at https://niwa.co.nz/climate-and-weather/nz-temperature-record/seven-station-series-temperature-data The data are provided when clicking 'The adjusted data' tab and 'NIWA 'seven-station' temperature series: annual data for mean temperature [XLS 23 KB]' file at the link above, with the '7-Stn Composite Temp' column taken and manually placed into a separate CSV file (.csv). Average annual Aotearoa New Zealand Temperature data for the period 2019-2025 were sourced from NIWA's Annual Climate Summaries for those years at https://niwa.co.nz/climate-and-weather/annual/annual-climate-summary-2019; https://niwa.co.nz/climate-and-weather/annual/annual-climate-summary-2020; https://niwa.co.nz/climate-and-weather/annual/annual-climate-summary-2021; https://niwa.co.nz/climate-and-weather/annual/annual-climate-summary-2022; https://niwa.co.nz/climate-and-weather/annual/annual-climate-summary-2023; https://niwa.co.nz/climate-and-weather/annual/annual-climate-summary-2024; https://niwa.co.nz/climate-and-weather/annual/annual-climate-summary-2025 under the 'Overview' sections and were placed in the same CSV file as the 1951-2018 data, which is referred to as 'Temperature_Data.csv' in the code.

Oceanic Nino Index (ONI) data were sourced from NOAA's National Weather Service Climate Prediction Center product 'Historical El Nino / La Nina episodes (1950-present) - Cold & Warm Episodes by Season' at https://web.archive.org/web/20260104051947/https://www.cpc.ncep.noaa.gov/products/analysis_monitoring/ensostuff/ONI_v5.php The table provided at the link above was manually placed into a CSV file (.csv), which is referred to as 'ONI_Data.csv' in the code.
