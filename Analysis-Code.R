#### INSTALL/IMPORT PACKAGES ####

#install.packages('ncdf4') # install ncdf4 package
library(ncdf4) # for importing netCDF data
#install.packages('clusterSim') # install clusterSim package
library(clusterSim) # for calculating Davies-Bouldin index (k-means)
#install.packages('sfsmisc') # install sfsmisc package
library(sfsmisc) # for scientific notation on axes
#install.packages('lattice') # install lattice package
library(lattice) # for plotting strip chart
#install.packages('mgcv') # install mgcv package
library(mgcv) # for generalised additive model (GAM)
#install.packages('MASS') # install MASS package
library(MASS) # for negative binomial model
#install.packages('dplyr') # install dplyr package
library(dplyr) # for test/train split (validation)
#install.packages('caret') # install caret package
library(caret) # for test/train split (validation)
#install.packages('segmented') # install segmented package
library(segmented) # for segmented model
#install.packages('BART') # install BART package
library(BART) # for BART model
library(scales) # for line opacity
#install.packages('tseries') # install tseries package
library(tseries) # for calculating autocorrelation
#install.packages('bayestestR') # install bayestestR package
library(bayestestR) # for calculating CIs



#### DATA SOURCES ####

# Precipitation data were sourced from the European Centre for Medium-Range 
# Weather Forecasts (ECMWF) Reanalysis v5 (ERA5) product 'ERA5 post-
# processed daily-statistics on single levels from 1940 to present' from the 
# Copernicus Climate Data Store at https://doi.org/10.24381/cds.4991cf48 
# using the following request parameters: 
# Variable: Total precipitation
# Year: 1951, 1952, ..., 2025
# Month: January, February, ..., December
# Day: 01, 02, ..., 31
# Daily statistic: Daily mean
# Time zone: UTC+00:00
# Frequency: 6-hourly
# Geographical area: North: -34°, West: 166°, South: -48°, East: 179°.
# Data were downloaded for each year as individual netCDF files (.nc) named
# after each year, e.g. '1951.nc'.
#
# Average annual Aotearoa New Zealand Temperature data for the period 1951-2018
# were sourced from the National Institute of Water and Atmospheric Research
# (NIWA) product "'Seven-station' series temperature data" at 
# https://niwa.co.nz/climate-and-weather/nz-temperature-record/seven-station-series-temperature-data
# The data are provided when clicking 'The adjusted data' tab and 
# 'NIWA 'seven-station' temperature series: annual data for mean temperature [XLS 23 KB]'
# file at the link above, with the '7-Stn Composite Temp' column taken and 
# manually placed into a separate CSV file (.csv).
# Average annual Aotearoa New Zealand Temperature data for the period 2019-2025
# were sourced from NIWA's Annual Climate Summaries for those years at
# https://niwa.co.nz/climate-and-weather/annual/annual-climate-summary-2019
# https://niwa.co.nz/climate-and-weather/annual/annual-climate-summary-2020
# https://niwa.co.nz/climate-and-weather/annual/annual-climate-summary-2021
# https://niwa.co.nz/climate-and-weather/annual/annual-climate-summary-2022
# https://niwa.co.nz/climate-and-weather/annual/annual-climate-summary-2023
# https://niwa.co.nz/climate-and-weather/annual/annual-climate-summary-2024
# https://niwa.co.nz/climate-and-weather/annual/annual-climate-summary-2025
# under the 'Overview' sections and were placed in the same CSV file as the
# 1951-2018 data.
#
# Oceanic Nino Index (ONI) data were sourced from NOAA's National Weather
# Service Climate Prediction Center product 'Historical El Nino / La Nina 
# episodes (1950-present) - Cold & Warm Episodes by Season' at
# https://web.archive.org/web/20260104051947/https://www.cpc.ncep.noaa.gov/products/analysis_monitoring/ensostuff/ONI_v5.php
# The table provided at the link above was manually placed into a CSV file 
# (.csv).



#### IMPORT AND PROCESS BASELINE (1951-1980) DATA ####

# create data frame for storing baseline daily precipitation data
dataFrameBase = data.frame(matrix(ncol = 3, nrow = 0))
colnames(dataFrameBase) = c('Year', 'Day', 'Precipitation')
# for each year of data from 1951-1980 (baseline)...
for (year in 1951:1980){
  # get next data file name (e.g., '1951.nc')
  fileName = paste(as.character(year),'.nc', sep = '')
  # open data file at file name
  dataFile = nc_open(fileName)
  # extract data from file
  data = ncvar_get(dataFile)
  # for each day of data...
  for (day in 1:365){
    # average data across space (latitude/longitude), i.e. zero-dimensional
    dataFrameBase[nrow(dataFrameBase)+1,] = c(year, day, mean(data[,,day]))
  }
}

# delete variables that are no longer needed
remove(fileName)
remove(data)
remove(dataFile)
remove(year)
remove(day)



#### EXTREME VALUE/ANOMALY DETECTION ####

# 99TH-PERCENTILE METHOD

# This method just uses the 99th-percentile of precipitation as the threshold 
# for defining extreme/anomalous precipitation.

# get value of extreme/anomaly threshold
anomalyPrecipitationThreshold = quantile(dataFrameBase$Precipitation, 0.99)[[1]]
# sort baseline precipitation data, ascending
precipitationBaseSorted = sort(dataFrameBase$Precipitation, decreasing = FALSE)
# get values of extreme/anomalous precipitation based on threshold
baselineAnomalies =
precipitationBaseSorted[which(precipitationBaseSorted > 
                              anomalyPrecipitationThreshold)]
# count number of extreme/anomalous precipitation observations at baseline
numBaselineAnomalies = length(baselineAnomalies)

# ALTERNATIVE K-MEANS METHOD

# This method applies the k-means algorithm to the precipitation data with the
# optimal value of k determined by the Davies-Bouldin index method. The final 
# non-overlapping cluster (i.e., the one with the largest precipitation values)
# defines extreme/anomalous precipitation; all values in this cluster are 
# extreme/anomalous, with the smallest value in the cluster defining the
# extreme/anomaly threshold.

# set random number seed for reproducibility
set.seed(123)
# create list for storing Davies-Bouldin index values across values of k
DBs = c()
# for each value of k from 2 to 10...
for (k in 2:10){
  # fit k-means clustering model to data (with k clusters)
  kMeansModel = kmeans(dataFrameBase$Precipitation, centers = k)
  # calculate Davies-Bouldin index for this value of k
  nextDB = index.DB(dataFrameBase$Precipitation, kMeansModel$cluster)$DB
  # add this to the list of Davies-Bouldin index values
  DBs = c(DBs, nextDB)
}
# find optimal value of k (the value that minimizes the Davies-Bouldin index)
optimalk = which.min(DBs)+1
# plot Davies-Bouldin index values for each value of k
plot(c(2:10), DBs, type = 'l', ylab = 'Davies-Bouldin Index', 
     xlab = 'Number of Clusters, k', pch = 20, col = 'grey', lwd = 2, 
     family = 'serif', main = '', cex.main = 1.6)
# add optimal k line
abline(v = optimalk, lty = 2, lwd = 2, col = 'black')

# fit k-means clustering model to data with optimal number of clusters
kMeansModel = kmeans(dataFrameBase$Precipitation, centers = optimalk)
# add cluster column to baseline precipitation data frame
dataFrameBase$Cluster = kMeansModel$cluster
# create list for storing minimum precipitation values of each cluster (edge 
# cases)
clusterMins = c()
# create list for storing number of observations in each cluster
clusterLengths = c()
# for each cluster...
for (c in 1:optimalk){
  # calculate minimum precipitation value in this cluster
  nextMin = min(dataFrameBase$Precipitation[dataFrameBase$Cluster == c])
  # add this to the list of cluster minimum precipitation values
  clusterMins = c(clusterMins, nextMin)
  # count number of observations in this cluster
  nextLength = length(dataFrameBase$Precipitation[dataFrameBase$Cluster == c])
  # add this to list of cluster lengths
  clusterLengths = c(clusterLengths, nextLength)
}
# the cluster indices that come from the kmeans() function are not in order of
# increasing precipitation, so it is necessary to create a data frame for the 
# clusters and sort it by precipitation value
kMeansResultsFrame = data.frame(c <- c(1:7), clusterMin <- clusterMins,
                                clusterLength = clusterLengths)
# sort clusters data frame by precipitation value
kMeansResultsFrame = kMeansResultsFrame[order(kMeansResultsFrame$clusterMin),]
# get minimum precipitation value (precipitation threshold) from final cluster 
altAnomalyPrecipitationThreshold =
kMeansResultsFrame$clusterMin[length(kMeansResultsFrame$clusterMin)]
# get size/length of final cluster (number of extreme values/anomalies)
altNumBaselineAnomalies =
kMeansResultsFrame$clusterLength[length(kMeansResultsFrame$clusterLength)]

# COMPARE 99TH-PERCENTILE and K-MEANS METHODS

cat('The thresholds defining an extreme/anomalous value for total', 
    'precipitation were', anomalyPrecipitationThreshold, 'for the',
    '99th-percentile method and', altAnomalyPrecipitationThreshold, 'for the',
    'k-means method. The numbers of extreme/anomalous precipitation days at', 
    'baseline were', numBaselineAnomalies, 'for the 99th-percentile method',
    'and', altNumBaselineAnomalies, 'for the k-means method.')

# delete variables that are no longer needed
remove(DBs)
remove(k)
remove(nextDB)
remove(optimalk)
remove(kMeansModel)
remove(c)
remove(clusterMins)
remove(clusterMin)
remove(nextMin)
remove(clusterLengths)
remove(nextLength)
remove(kMeansResultsFrame)



#### IMPORT AND PROCESS PRESENT (1996-2025) DATA ####

# create data frame for storing present daily precipitation data
dataFramePres = data.frame(matrix(ncol = 3, nrow = 0))
colnames(dataFramePres) = c('Year', 'Day', 'Precipitation')
# for each year of data from 1996-2025 (present)...
for (year in 1996:2025){
  # get next data file name (e.g., '1996.nc')
  fileName = paste(as.character(year),'.nc', sep = '')
  # open data file at file name
  dataFile = nc_open(fileName)
  # extract data from file
  data = ncvar_get(dataFile)
  # for each day of data...
  for (day in 1:365){
    # average data across space (latitude/longitude), i.e. zero-dimensional
    dataFramePres[nrow(dataFramePres)+1,] = c(year, day, mean(data[,,day]))
  }
}

# delete variables that are no longer needed
remove(fileName)
remove(data)
remove(dataFile)
remove(year)
remove(day)



### APPLY EXTREME/ANOMALY THRESHOLD TO PRESENT (1996-2025) DATA ####

# sort present precipitation data, ascending
precipitationPresSorted = sort(dataFramePres$Precipitation, decreasing = FALSE)
# get values of extreme/anomalous precipitation based on threshold
presAnomalies = 
precipitationPresSorted[which(precipitationPresSorted > 
                              anomalyPrecipitationThreshold)]
# count number of extreme/anomalous precipitation observations at present
numPresAnomalies = length(presAnomalies)

cat('The number of extreme/anomalous zero-dimensional daily total',
    'precipitation events at baseline is', numBaselineAnomalies, '. The number',
    'of extreme/anomalous zero-dimensional daily total precipitation events at',
    'present is', numPresAnomalies)



### MANN-WHITNEY U TEST FOR DIFFERENCE IN PRECIPITATION DISTRIBUTION

# Mann-Whitney U test for difference in distribution (Null hypothesis: baseline 
# and present precipitation data have the same distribution):
wilcox.test(precipitationBaseSorted, precipitationPresSorted)



#### VISUALISE BASELINE AND PRESENT PRECIPITATION ####

# FIGURE 1 (HISTOGRAMS):

# plot histogram of baseline precipitation
hist(precipitationBaseSorted, col=NULL, border = 'dodgerblue4', 
     ylab = 'Density', 
     xlab = 'Zero-Dimensional Daily Mean Total Precipitation (m)', lwd = 1, 
     xlim = c(0, 7e-4), ylim = c(0, 6500), freq = FALSE, main = '', xaxt = 'n',
     family = 'serif', cex.lab = 1.3)
# fit gamma distribution to baseline precipitation data
alphaHatBaseline = (mean(precipitationBaseSorted)/
                    sd(precipitationBaseSorted))^2
lambdaHatBaseline = mean(precipitationBaseSorted)/
                    (sd(precipitationBaseSorted)^2)
gamFitBaseline = dgamma(precipitationBaseSorted, alphaHatBaseline, 
                        lambdaHatBaseline)
# add gamma fit line to plot
lines(precipitationBaseSorted, gamFitBaseline, type = 'l', col = 'dodgerblue4')
# plot histogram of present precipitation
hist(precipitationPresSorted, col=NULL, border = 'cyan3', lwd = 1, 
     freq = FALSE, add = TRUE, breaks = 17) # breaks chosen to match baseline
# fit gamma distribution to present precipitation data
alphaHatPres = (mean(precipitationPresSorted)/sd(precipitationPresSorted))^2
lambdaHatPres = mean(precipitationPresSorted)/(sd(precipitationPresSorted)^2)
gamFitPres = dgamma(precipitationPresSorted, alphaHatPres, lambdaHatPres)
# add gamma fit line to plot
lines(precipitationPresSorted, gamFitPres, type = 'l', col = 'cyan3')
# add tick labels with scientific notation (x-axis only)
eaxis(1, cex.axis = 1.1, family = 'serif')
# set font
op = par(family = 'serif')
# add group means to plot as vertical lines
abline(v = mean(precipitationBaseSorted), lty = 2, lwd = 1, col = 'dodgerblue4')
abline(v = mean(precipitationPresSorted), lty = 2, lwd = 1, col = 'cyan3')
# add legend to plot
legend('topright', lty = c(1, 1, 2), 
       col = c('dodgerblue4', 'cyan3', 'cyan3'), 
       legend=c('Baseline (1951–1980)',
                'Present   (1996–2025)',
                'Averages'),
       pt.cex = 0.5, cex = 1.2, seg.len = c(1.3, 1.3, 1.3), bty = 'n', 
       x.intersp = 0.4)

# FIGURE 2 (BOX PLOTS AND STRIP CHARTS):

# create plot object
plot(1, type = 'n', 
     xlab = 'Zero-Dimensional Daily Mean Total Precipitation (m)', 
     cex.lab = 1.3, ylab = '', ylim = c(0.45, 2.15), 
     xlim=c(0, 0.00115), xaxt = 'n', yaxt = 'n', family = 'serif', 
     main = '', cex.main = 1.6)
# plot box plot of baseline precipitation
boxplot(precipitationBaseSorted, horizontal = TRUE, outline = FALSE, 
        outwex = 0.5, boxwex = 0.5, add = TRUE, xaxt = 'n', col = 'white', 
        border = 'dodgerblue4')
# plot strip chart of baseline precipitation
stripchart(precipitationBaseSorted, method='jitter', pch = 19, cex = 0.3, 
           col = 'dodgerblue4', add = TRUE, at = 1.7, xaxt = 'n')
# plot box plot of present precipitation
boxplot(precipitationPresSorted, horizontal = TRUE, outline = FALSE, 
        outwex = 0.5, boxwex = 0.5, add = TRUE, xaxt = 'n', col = 'white', 
        border = 'cyan3', at = 0.65)
# plot strip chart of present precipitation
stripchart(precipitationPresSorted, method='jitter', pch = 19, cex = 0.3, 
           col = 'cyan3', add = TRUE, at = 1.35, xaxt = 'n')
# add 99th-percentile method threshold to plot as vertical line
abline(v = anomalyPrecipitationThreshold, lty = 2, lwd = 2, col = 'black')
# add k-means method threshold to plot as vertical line
abline(v = altAnomalyPrecipitationThreshold, lty = 2, lwd = 2, col = 'darkgrey')
# add text to plot indicating extreme and non-extreme regions
text(2.5e-4, 2.05, 'Non-Extreme', family = 'serif', cex = 1.1)
# add text to plot indicating extreme and non-extreme regions
text(8.7e-4, 2.05, 'Extreme', family = 'serif', cex = 1.1)
# add tick labels with scientific notation (x-axis only)
eaxis(1, cex.axis = 1.1, family = 'serif')
# set font
op = par(family = 'serif')
# add legend to plot
legend('bottomright', pch = c(19, 19, NA, NA), lty = c(NA, NA, 2, 2),
       lwd = c(NA, NA, 2, 2),
       col = c('dodgerblue4', 'cyan3', 'black', 'darkgrey'), 
       legend=c('Baseline (1951–1980)', 
                'Present   (1996–2025)', 
                '99th-Percentile Threshold',
                'K-Means Threshold'),
       pt.cex = 0.5, cex = 1.1, seg.len = c(1.5, 1.5, 1.5), bty = 'n', 
       x.intersp = 0.4)



#### CALCULATE BASELINE/PRESENT EXRTREME/ANOMALOUS PRECIP. PROBABILITY ####

# This method uses the gamma distributions that were fitted to the baseline 
# precipitation data and the present precipitation data respectively (above) to 
# calculate the probability of observing a precipitation value at least as large 
# as the extreme/anomaly threshold.

# Pr(precipitation >= threshold) from baseline gamma CDF
probAnomalyBaseline = pgamma(anomalyPrecipitationThreshold, alphaHatBaseline, 
                             lambdaHatBaseline, lower.tail=F)
# Pr(precipitation >= threshold) from present gamma CDF
probAnomalyPres = pgamma(anomalyPrecipitationThreshold, alphaHatPres, 
                         lambdaHatPres, lower.tail=F)

cat('The probability of observing an extreme/anomalous zero-dimensional daily',
    'total precipitation at baseline is', probAnomalyBaseline, 
    '. The probability of observing an extreme/anomalous zero-dimensional',
    'daily total precipitation at present is', probAnomalyPres)

# delete variables that are no longer needed
remove(alphaHatBaseline)
remove(lambdaHatBaseline)
remove(alphaHatPres)
remove(lambdaHatPres)
remove(gamFitBaseline)
remove(gamFitPres)
remove(probAnomalyBaseline)
remove(probAnomalyPres)



#### BOOTSTRAP INFERENCE ####

# This method re-samples from the baseline precipitation data to generate a
# prediction interval. This interval describes the range of values for the 
# number of days of extreme/anomalous precipitation for which there is a 95% 
# chance that the number of days of extreme/anomalous precipitation for a future 
# 30-year window will fall within the range, given the baseline distribution of 
# precipitation data.

# set number of bootstrap sample sets
B = 1e5
# create list for storing bootstrap statistics
bootstrapTs = c() 
# for each desired bootstrap sample set...
for (k in 1:B){
  # generate bootstrap sample set
  nextBootstrapSampleSet = sample(precipitationBaseSorted, 
                                  length(precipitationBaseSorted), 
                                  replace = TRUE)
  # sort bootstrap sample set
  nextBootstrapSampleSetSorted = sort(nextBootstrapSampleSet, 
                                      decreasing = FALSE)
  # find num. of extreme/anomalous precipitation events in bootstrap sample set
  nextBootstrapT = length(nextBootstrapSampleSetSorted[which(
    nextBootstrapSampleSetSorted > anomalyPrecipitationThreshold)])
  # add this to list of bootstrap statistics 
  bootstrapTs = c(bootstrapTs, nextBootstrapT)
}
# set alpha level for (1-alpha)*100% prediction interval
alpha = 0.05
# calculate prediction interval (bootstrap statistic quantiles)
T_UL = quantile(bootstrapTs, c(alpha/2, 1-alpha/2))
T_L = T_UL[[1]] # lower
T_U = T_UL[[2]] # upper
predictionInterval = paste('(', as.character(T_L), '–', as.character(T_U), ')', 
                           sep = '')
cat('The prediction interval for the number of extreme/anomalous precipitation',
    'events at baseline is', predictionInterval)

# delete variables that are no longer needed
remove(B)
remove(bootstrapTs)
remove(k)
remove(nextBootstrapSampleSet)
remove(nextBootstrapSampleSetSorted)
remove(nextBootstrapT)
remove(alpha)
remove(T_UL)
remove(T_L)
remove(T_U)



#### IMPORT AND PROCESS REMAINING DATA (1981-1995) ####

# create data frame for storing present daily precipitation data
dataFrameMid = data.frame(matrix(ncol = 3, nrow = 0))
colnames(dataFrameMid) = c('Year', 'Day', 'Precipitation')
# for each year of data from 1981-1995 (present)...
for (year in 1981:1995){
  # get next data file name (e.g., '1981.nc')
  fileName = paste(as.character(year),'.nc', sep = '')
  # open data file at file name
  dataFile = nc_open(fileName)
  # extract data from file
  data = ncvar_get(dataFile)
  # for each day of data...
  for (day in 1:365){
    # average data across space (latitude/longitude), i.e. zero-dimensional
    dataFrameMid[nrow(dataFrameMid)+1,] = c(year, day, mean(data[,,day]))
  }
}

# delete variables that are no longer needed
remove(fileName)
remove(data)
remove(dataFile)
remove(year)
remove(day)



#### COMBINE AND PROCESS DATA ACROSS ALL YEARS ####

# for each data set (baseline, present, and between) create a new variable whose
# value is 1 if the day had extreme/anomalous precipitation, and whose value is 
# 0 if the day did not have extreme/anomalous precipitation
dataFrameBase$AnomPrecipitation =  with(dataFrameBase, 
                                    ifelse(Precipitation > 
                                           anomalyPrecipitationThreshold, 1, 0))
dataFrameMid$AnomPrecipitation =  with(dataFrameMid, 
                                    ifelse(Precipitation > 
                                           anomalyPrecipitationThreshold, 1, 0))
dataFramePres$AnomPrecipitation =  with(dataFramePres, 
                                    ifelse(Precipitation > 
                                           anomalyPrecipitationThreshold, 1, 0))
# count the number of extreme/anomalous precipitation days in a given year by 
# summing across the year (annualise the extreme/anomaly count)
dataFrameBaseAnnualised = aggregate(AnomPrecipitation~Year, dataFrameBase, sum)
dataFrameMidAnnualised = aggregate(AnomPrecipitation~Year, dataFrameMid, sum)
dataFramePresAnnualised = aggregate(AnomPrecipitation~Year, dataFramePres, sum)
# combine annualised precipitation data across data sets (1951-2025)
combinedData = rbind(dataFrameBaseAnnualised, dataFrameMidAnnualised, 
                     dataFramePresAnnualised)



#### EXTREME WEATHER TEMPORAL TREND MODELS ####

# These methods fit linear models of annual extreme/anomaly count against time
# to investigate trends in extreme/anomalous precipitation over time.

# SIMPLE CORRELATION

correl = cor(combinedData$Year, combinedData$AnomPrecipitation)
cat('The correlation coefficient between annual extreme/anomaly count and',
    'year is:', correl)

# GENERALISED ADDITIVE TEMPORAL TREND MODEL (POISSON)

# fit generalised additive model with poisson distribution and canonical log 
# link function to data (log-linear model)
fit = gam(AnomPrecipitation ~ s(Year), family = poisson(link = 'log'), 
          data = combinedData)
# view model summary
summary(fit)
# according to the summary output, the edf=1, therefore a spline term is not 
# necessary and a parametric coefficient will suffice.

# GENERALISED LINEAR TEMPORAL TREND MODEL (POISSON)

# fit generalised linear model with poisson distribution and canonical log link 
# function to data (log-linear model)
fit = glm(AnomPrecipitation ~ Year, family = poisson(link = 'log'), 
          data = combinedData)
# view model summary
summary(fit)
# check for overdispersion; the dispersion parameter is estimated by the ratio 
# between the model deviance and the residual degrees of freedom
dispersionHat = fit$deviance/fit$df.residual
cat('The estimate for the dispersion value of the Poisson model is', 
    dispersionHat)
# since the dispersion is >1, there is evidence of overdispersion. Thus, a model
# that accounts for overdispersion should be used. Two options are quasi-Poisson
# and negative binomial. Both are fitted below and compared via the root mean 
# square error (RMSE) predictive criterion in 10-fold cross-validation. Note 
# that traditional model comparison techniques such as AIC, BIC, and 
# likelihood-ratio test are not valid because they depend on a distribution, 
# and quasi-Poisson does not have a distribution. See 
# https://doi.org/10.1890/07-0043.1 for details

# QUASI-POISSON AND NEGATIVE BINOMIAL TEMPORAL TREND MODELS

quasiPoissonTestRMSEs = c() # RMSEs for testing split for quasi-Poisson model
negBinomialTestRMSEs = c() # RMSEs for testing split for negative binomial model
# cross-validation loop. For each loop...
for (j in 1:10){
  # set seed for reproducibility
  set.seed(j)
  # split data into training and testing sets (80-20 split)
  testTrainSplit = combinedData$AnomPrecipitation %>% 
    createDataPartition(p = 0.8, list = FALSE)
  
  # training split
  train = combinedData[testTrainSplit,]
  # fit quasi-Poisson model to training data
  nextQuasiPoissonFit = glm(AnomPrecipitation ~ Year, 
                            family = quasipoisson(link = 'log'), 
                            data = train)
  # fit negative binomial model to training data
  nextNegBinomialFit = glm.nb(AnomPrecipitation ~ Year, 
                              link = log, 
                              data = train)
  
  # testing split
  test = combinedData[-testTrainSplit,]
  # calculate RMSE for quasi-Poisson model for this loop
  quasiPoissonTestRMSEs[j] = sqrt(mean((test$AnomPrecipitation - 
                             predict(nextQuasiPoissonFit, newdata = test))^2))
  # calculate RMSE for negative binomial model for this loop
  negBinomialTestRMSEs[j] = sqrt(mean((test$AnomPrecipitation - 
                            predict(nextNegBinomialFit, newdata = test))^2))
  
}

# calculate mean validated RMSEs for each model
meanValidatedQuasiPoissonTestRMSE = mean(quasiPoissonTestRMSEs)
meanValidatedNegativeBinomialTestRMSE = mean(negBinomialTestRMSEs)
# print results
cat('The quasi-Poisson model mean validated RMSE is', 
    meanValidatedQuasiPoissonTestRMSE, 
    '. The negative binomial model mean validated RMSE is', 
    meanValidatedNegativeBinomialTestRMSE, 
    '. These compare to the range of the outcome:', 
    min(combinedData$AnomPrecipitation), '–', 
    max(combinedData$AnomPrecipitation))
# the two mean validated RMSEs are practically identical. The negative binomial
# model is selected because its mean validated RMSE is technically lower, but 
# only at the fifth decimal place...

# fit negative binomial model to whole dataset
fit = glm.nb(AnomPrecipitation ~ Year, link = log, data = combinedData)
# print negative binomial model summary
summary(fit)
# estimate change in number of days of extreme/anomalous precipitation over time
cat('The annual number of days of extreme/anomalous precipitation increased by', 
    (exp(fit$coefficients[[2]]) - 1)*100, 
    '% (95% CI:',  (exp(confint(fit, level = 0.95)[2, 1]) - 1)*100, 
    '–', (exp(confint(fit, level = 0.95)[2, 2]) - 1)*100, '%)',
    'per year from 1951-2025.')

# NEGATIVE BINOMIAL TEMPORAL TREND MODEL DIAGNOSTICS

# plot diagnostics side-by-side
par(mfrow=c(1,2))
# residuals plot:
# get studentized residuals
r_stu = rstudent(fit)
# Shapiro-Wilk test for normality using studentized residuals (Null hypothesis: 
# residuals are normally distributed):
shapiro.test(r_stu)
# plot residuals vs predicted values
plot(fit$fitted.values, r_stu, type = 'p', xlab = 'Predicted Values',
     ylab = 'Studentized Residuals', pch = 20, col = 'black', ylim=c(-5, 5))
abline(h=3, lty = 'dashed', col = 'gray') # add line indicating large residual
abline(h=0, lty = 'dashed', col = 'gray') # add zero line
abline(h=-3, lty = 'dashed', col = 'gray') # add line indicating large residual
# QQ plot:
qqnorm(r_stu, xlab = 'Theoretical Values', ylab = 'Sample Values', main = '', 
       pch = 20)
qqline(r_stu, col = 'gray', lwd = 2, lty = 2)
# reset to single plots
par(mfrow=c(1,1))
# Leverage and influence:
# number of non-intercept parameters
p = length(fit$coefficients)-1
# sample size
n = length(combinedData$AnomPrecipitation)
# threshold of 'large' leverage values
largeLev = (2*(p+1))/n
# leverage of observations
leverage = hatvalues(fit)
# identify which values have large leverage in the negative binomial model
largeLevIndicator = which(leverage > largeLev)
# threshold of 'large' influence values
largeInf = 1
# influence (Cook's distance) of observations
influence = cooks.distance(fit)
# identify which values have large influence in the negative binomial model
largeInfIndicator = which(influence > largeInf)
cat('There were', length(largeLevIndicator), 'observations with large leverage',
    'which were', leverage[largeLevIndicator], 'compared to the large leverage',
    'threshold of', largeLev, '. There were ', length(largeInfIndicator), 
    'observations with large influence.')

# SEGMENTED NEGATIVE BINOMIAL TEMPORAL TREND MODEL

# To further investigate time trends, fit a piecewise model based on the 
# negative binomial model but with two segments. The starting value for the 
# break-point is 1981, i.e. at the end of the baseline period.
# fit segmented/piecewise model
fitSeg = segmented.glm(fit, seg.Z = ~Year, fixed.psi=c(1981), it.max=0, npsi = 1, 
                       control=seg.control(it.max=0, fix.npsi = TRUE))
# compare AICs for segmented and non-segmented negative binomial models
cat('The segmented model AIC is', fitSeg$aic, 
    '. The non-segmented model AIC is', fit$aic, '.')
# the AIC of the non-segmented model is slightly lower, therefore the segmented 
# model provides no improvement over the non-segmented model given its
# additional complexity

# FIGURE 3 (TEMPORAL TRENDS, NEGATIVE BINOMIAL MODEL):

# get neg. binomial model predictions on link scale with confidence intervals
predictions = predict(fit, newdata = combinedData, type = 'link', 
                      se.fit = TRUE)
# back-transform predictions to response scale and add to combined data frame
combinedData$Pred = exp(predictions$fit)
# set alpha level for (1-alpha)*100% confidence interval
alpha = 0.05
# calculate Wald-type 95% confidence intervals on link scale then back-transform
# to response scale and add to combined data frame
combinedData$PredUpper = exp(predictions$fit + 
                               (qnorm(1-alpha/2) * predictions$se.fit))
combinedData$PredLower = exp(predictions$fit - 
                               (qnorm(1-alpha/2) * predictions$se.fit))
# create plot object
plot(1, type = 'n', xlab = '', ylim = c(0, 13), xlim = c(1949, 2026),
     ylab = 'Days with Extreme 0-D Daily Mean Tot. Precip.', 
     cex.lab = 1.3, family = 'serif', cex.main = 1.6, xaxt = 'n', yaxt = 'n')
# plot data
points(combinedData$Year, combinedData$AnomPrecipitation, pch = 20)
# plot negative binomial fit
lines(combinedData$Year, combinedData$Pred, col = 'blue', lwd = 2)
# plot confidence interval for negative binomial fit
polygon(c(combinedData$Year, rev(combinedData$Year)), 
        c(combinedData$PredUpper, rev(combinedData$PredLower)), 
        col = rgb(0, 0, 1, 0.2), border = FALSE)
# add tick labels with scientific notation (x-axis only)
eaxis(1, cex.axis = 1.3, family = 'serif', at = seq(1950, 2025, 15))
# add tick labels with scientific notation (y-axis only)
eaxis(2, cex.axis = 1.3, family = 'serif')
# set font
op = par(family = 'serif')
# add legend to plot
legend(x = 1945, y = 13.5, lty = c(1), col = c('blue'),
       legend=c('Model Fit'),
       pt.cex = 1.4, cex = 1.2, seg.len = c(1.3), bty = 'n', x.intersp = 0.4)
legend(x = 1945.8, y = 12.6, col = c(rgb(0, 0, 1, 0.2)),
       legend=c('95% Confidence Interval'),
       pt.cex = 1.4, cex = 1.2, fill = c(rgb(0, 0, 1, 0.2)), 
       density = c(200), border = c(rgb(0, 0, 1, 0.2)),
       seg.len = c(1.3), bty = 'n', x.intersp = 0.65)

# delete variables that are no longer needed
remove(dataFrameBase)
remove(dataFrameBaseAnnualised)
remove(dataFrameMid)
remove(dataFrameMidAnnualised)
remove(dataFramePres)
remove(dataFramePresAnnualised)
remove(correl)
remove(fit)
remove(fitSeg)
remove(dispersionHat)
remove(alpha)
remove(p)
remove(influence)
remove(largeInf)
remove(largeInfIndicator)
remove(leverage)
remove(largeLev)
remove(largeLevIndicator)
remove(n)
remove(r_stu)
remove(j)
remove(train)
remove(test)
remove(testTrainSplit)
remove(nextQuasiPoissonFit)
remove(nextNegBinomialFit)
remove(quasiPoissonTestRMSEs)
remove(negBinomialTestRMSEs)
remove(meanValidatedQuasiPoissonTestRMSE)
remove(meanValidatedNegativeBinomialTestRMSE)
remove(predictions)



#### IMPORT AND PROCESS TEMPERATURE AND ONI DATA ####

# import temperature data
warmingData = read.csv('Temperature_Data.csv')
# calculate average annual temperature for the baseline period 1951-1980
baselineTemp = mean(warmingData[warmingData$Year >= 1951 & 
                    warmingData$Year <= 1980, ]$Seven_Stn_Composite_Temp)
# calculate the temperature anomaly
# (annual avg. temp in each year - avg. annual temp. for period 1951-1980)
combinedData$Tanom = warmingData$Seven_Stn_Composite_Temp - baselineTemp
# import ONI data
ONIdata = read.csv('ONI_Data.csv')
# calculate annual average Oceanic Nino Index (ONI) and add to combined data 
# frame
combinedData$ONI = (ONIdata$JFM + ONIdata$AMJ + ONIdata$JAS + ONIdata$OND)/4



#### BAYESIAN ADDITIVE REGRESSION TREE (BART) WARMING ASSOCIATION MODEL ####

# This method investigates the possible association between warming/temperature 
# anomaly and annual extreme/anomalous precipitation count, adjusting for 
# potential confounding by El Nino/La Nina events. A BART model is fitted to the 
# data and used to predict the annual extreme/anomalous precipitation count when 
# the temperature anomaly is 0°C compared to 1.4°C (the maximum temperature
# anomaly observed, which was observed in 2022) across a range of values of ONI.

# FIT BART MODEL

# (this can be slow)
BARTfit = mc.gbart(x.train =  combinedData[c('Tanom', 'ONI')],
                   y.train = combinedData$AnomPrecipitation,
                   type = 'wbart',
                   mc.cores = 5, # number of chains
                   ntree = 1000, # number trees in sum
                   ndpost = 10000, # number of posterior draws (all chains)
                   nskip = 1000, # burn-in
                   seed = 1234)

# BART MODEL DIAGNOSTICS

# trace plot
plot(c(1:length(BARTfit$sigma[,1])), BARTfit$sigma[,1], type = 'l',
     ylab = 'Sigma', xlab = 'Draw', pch = 20, 
     col = alpha('red', 0.4), lwd = 2, family = 'serif', main = '', 
     cex.main = 1.6, ylim = c(1.5, 3.3), xlim = c(0, 3000))
lines(c(1:length(BARTfit$sigma[,2])), BARTfit$sigma[,2], 
      col = alpha('orange', 0.4), lwd = 1)
lines(c(1:length(BARTfit$sigma[,3])), BARTfit$sigma[,3], 
      col = alpha('yellow', 0.4), lwd = 1)
lines(c(1:length(BARTfit$sigma[,4])), BARTfit$sigma[,4], 
      col = alpha('green', 0.4), lwd = 1)
lines(c(1:length(BARTfit$sigma[,5])), BARTfit$sigma[,5], 
      col = alpha('blue', 0.4), lwd = 1)

# calculate autocorrelations
acf(BARTfit$sigma[,1], pl = FALSE)
acf(BARTfit$sigma[,2], pl = FALSE)
acf(BARTfit$sigma[,3], pl = FALSE)
acf(BARTfit$sigma[,4], pl = FALSE)
acf(BARTfit$sigma[,5], pl = FALSE)

# calculate RMSE for BART model
RMSEB = sqrt( sum((BARTfit$yhat.train.mean - combinedData$AnomPrecipitation)^2)/
             length(combinedData$AnomPrecipitation) )
cat('The BART model RMSE is', RMSEB, 
    'This compares to the range of the outcome:', 
    min(combinedData$AnomPrecipitation), '–', 
    max(combinedData$AnomPrecipitation))

# ESTIMATE EFFECT OF WARMING ON EXTREME/ANOMALOUS PRECIPITATION FOR EACH ONI

# create list of ONI values (minimum observed, maximum observed, and zero)
ONIvals = c(min(ONIdata$JFM, ONIdata$AMJ, ONIdata$JAS, ONIdata$OND),
            max(ONIdata$JFM, ONIdata$AMJ, ONIdata$JAS, ONIdata$OND), 
            0)
# for each ONI value...
for (ONIval in ONIvals){
  # predict number of days of extreme/anomalous precipitation for 0°C warming
  zeroCelsiusPreds = predict(BARTfit, 
                             newdata = matrix(c(0, ONIval), 
                             nrow=1, 
                             byrow = TRUE))
  # predict number of days of extreme/anomalous precipitation for 1.4°C warming
  oneCelsiusPreds = predict(BARTfit, 
                            newdata = matrix(c(max(combinedData$Tanom), ONIval),
                            nrow=1, 
                            byrow = TRUE))
  # calculate CIs from predictions for 0°C and 1.4°C warming
  zeroCelsiusCI = ci(zeroCelsiusPreds, method = 'HDI')
  oneCelsiusCI = ci(oneCelsiusPreds, method = 'HDI')
  cat('For fixed ONI =', ONIval,
      ', the predicted number of days of extreme/anomalous precipitation for',
      '0°C warming is', 
      mean(zeroCelsiusPreds), '95% HDI: (', zeroCelsiusCI[[2]], '—', 
      zeroCelsiusCI[[3]], 
      '). The predicted number of days of extreme/anomalous precipitation for',
      '1.4°C warming is', 
      mean(oneCelsiusPreds), '95% HDI: (', oneCelsiusCI[[2]], '—', 
      oneCelsiusCI[[3]], 
      '). Thus, the model estimates an increase in the number of days of',
      'extreme/anomalous precipitation by', 
      (((mean(oneCelsiusPreds)/mean(zeroCelsiusPreds))-1)*100)/
        max(combinedData$Tanom),
      '% per Celsius of warming.')
}

# FIGURE 4 (BART MODEL PREDICTIONS):

# sort 0°C predictions, ascending
zeroCelsiusPredsSorted = sort(zeroCelsiusPreds, decreasing = FALSE)
# plot histogram of predictions for 0°C of warming
hist(zeroCelsiusPredsSorted, col = NULL, border = 'deepskyblue', 
     ylab = 'Density', 
     xlab = 'Predicted Days with Extreme 0-D Daily Mean Tot. Precip.', lwd = 1, 
     xlim = c(-0.5, 9), ylim = c(0, 0.55), freq = FALSE, main = '', xaxt = 'n',
     yaxt = 'n', family = 'serif', cex.lab = 1.5)
# add density line
lines(density(zeroCelsiusPredsSorted, bw = 0.3)$x, 
      density(zeroCelsiusPredsSorted, bw = 0.3)$y,
      col = 'deepskyblue', lwd = 1)
# sort 1.4°C predictions, ascending
oneCelsiusPredsSorted = sort(oneCelsiusPreds, decreasing = FALSE)
# plot histogram of predictions for 1.4°C of warming
hist(oneCelsiusPredsSorted, col = NULL, border = 'blue', lwd = 1, 
     freq = FALSE, add = TRUE, breaks = 17) # breaks chosen to match 0°C hist.
# add density line
lines(density(oneCelsiusPredsSorted, bw = 0.3)$x, 
      density(oneCelsiusPredsSorted, bw = 0.3)$y,
      col = 'blue', lwd = 1)
# add tick labels with scientific notation (x-axis only)
eaxis(1, cex.axis = 1.3, family = 'serif')
# add tick labels with scientific notation (y-axis only)
eaxis(2, cex.axis = 1.3, family = 'serif')
# set font
op = par(family = 'serif')
# plot mean predictions as vertical lines
abline(v = mean(zeroCelsiusPreds), lty = 2, lwd = 1, col = 'deepskyblue')
abline(v = mean(oneCelsiusPreds), lty = 2, lwd = 1, col = 'blue')
# plot CIs
rect(xleft = zeroCelsiusCI[[2]], xright = zeroCelsiusCI[[3]], 
     ybottom = par("usr")[3], ytop = par("usr")[4], border = NA, 
     col = adjustcolor('deepskyblue', alpha = 0.1))
rect(xleft = oneCelsiusCI[[2]], xright = oneCelsiusCI[[3]], 
     ybottom = par("usr")[3], ytop = par("usr")[4], border = NA, 
     col = adjustcolor('blue', alpha = 0.1))
# add legend to plot
legend(x = -1.0, y = 0.58, lty = c(1, 1, 2), 
       col = c('deepskyblue', 'blue', 'deepskyblue'), 
       legend=c('0°C Warming', 
                '1.4°C Warming', 
                'Averages'),
       pt.cex = 0.5, cex = 1.2, seg.len = c(1.3, 1.3, 1.3), bty = 'n', 
       x.intersp = 0.4)
legend(x = -0.92, y = 0.45, 
       col = c('lightblue1'), 
       legend=c('95% Credible \nIntervals'),
       pt.cex = 0.5, cex = 1.2, fill = c('lightblue1'), 
       border = c('lightblue1'), seg.len = c(0), bty = 'n', x.intersp = 0.68)

# delete variables that are no longer needed
remove(ONIdata)
remove(warmingData)
remove(baselineTemp)
remove(BARTfit)
remove(oneCelsiusPreds)
remove(zeroCelsiusPreds)
remove(RMSEB)
remove(ONIval)
remove(ONIvals)



#### ALTERNATIVE NEGATIVE BINOMIAL WARMING ASSOCIATION MODEL ####

# a negative binomial model is fitted as an alternative to BART for 
# investigating a possible association between extreme weather and warming. 
# This is to determine whether the warming findings are model-dependent.

# FIT ALTERNATIVE NEGATIVE BINOMIAL WARMING ASSOCIATION MODEL

# fit negative binomial model to data
altToBARTfit = glm.nb(AnomPrecipitation ~ Tanom*ONI, link = log, 
                      data = combinedData)
# calculate RMSE for negative binomial model
RMSEalt = sqrt( sum((altToBARTfit$fitted.values - 
                     combinedData$AnomPrecipitation)^2)/
                  length(combinedData$AnomPrecipitation) )
cat('The alternative negative binomial model for warming association RMSE is', 
    RMSEalt)
# print alternative negative binomial model summary
summary(altToBARTfit)
# estimate change in number of days of extreme/anomalous precipitation as 
# temperature increased
cat('The annual number of days of extreme/anomalous precipitation increased by', 
    (exp(altToBARTfit$coefficients[[2]]) - 1)*100, 
    '% (95% CI:',  (exp(confint(altToBARTfit, level = 0.95)[2, 1]) - 1)*100, 
    '–', (exp(confint(altToBARTfit, level = 0.95)[2, 2]) - 1)*100, '%)',
    'per Celcius of warming.')

# ALTERNATIVE NEGATIVE BINOMIAL WARMING ASSOCIATION MODEL DIAGNOSTICS

# plot diagnostics side-by-side
par(mfrow=c(1,2))
# residuals plot:
# get studentized residuals
r_stuAlt = rstudent(altToBARTfit)
# Shapiro-Wilk test for normality using studentized residuals (Null hypothesis: 
# residuals are normally distributed):
shapiro.test(r_stuAlt)
# plot residuals vs predicted values
plot(altToBARTfit$fitted.values, r_stuAlt, type = 'p', 
     xlab = 'Predicted Values', ylab = 'Studentized Residuals', pch = 20, 
     col = 'black', ylim=c(-5, 5))
abline(h=3, lty = 'dashed', col = 'gray') # add line indicating large residual
abline(h=0, lty = 'dashed', col = 'gray') # add zero line
abline(h=-3, lty = 'dashed', col = 'gray') # add line indicating large residual
# QQ plot:
qqnorm(r_stuAlt, xlab = 'Theoretical Values', ylab = 'Sample Values', main = '', 
       pch = 20)
qqline(r_stuAlt, col = 'gray', lwd = 2, lty = 2)
# reset to single plots
par(mfrow=c(1,1))

# delete variables that are no longer needed
remove(altToBARTfit)
remove(RMSEalt)
remove(r_stuAlt)
remove(op)