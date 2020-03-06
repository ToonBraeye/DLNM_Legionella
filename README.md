# DLNM_Legionella
A short simulation study that creates and analyses data on Legionnaires' disease and meteorological variables

# Research question
Effect of meteorological variables during the risk period (10-2 days prior to disease onset) on the Legionnaires' disease incidence

# Data Input
The code requires two datasets:
1.meteoProvS.Rdata: standardized/transformed meteorological variables (temperature, wind speed, relative humidity, precipitation) by province by day
2.pop_prov_year: population totals by province by year 

The legionnaires' disease incidence data is a function of the meteorological data: negative effect of wind speed on day 4-6, positive effect of relative humidity on day 5 and precipitation on day 6

# Analysis
The code includes two case-crossover designs:
The case-crossover design eliminates time-varying confounding by comparing values to values measured on corresponding days (same day of the year in other study years).

1. Single day model: A model is fitted for each day in the risk period (e.g. association between disease incidence and meteorological variables six days prior). The meteorological variables are modelled as factors (quantiles) to allow for non-linear associations.
2. DLNM: A single model fits trends from the complete risk period to the observed Legionnaires' disease incidence. The crossbases are piecewise-linear for the lag-dimension and linear for the meteorological variables (with the exception of temperature, which is included as a cubic spline).

# Results presentation
1. Single day model: a single boxplot with all the significant coefficients obtained from the single day models
2. DLNM: A graph in which the results of the DLNM are presented accumulated over the risk period and one in which results by lag (day) and variable are presented. 
