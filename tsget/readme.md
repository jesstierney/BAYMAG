# TSget

Description:
This folder contains some utility functiosn to grab modern estimates of sea-surface temperature and sea-surface salinity -- as well as Late Quaternary estimates of salinity and pH -- for use with the BAYMAG Mg/Ca models.

Use: 
**get_ann.m** will simply get mean annual values for any latitude, longitude location.

**get_sea.m** will get average values for months that fall within a given foraminiferal species' "optimal" temperature range. These ranges come from kernel density estimates based on the global sediment trap data compiled by Zaric et al., 2005 (https://doi.org/10.1016/j.marmicro.2005.01.002). For details, see Malevich et al., 2019, Paleoceanography & Paleoclimatology https://doi.org/10.1029/2019PA003576.

**get_phsal.m** will return estimates of pH and salinity back in time. Inputs are the age vector of a time series, and modern values for pH and salinity for the location of the time series. The function then scales these modern values back in time using the antarctic ice core CO2 record (for pH) and a reconstruction of sea level (for salinity). The scaling assumes Last Glacial Maximum anomalies of 0.13 units for pH, and 1.1 units for salinity. The outputs can then be used to invert Mg/Ca data for temperature using **baymag_predict.m.**

IMPORTANT NOTE: The sea level reconstruction only extends back to 240,000 years, so get_phsal.m can only be used on records shorter than this.

Figure 8 in the BAYMAG manuscript shows some applications using this method.

Underlying data: 

**get_ann.m** and **get_sea.m** grab surface temperature and salinity data from the World Ocean Atlas 2013 gridded product (Boyer et al., 2013; https://www.nodc.noaa.gov/OC5/indprod.html).

**get_phsal.m** uses the antarctice ice core record of Bereiter et al., (2015) (https://doi.org/10.1002/2014GL061957) and the sea level reconstruction of Spratt and Lisiecki (2016) (https://doi.org/10.5194/cp-12-1079-2016).

Function created by Dr. Jessica Tierney, The University of Arizona (2019)



