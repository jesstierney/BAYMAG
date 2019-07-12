TSget

Description:
A function to grab modern estimates of sea-surface temperature and sea-surface salinity, for use with the BAYMAG Mg/Ca models.

Use: get_ann.m will simply get mean annual values for any latitude, longitude location. get_sea.m will get average values for months that fall within a given foraminiferal species' "optimal" temperature range. These ranges come from kernel density estimates based on the global sediment trap data compiled by Zaric et al., 2005 (https://doi.org/10.1016/j.marmicro.2005.01.002). For details, see Malevich et al., 2019, Paleoceanography & Paleoclimatology https://doi.org/10.1029/2019PA003576.

Underlying data: The functions grab surface temperature and salinity data from the World Ocean Atlas 2013 gridded product (Boyer et al., 2013; https://www.nodc.noaa.gov/OC5/indprod.html).

Function created by Dr. Jessica Tierney, The University of Arizona (2019)



