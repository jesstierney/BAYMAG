# OMGpH

Description:
A function to grab modern estimates of surface pH and deep calcite Omega, for use with the BAYMAG Mg/Ca models.

Use: Give OMGpH the latitude, longitude, and water depth (in meters, positive value) of your core site. It will return surface water (0m) pH and omega at your water depth. You can then use these values in BAYMAG.

Underlying data: The function grabs pH and calcite Omega from GLODAPv2 (Lauvset et al., 2016, Earth System Data, https://doi.org/10.5194/essd-8-325-2016). However, GLODAPv2 lacks coverage in the Gulf of Mexico, so for sites in the GoM pH and Omega is estimated from observations taken on the GOMECC-2 cruise (https://www.aoml.noaa.gov/ocd/gcc/GOMECC2/). CO2SYS v1.1 (Van Heuven et al., 2011, https://cdiac.ess-dive.lbl.gov/ftp/co2sys/CO2SYS_calc_MATLAB_v1.1/) was used to compute pH and calcite Omega from measurements of alkalinity, DIC, salinity, temperature, pressure, silicate, and phosphate.  Resulting values of pH and Omega were averaged by depth, and smoothed with a smoothing spline to create depth profiles. These are stored in gom.mat.

Notes: GLODAP is a gridded product that interpolates values to locations that may not have any data. The Gulf of Mexico values are based on data from one cruise. So keep in mind that these are truly "estimates".

Function created by Dr. Jessica Tierney, The University of Arizona (2019)



