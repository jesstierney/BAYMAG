function [salinity,ph]=get_phsal(age,sal,ph_in)
%simple function to get idealized changes in salinity and pH in the Late
%Quaternary, for use with the BAYMAG prediction model. Modified from the
%technique in MgCaRB v1 by David Evans (2019).
%
%INPUTS:
%age = age in yr BP - must be younger than 240 ka.
%sal = modern observed salinity
%ph = modern observed ph
load('antco2.mat','age_gas_calBP','co2_ppm');
load('sealevel.mat','seaLevel');

if max(age) > max(seaLevel(:,1).*1000)
    error('Your time series is too long to use get_phsal.m. Sorry!')
else
end
   
%assumptions about 2-sigma bounds of change in salinity and ph:
lgm_sal = 0.6; %psu
lgm_ph = 0.1; %ph total
%years BP for preindustrial (1750 AD)
pi1=150;
ind_pi_ph=findnearest(pi1,age_gas_calBP);
ind_pi_sal=findnearest(pi1,seaLevel(:,1).*1000);

sal_scaled = normalize(seaLevel(:,2)).*lgm_sal/2;
sal_scaled = sal_scaled - sal_scaled(ind_pi_sal);
salinity = sal - interp1(seaLevel(:,1).*1000,sal_scaled,age,'linear','extrap');
ph_scaled = normalize(co2_ppm).*lgm_ph/2;
ph_scaled = ph_scaled - ph_scaled(ind_pi_ph);
ph = ph_in - interp1(age_gas_calBP,ph_scaled,age,'linear','extrap');
