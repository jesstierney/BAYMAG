function [t_ann,s_ann]=get_ann(lat,lon)
% function [t_ann,s_ann]=get_ann(lat,lon)
%
% Simple function to grab annual average temperature and salinity from the
% World Ocean Atlas 2013 for any given location
%
% ----- Input -----
% lat = latitude in decimal degrees (scalar or N x 1 vector)
% lon = longitude in decimal degrees (scalar or N x 1 vector, -180 to 180)
%
% ----- Output -----
% t_ann = mean annual temperature
% s_ann= mean annual salinity
%
% ----- Dependencies -----
% 1) woa13_annual.mat - .mat file containing WOA13 mean annual SST and SSS.
% 2) EarthChordDistances_2.m
%
% function created by Dr. Jessica Tierney, The University of Arizona (2019)

%% get Nobs
Nobs=length(lat);
%put together entered locations
coords=[lon lat];

load('woa13_annual.mat','lon_f','lat_f','woa13_t_ann','woa13_s_ann');

%vectorize locations
[A,B] = meshgrid(lon_f,lat_f);
c=cat(2,A',B');
locs=reshape(c,[],2);

%reorganize into vector format
n_lon=length(lon_f);
n_lat=length(lat_f);
woa13_vec_t=reshape(woa13_t_ann,n_lon*n_lat,1);
woa13_vec_s=reshape(woa13_s_ann,n_lon*n_lat,1);
%pull out locations, data where there are obs ssts.
locs_obs=locs(~isnan(woa13_vec_t),:);
woa_obs_t=woa13_vec_t(~isnan(woa13_vec_t),:);
woa_obs_s=woa13_vec_s(~isnan(woa13_vec_s),:);

%% cycle through each lat, lon, depth and extract values
%preallocate vector for seasT and dists
t_ann=NaN(Nobs,1);
s_ann=NaN(Nobs,1);
%now loop across coords to get data.
for i=1:Nobs
    tloc=coords(i,:);
    %find closest
    [~,imin]=min(EarthChordDistances_2(locs_obs,tloc));
    t_ann(i)=woa_obs_t(imin);
    s_ann(i)=woa_obs_s(imin);
end