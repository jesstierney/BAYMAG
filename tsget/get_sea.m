function [t_sea,s_sea,months]=get_sea(lat,lon,foram)
% function [t_sea,s_sea,months]=get_sea(lat,lon,foram)
%
% Grabs monthly temperatures and salinities from the World Ocean Atlas 2013
% climatologies for a given location, and averages values within the
% optimal (min, max) temperature growth range of the given foraminiferal
% species. If fewer than 3 months match the growth range, this function
% takes the 3 months closest to the optimal temperature range. Foram growth
% is based on sediment trap foram fluxes in Zaric et al 2005, estimated by
% KDE in Malevich et al., 2019, Paleoceanography & Paleoclimatology,
% https://doi.org/10.1029/2019PA003576.
%
% ----- Input -----
% lat = latitude in decimal degrees (scalar or N x 1 vector)
% lon = longitude in decimal degrees (scalar or N x 1 vector, -180 to 180)
% foram = species of foram. string or N x 1 string array - must
% match the size of lat and lon and contain one of the following values:
%
% "ruber"
% "bulloides"
% "incompta"
% "pachy"
% "sacculifer"
%
% ----- Output -----
% t_sea = average temperature for months in optimal growth range
% s_sea = average salinity for months in optimal growth range
% months = cell array of months that were averaged
%
% ----- Dependencies -----
% 1) woa13_monthly.mat - .mat file containing monthly WOA13 values of SST
% and SSS.
% 2) foramseason.csv - .csv file defining min/max values for each species
% 3) read_foramseason.m - function to import foramseason.csv
% 4) EarthChordDistances_2.m
%
% function created by Dr. Jessica Tierney, The University of Arizona (2019)

%% get Nobs
Nobs=length(lat);
%put together entered locations
coords=[lon lat];
%read in optimal t ranges from csv file.
read_foramseason
min_month = 3;  %Min number of months before we look outside optimal temp range.

load('woa13_monthly.mat','lon_f','lat_f','woa13_t_mon','woa13_s_mon');

%vectorize locations
[A,B] = meshgrid(lon_f,lat_f);
c=cat(2,A',B');
locs=reshape(c,[],2);

%reorganize into vector format
n_lon=length(lon_f);
n_lat=length(lat_f);
woa13_vec_t=reshape(woa13_t_mon,n_lon*n_lat,12);
woa13_vec_s=reshape(woa13_s_mon,n_lon*n_lat,12);
%pull out locations, data where there are obs ssts.
locs_obs=locs(~isnan(woa13_vec_t(:,1)),:);
woa_obs_t=woa13_vec_t(~isnan(woa13_vec_t(:,1)),:);
woa_obs_s=woa13_vec_s(~isnan(woa13_vec_s(:,1)),:);

%grab min and max for the given species
optimin=NaN(Nobs,1);
optimax=NaN(Nobs,1);
for i=1:Nobs
    ind=strcmp(foramseason.foram,foram(i));
    %error message if I don't recognize the foram
    if sum(ind)==0
    error('foram type not recognized') 
    else
    end
    optimin(i) = foramseason.min(ind);
    optimax(i) = foramseason.max(ind);
end

%% cycle through each lat, lon, depth and extract values
%preallocate vector for seasT and dists
t_sea=NaN(Nobs,1);
s_sea=NaN(Nobs,1);
grid_lon=NaN(Nobs,1);
grid_lat=NaN(Nobs,1);
indm=NaN(Nobs,12);
month_id={'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug',...
'Sep','Oct','Nov','Dec'};
months=cell(Nobs,1);
%now loop across coords to get data.
for i=1:Nobs
    tloc=coords(i,:);
    %find closest
    [~,imin]=min(EarthChordDistances_2(locs_obs,tloc));
    woa_now_t=woa_obs_t(imin,:);
    woa_now_s=woa_obs_s(imin,:);
    %get months in range
    indt=woa_now_t >= optimin(i) & woa_now_t <= optimax(i);
    indm(i,:)=double(indt);
    gots_t=woa_now_t(indt);
    gots_s=woa_now_s(indt);
    %get months outside of range
    nots_t=woa_now_t(~indt);
    nots_s=woa_now_s(~indt);
    while length(gots_t) < min_month
       diffs(1,:)=abs(nots_t-optimin(i));
       diffs(2,:)=abs(nots_t-optimax(i));
       %find the closest value
       closest=min(diffs,[],[1 2]);
       %get its location
       [~,y]=find(diffs==closest);
       %update gots
       gots_t = [gots_t nots_t(y)];
       gots_s = [gots_s nots_s(y)];
       %update indm
       mt=ismember(woa_now_t,nots_t(y));
       %update nots
       nots_t(y)=[];
       nots_s(y)=[];
       indm(i,mt)=1;
       clear diffs
    end
    %average months for seasonal T
    t_sea(i)=mean(gots_t);
    s_sea(i)=mean(gots_s);
    %grid locs, you may output these if you like.
    grid_lon(i)=locs_obs(imin,1);
    grid_lat(i)=locs_obs(imin,2);
    %use indm to get months
    months{i}=month_id(logical(indm(i,:)));
end