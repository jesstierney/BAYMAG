function [value,label,dists] = nearestgriddedvalue_mg(locs_lat,locs_lon,species)
% this function finds the average mean annual or seasonal SST for a given
% foraminiferal species group at a given location. Used by baymag_predict.m.
% INPUT: 
%  locs_lat = N x 1 vector of latitudes of sites in decimal degrees
%  locs_lon = N x 1 vector longitudes of sites in decimal degrees
%  species: seven options:
% 'ruber' = Globigerinoides ruber, white or pink
% 'bulloides' = Globigerina bulloides
% 'sacculifer' = Globigerinoides sacculifer
% 'pachy' = Neogloboquadrina pachyderma
% 'incompta' = Neogloboquadrina incompta
% 'all' = pooled calibration, annual SST
% 'all_sea' = pooled calibration, seasonal SST
%OUTPUT:
% value = nearest average SST value
% label = seasonality (list of months or 'Annual')
% dist = distance between site and nearest gridded value in km
% 
% function created by Dr. Jessica Tierney, The University of Arizona (2019)

%% put locs data together
coords = [locs_lon locs_lat];
%load gridded data
lat=double(ncread('foram_seasons.nc','lat'));
lon=double(ncread('foram_seasons.nc','lon'));
%load appropriate species
if  strfind(species,char('ruber'))==1
    field=double(ncread('foram_seasons.nc','ruber_t'));
    field_labels=ncread('foram_seasons.nc','ruber_months');
elseif strfind(species,char('bulloides'))==1
    field=double(ncread('foram_seasons.nc','bulloides_t'));
    field_labels=ncread('foram_seasons.nc','bulloides_months');
elseif strfind(species,char('sacculifer'))==1
    field=double(ncread('foram_seasons.nc','sacculifer_t'));
    field_labels=ncread('foram_seasons.nc','sacculifer_months');
elseif strfind(species,char('pachy'))==1
    field=double(ncread('foram_seasons.nc','pachy_t'));
    field_labels=ncread('foram_seasons.nc','pachy_months');
    elseif strfind(species,char('incompta'))==1
    field=double(ncread('foram_seasons.nc','incompta_t'));
    field_labels=ncread('foram_seasons.nc','incompta_months');
else
    %field for all-species model. labels are a dummy value.
    field=double(ncread('foram_seasons.nc','ann_t'));
    field_labels=ncread('foram_seasons.nc','ruber_months');
end

%reshape grid into vector form:
field_vec=reshape(field,size(field,1)*size(field,2),1);
%reshape labels into vector form:
field_vec_labels=reshape(field_labels,size(field_labels,1)*size(field_labels,2),size(field_labels,3));

%vectorize grid locations
[A,B] = meshgrid(lon,lat);
c=cat(2,A',B');
locs_vec=reshape(c,[],2);

%take subset without NaNs
inder_g=~isnan(field_vec);
locs_obs=locs_vec(inder_g,:);
vec_obs=field_vec(inder_g);
labels_obs=field_vec_labels(inder_g,:);

%set a km cutoff in case something bad happens:
max_dist=500;

%note call to EarthChordDistances_2
dists_all=EarthChordDistances_2(locs_obs,coords);
%get the min distance and index
[dmin,imin]=min(dists_all);
dists=dmin;  
value=vec_obs(imin);
label_t=logical(labels_obs(imin,:));
months={'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
label=months(label_t);
if size(label,2)==12
    %label as annual if all months are there
   label='Annual';
elseif strfind(species,char('all'))==1
    %label as annual if all species model is used
   label='Annual';
elseif strfind(species,char('all_sea'))==1
    %label as annual if all_sea species model is used
    label='Annual';
elseif dists>max_dist
    %print a special warning if the site is far from the ocean
   label='WARNING: site is more than 500 km from an ocean grid cell. Check your lat/lon.';
else
end