function output=baymag_predict(lat,lon,age,mg,omega,salinity,pH,clean,species,pstd,varargin)
% BAYMAG inverse model to predict SSTs from given MgCa, salinity, omega, pH,
% cleaning method, and species group.
% This code uses STAN: http://mc-stan.org/ for MCMC inference
% STAN is accessed through the Matlab wrapper MatlabStan, available here:
% https://github.com/brian-lau/MatlabStan
% You need to install STAN and MatlabStan on your system before using this code.
% 
% INPUTS:
% lat = latitude of site (used to pull modern SST/SSS)
% lon = longitude of site (used to pull modern SST/SSS)
% age = scalar or N x 1 vector of ages. Needs to be Ma for seawater
% correction! otherwise units don't matter. 
% mg = scalar or N x 1 vector of MgCa values
% omega = scalar or N x 1 vector of bottom water saturation state
% salinity = scalar or N x 1 vector of salinity (psu)
% pH = scalar or N x 1 vector of pH (total scale). If you are using a
% species not sensitive to pH you can enter a dummy value.
% clean = scalar to describe cleaning technique:
%  1 = reductive 
%  0 = oxidative
%   values between 0 and 1 are allowed and will be treated as a mix of
%   cleaning methods.
% pstd = prior standard deviation in degree C. Suggested values: 5-10C.
% species = string of target species. seven options:
% 'ruber' = Globigerinoides ruber, white or pink
% 'bulloides' = Globigerina bulloides
% 'sacculifer' = Globigerinoides sacculifer
% 'pachy' = Neogloboquadrina pachyderma
% 'incompta' = Neogloboquadrina incompta
% 'all' = pooled calibration, annual SST
% 'all_sea' = pooled calibration, seasonal SST
% varargin = optional arguments, right now there is just one allowed:
%  1: a scalar to choose whether to account for changes in mgca of seawater.
%  For this to work properly your ages need to be in units of *millions of years*
%  (Ma)! If this is not entered then no seawater correction is applied.
%    0 = do not include a sw term
%    1 = include a sw term
%
% NOTES:
% 1) prior mean is set automatically using the Anand03 relationship.
% 2) No NaNs allowed!
% 3) incompta and pachy data are calibrated with the pooled pachy/incompta
%    regression.
%
% OUTPUTS in structure output:
% rhat = N x 1 convergence statistic. Should be close to 1.
% neff_ratio = N x 1 ratio of effective sample size to total sample size.
% Should be above 0.1.
% SST = 2.5%, 50%, and 97.5% levels of predictions
% ens = N x 2000 ensemble of predicted SSTs
% modernSST = modern, seasonally-averaged SST at site from World Ocean Atlas 2013
% seasonality = estimate of seasonality of derived SSTs at site, based on modern
% sediment trap flux data
%
% function created by Dr. Jessica Tierney, The University of Arizona (2019)

%% deal with optional arguments
ng=nargin;
if ng==11
    sw=varargin{1};
elseif ng==10
    sw=0;
else
    error('You entered too many or too few arguments');
end
%% locate where the model (.stan) files are
path_ind=which('baymag_predict');
path_ind=path_ind(1:end-16);
%% prepare data and parameters
species_list = {'ruber','bulloides','sacculifer','pachy','incompta','all','all_sea'};
species_list_model = {'ruber','bulloides','sacculifer','pachy'};
%check that you have an allowable species
if ~ismember(species,species_list)
    error('Species not recognized');
end
%check for incompta. Change to pachy but retain for nearestgrid function.
if ismember(species,{'incompta'})
    species2='incompta';
    species='pachy';
else
    species2=species;
end    
id = (1:1:4);

%define dimensions
Nobs=length(mg);
%transform variables and vectorize
omega=(omega.^-2).*ones(Nobs,1);
clean=clean.*ones(Nobs,1);
salinity=salinity.*ones(Nobs,1);
pH=pH.*ones(Nobs,1);
%check for NaNs
dats=[age mg omega clean salinity pH];
if sum(isnan(dats(:)))>0
    error('You got NaNs! Get rid of yer NaNs')
end
%load appropriate model
if  strfind(species,char('all'))==1
    params = load('pooled_model_params.mat');
    id = 1;
elseif strfind(species,char('all_sea'))==1
    params = load('pooled_sea_model_params.mat');
    id = 1;
else
    params = load('species_model_params.mat');
    %grab id location for species.
    id = id(ismember(species_list_model,species));
end

%subsample posterior to 200 draws
subs=5;
betaT=params.betaT(1:subs:end);
betaC=params.betaC(1:subs:end);
betaO=params.betaO(1:subs:end);
betaS=params.betaS(1:subs:end);
betaP=params.betaP(1:subs:end);
%for sigma and alpha
sigma=params.sigma(1:subs:end,id);
alpha=params.alpha(1:subs:end,id);

%define number of parameters sampled.
Mparams=length(betaT);

%get mgca sw estimates
if sw==1
    load('mgsw_iters.mat','xt','mg_smooth');
    %subsample to 200 draws
    mg_smooth=mg_smooth(:,1:subs:end);
    mg_mod=mg_smooth(1,:);
    mgsw=interp1(xt,mg_smooth,age,'linear','extrap');
    %ratio to modern value, convert to log units
    mgsw=log(mgsw./repmat(mg_mod,Nobs,1));
else
end

%assign prior mean - use Anand approximation, sw corrected if needed.
if sw==1
    mgsw_est=median(exp(mgsw),2);
    prior_mu=(log(mg)-log(.38.*mgsw_est))./.09;
else
    prior_mu=(log(mg)-log(.38))./.09;
end
%assign prior standard deviation
prior_sig=pstd;

if sw==1
%put data into struct
mg_dat = struct('N',Nobs,'M',Mparams,'omega',omega,'clean',clean,'s',salinity,'ph',pH,...
    'mg',log(mg),'betaT',betaT,'betaO',betaO,'betaC',betaC,'betaS',betaS,'betaP',betaP,...
    'sigma',sigma,'alpha',alpha,'prior_mu',prior_mu,'prior_sig',prior_sig,'id',id,'mgsw',mgsw);
else
%put data into struct
mg_dat = struct('N',Nobs,'M',Mparams,'omega',omega,'clean',clean,'s',salinity,'ph',pH,...
    'mg',log(mg),'betaT',betaT,'betaO',betaO,'betaC',betaC,'betaS',betaS,'betaP',betaP,...
    'sigma',sigma,'alpha',alpha,'prior_mu',prior_mu,'prior_sig',prior_sig,'id',id);
end

%initialize at prior mean because it improves convergence
mg_init = struct('t',repmat(prior_mu,1,Mparams));

%stan settings. The absolute minimum so that the model runs fast. Increase
%warmup and iters if you are not getting good convergence.
chains=4;
warmup=150;
iters=100;
thin=4;
%% run Stan, time the run, use verbose option to see progress.
% Note: the first time you run this, it will take longer because Stan will
% have to compile the model. Once the models are compiled subsequent runs
% are faster.
if sw==1
    file_name=strcat(path_ind,'mgpred_sw.stan');
else
    file_name=strcat(path_ind,'mgpred.stan');
end

tic
    fit = stan('file',file_name,'data',mg_dat,'chains',chains,'warmup',warmup,'iter',iters,'thin',thin,'init',mg_init,'verbose',true);
fit.block();
toc
%% extract and analyze data
samples=NaN(chains,iters/thin,Nobs,Mparams);
for i=1:chains
    samples(i,:,:,:)=fit.sim.samples(i).t;
end
%reshape for Rhat, Neff calc
samples=permute(samples,[1 3 2 4]);
samples_r=reshape(samples,chains,Nobs,iters/thin*Mparams);

output.rhat=NaN(Nobs,1);
neff=NaN(Nobs,1);
%note call to ChainConvergence
for i=1:Nobs
    [output.rhat(i),neff(i)]=ChainConvergence(squeeze(samples_r(:,i,:)),chains);
end
output.neff_ratio=neff./(chains*iters/thin*Mparams);
sst=permute(samples_r,[2 1 3]);
sst=reshape(sst,Nobs,chains*iters/thin*Mparams);
%save 2000 iterations in output
output.ens=sst(:,1:10:end);
pers3=[.025 .5 .975].*size(sst,2);
sst_s=sort(sst,2);
output.SST=sst_s(:,pers3);

%% extract modern sst and seasonality
% note call to nearestgriddedvalue_mg
[output.modernSST,output.seasonality] = nearestgriddedvalue_mg(lat,lon,species2);
% also print out the species used
output.species = species2;
%% some sanity check plots
%plot prior and posterior.
f1=figure(1); clf;
set(f1,'pos',[50 700 400 400]);
xt=(-2:.1:40)';
prior=normpdf(xt,mean(prior_mu),pstd);
post=ksdensity(output.ens(:),xt);
pr=plot(xt,prior,'k-','linewidth',1); hold on;
pt=plot(xt,post,'b-','linewidth',1);
legend([pr pt],'Prior','Posterior');

%plot inferred SSTs with 95% CI. note this does not plot on age (in case
%ages are non-sequential)
f2=figure(2); clf;
set(f2,'pos',[550 700 500 400]);
 p1=plot(output.SST(:,2),'color','k','linewidth',2);
 hold on;
 p2=plot(output.SST(:,1),'color',[.6 .6 .6],'linewidth',1);
 plot(output.SST(:,3),'color',[.6 .6 .6],'linewidth',1);
 legend([p1 p2],'Median','2.5 and 97.5');