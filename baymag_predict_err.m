function output=baymag_predict_err(age,mg,omega_m,omega_s,salinity_m,salinity_s,pH_m,pH_s,clean,species,pstd,varargin)
% function output=baymag_predict_err(age,mg,omega,salinity,pH,clean,species,pstd,varargin)
% 
% BAYMAG inverse model to predict SSTs from given MgCa, salinity, omega, pH,
% cleaning method, and species group.
% This version, baymag_predict_err, allows the user to input prior means and standard
% deviations for omega, salinity, and pH.
% This code uses STAN: http://mc-stan.org/ for MCMC inference
% STAN is accessed through the Matlab wrapper MatlabStan, available here:
% https://github.com/brian-lau/MatlabStan
% You need to install STAN and MatlabStan on your system before using this code.
% 
% ----- Inputs -----
% age = scalar or N x 1 vector of ages. Needs to be Ma for seawater
% correction! otherwise units don't matter. 
% mg = scalar or N x 1 vector of MgCa values
% omega_m = prior mean for bottom water saturation state (scalar or N x 1
% vector)
% omega_s = prior standard deviation for bottom water saturation state (scalar or N x 1
% vector)
% salinity_m = prior mean for salinity (scalar or N x 1 vector)
% salinity_s = prior standard deviation for salinity (scalar or N x 1 vector)
% pH_m = prior mean for pH (scalar or N x 1 vector). If you are using a
% species not sensitive to pH you can enter a dummy value.
% pH_s = prior standard deviation for pH (scalar or N x 1 vector). If you are using a
% species not sensitive to pH you can enter a dummy value.
% clean = scalar or N x 1 vector to describe cleaning technique:
%  1 = reductive 
%  0 = oxidative
%   values between 0 and 1 are allowed and will be treated as a mix of
%   cleaning methods.
% species = string of target species. seven options:
% 'ruber' = Globigerinoides ruber, white or pink
% 'bulloides' = Globigerina bulloides
% 'sacculifer' = Globigerinoides sacculifer
% 'pachy' = Neogloboquadrina pachyderma
% 'incompta' = Neogloboquadrina incompta
% 'all' = pooled calibration, annual SST
% 'all_sea' = pooled calibration, seasonal SST
% pstd = prior standard deviation in degree C. Suggested values: 5-10C.
% varargin = optional arguments, right now there is just one:
%  1: a scalar to choose whether to account for changes in mgca of seawater.
%  For this to work properly your ages need to be in units of *millions of years*
%  (Ma)! If this is not entered then no seawater correction is applied.
%    0 = do not include a sw term
%    1 = include a sw term - original values from the 2019 paper
%    2 = include a sw term - v2 values incl Na/Ca (Rosenthal et al 2022) 
%  2: a value for H, the non-linear power component of the Mg/Casw to
%  Mg/Caforam relationship, to use in the seawater correction. Default is 1
%  (no non-linearity).
%  3: input different bayesian parameters. Should be a string array of size
%  3 x 1, with each column containing the parameters filename.
%
% ----- Notes -----
% 1) prior mean is set automatically using the Anand03 relationship.
% 2) No NaNs allowed!
% 3) incompta and pachy data are calibrated with the pooled pachy/incompta
%    regression.
% 4) Known bug - the stan code wants at least two values for Mg/Ca.
%    This has to do with how stan treats scalars vs vectors. So if you only
%    want to calibrate a single value of Mg/Ca, please input a second dummy
%    value.
%
% ----- Output -----
% output: a structure containing the following:
% rhat = N x 1 convergence statistic. Should be close to 1.
% neff_ratio = N x 1 ratio of effective sample size to total sample size.
% Should be above 0.1.
% SST = 2.5%, 50%, and 97.5% levels of predictions
% ens = N x 2000 ensemble of predicted SSTs
%
% ----- Dependencies -----
% 1) pooled_model_params.mat, pooled_sea_model_params.mat,
% species_model_params.mat: .mat files with default Bayesian parameters
% 2) mgsw_iters.mat: .mat file with Gaussian smoothed estimates of Mg/Ca of
% seawater through time.
% 3) mgpred_sw.stan and mgpred.stan: stan code to run the models.
% 4) MatlabStan package (see note above)
% 5) mcstan (see note above)
% 6) ChainConvergence.m: code to calculate Rhat and Neff
%
% function created by Dr. Jessica Tierney, The University of Arizona (2019)

%% deal with optional arguments
ng=nargin;
if ng==14
    sw=varargin{1};
    H=varargin{2};
    bayes=varargin{3};      
elseif ng==13
    sw=varargin{1};
    H=varargin{2};
    bayes=["pooled_model_params.mat";"pooled_sea_model_params.mat";"species_model_params.mat"];
elseif ng==12
    sw=varargin{1};
    H=1;
    bayes=["pooled_model_params.mat";"pooled_sea_model_params.mat";"species_model_params.mat"];
elseif ng==11
    sw=0;
    H=1;
    bayes=["pooled_model_params.mat";"pooled_sea_model_params.mat";"species_model_params.mat"];
else
    error('You entered too many or too few arguments');
end
%% locate where the model (.stan) files are
path_ind=which('baymag_predict');
path_ind=path_ind(1:end-16);
%% prepare data and parameters
%ensure everything is column vectors.
age=age(:);
mg=mg(:);
salinity_m=salinity_m(:);
salinity_s=salinity_s(:);
omega_m=omega_m(:);
omega_s=omega_s(:);
pH_m=pH_m(:);
pH_s=pH_s(:);
clean=clean(:);
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
omega=(omega_m.^-2).*ones(Nobs,1);
omega_s=mean([((omega_m-omega_s).^-2)-omega_m omega_m - ((omega_m+omega_s).^-2)],2).*ones(Nobs,1);
clean=clean.*ones(Nobs,1);
salinity=salinity_m.*ones(Nobs,1);
salinity_s=salinity_s.*ones(Nobs,1);
pH=pH_m.*ones(Nobs,1);
pH_s=pH_s.*ones(Nobs,1);
%check for NaNs
dats=[age mg omega clean salinity pH];
if sum(isnan(dats(:)))>0
    error('You got NaNs! Get rid of yer NaNs')
end
%load appropriate model
if  strcmp(species,'all')
    params = load(bayes(1));
    id = 1;
elseif strcmp(species,'all_sea')
    params = load(bayes(2));
    id = 1;
else
    params = load(bayes(3));
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
    %ratio to modern value, convert to log units, include H if given
    mgsw=log((mgsw./repmat(mg_mod,Nobs,1)).^H);
elseif sw==2
    %load MgCa_sw parameters
    load('mgsw_iters_v2.mat','xt','mg_smooth');
    %hold onto modern value
    mg_mod=mg_smooth(1,:);
    %interpolate to given ages
    mgsw=interp1(xt,mg_smooth,age);
    %ratio to modern value, convert to log units, include H if given
    mgsw=log((mgsw./repmat(mg_mod,Nobs,1)).^H);
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

%put data into struct
if sw==1
mg_dat = struct('N',Nobs,'M',Mparams,'omega_m',omega,'omega_s',omega_s,'clean',clean,'s_m',salinity,'s_s',salinity_s,'ph_m',pH,... %,
    'ph_s',pH_s,'mg',log(mg),'betaT',betaT,'betaO',betaO,'betaC',betaC,'betaS',betaS,'betaP',betaP,...
    'sigma',sigma,'alpha',alpha,'prior_mu',prior_mu,'prior_sig',prior_sig,'id',id,'mgsw',mgsw);
else
mg_dat = struct('N',Nobs,'M',Mparams,'omega_m',omega,'omega_s',omega_s,'clean',clean,'s_m',salinity,'s_s',salinity_s,'ph_m',pH,... %,
    'ph_s',pH_s,'mg',log(mg),'betaT',betaT,'betaO',betaO,'betaC',betaC,'betaS',betaS,'betaP',betaP,...
    'sigma',sigma,'alpha',alpha,'prior_mu',prior_mu,'prior_sig',prior_sig,'id',id);
end
%initialize to improve convergence
mg_init = struct('t',repmat(prior_mu,1,Mparams),'s',repmat(salinity,1,Mparams),'ph',repmat(pH,1,Mparams),'omega',repmat(omega,1,Mparams));

%stan settings. Designed for fast convergence, can be changed if necessary.
chains=4;
warmup=150;
iters=100;
thin=4;
%% run Stan, time the run, use verbose option to see progress.
% Note: the first time you run this, it will take longer because Stan will
% have to compile the model. Once the models are compiled subsequent runs
% are faster.
if sw==1
    file_name=strcat(path_ind,'mgpred_werr_sw.stan');
else
    file_name=strcat(path_ind,'mgpred_werr.stan');
end

tic
    fit = stan('file',file_name,'data',mg_dat,'chains',chains,'warmup',warmup,'iter',iters,'thin',thin,'init',mg_init,'verbose',true);
fit.block();
toc
%% extract and analyze data
samples_t=NaN(chains,iters/thin,Nobs,Mparams);
samples_s=NaN(chains,iters/thin,Nobs,Mparams);
samples_o=NaN(chains,iters/thin,Nobs,Mparams);
samples_p=NaN(chains,iters/thin,Nobs,Mparams);
for i=1:chains
    samples_t(i,:,:,:)=fit.sim.samples(i).t;
    samples_s(i,:,:,:)=fit.sim.samples(i).s;
    samples_o(i,:,:,:)=fit.sim.samples(i).omega;
    samples_p(i,:,:,:)=fit.sim.samples(i).ph;
end
%reshape for Rhat, Neff calc
samples_t=permute(samples_t,[1 3 2 4]);
samples_s=permute(samples_s,[1 3 2 4]);
samples_o=permute(samples_o,[1 3 2 4]);
samples_p=permute(samples_p,[1 3 2 4]);
samples_r_t=reshape(samples_t,chains,Nobs,iters/thin*Mparams);
samples_r_s=reshape(samples_s,chains,Nobs,iters/thin*Mparams);
samples_r_o=reshape(samples_o,chains,Nobs,iters/thin*Mparams);
samples_r_p=reshape(samples_p,chains,Nobs,iters/thin*Mparams);

output.rhat_t=NaN(Nobs,1);
output.rhat_s=NaN(Nobs,1);
output.rhat_o=NaN(Nobs,1);
output.rhat_p=NaN(Nobs,1);
neff_t=NaN(Nobs,1);
neff_s=NaN(Nobs,1);
neff_o=NaN(Nobs,1);
neff_p=NaN(Nobs,1);
%note call to ChainConvergence
for i=1:Nobs
    [output.rhat_t(i),neff_t(i)]=ChainConvergence(squeeze(samples_r_t(:,i,:)),chains);
    [output.rhat_s(i),neff_s(i)]=ChainConvergence(squeeze(samples_r_s(:,i,:)),chains);
    [output.rhat_o(i),neff_o(i)]=ChainConvergence(squeeze(samples_r_o(:,i,:)),chains);
    [output.rhat_p(i),neff_p(i)]=ChainConvergence(squeeze(samples_r_p(:,i,:)),chains);
end
output.neff_ratio_t=neff_t./(chains*iters/thin*Mparams);
output.neff_ratio_s=neff_s./(chains*iters/thin*Mparams);
output.neff_ratio_o=neff_o./(chains*iters/thin*Mparams);
output.neff_ratio_p=neff_p./(chains*iters/thin*Mparams);

t_out=permute(samples_r_t,[2 1 3]);
t_out=reshape(t_out,Nobs,chains*iters/thin*Mparams);

s_out=permute(samples_r_s,[2 1 3]);
s_out=reshape(s_out,Nobs,chains*iters/thin*Mparams);

omega_out=permute(samples_r_o,[2 1 3]);
omega_out=reshape(omega_out,Nobs,chains*iters/thin*Mparams);

ph_out=permute(samples_r_p,[2 1 3]);
ph_out=reshape(ph_out,Nobs,chains*iters/thin*Mparams);

%save 2000 iterations in output
output.ens.t=t_out(:,1:10:end);
output.ens.s=s_out(:,1:10:end);
output.ens.omega=omega_out(:,1:10:end);
output.ens.ph=ph_out(:,1:10:end);

pers3=[.025 .5 .975].*size(t_out,2);
sst_s=sort(t_out,2);
sss_s=sort(s_out,2);
om_s=sort(omega_out,2);
ph_s=sort(ph_out,2);

output.SST=sst_s(:,pers3);
output.SSS=sss_s(:,pers3);
output.Omega=om_s(:,pers3);
output.Ph=ph_s(:,pers3);
%also print out the species used
output.species = species2;
%% some sanity check plots
%plot prior and posterior.
f1=figure(1); clf;
set(f1,'pos',[50 700 400 400]);
xt=(-2:.1:40)';
prior=normpdf(xt,mean(prior_mu),pstd);
post=ksdensity(output.ens.t(:),xt);
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