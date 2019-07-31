# BAYMAG

The BAYMAG package is a set of Bayesian hierarchical models for Mg/Ca in planktic foraminifera. The functions here can be used to both predict Mg/Ca from T, S, pH, Omega, and cleaning method, and predict T given constraints on the other four predictor variables. For a full discussion of this model, please see Tierney et al., 2019, Paleoceanography & Paleoclimatology (in prep).

A quick guide to basic use:

To model Mg/Ca from T, S, pH, Omega, and cleaning values:
Use baymag_forward.m. This function calculates Mg/Ca values using posterior draws from the calibration model (stored in the params.mat files). It has the option to add a normal prior if you would like to place some restrictions on posterior Mg/Ca, and to account for changing Mg/Ca of seawater (for longer geological timescales). baymag_forward has no dependent functions but does need the params.mat files.

To model SST from Mg/Ca, S, pH, Omega, and cleaning values:
Use baymag_predict.m. This function calculates SST given Mg/Ca and inputs of S, pH, Omega, and cleaning through Bayesian inference. It has an option to correct for Mg/Ca of seawater changes. The Bayesian models are written in Stan (mgpred.stan and mgpred_sw.stan). You will need to install both Stan and MatlabStan in order to use this function. Download the latest version of the command line version of Stan here: https://github.com/stan-dev/cmdstan/releases. BAYMAG has been tested up to v2.19.1. MatlabStan is here: https://github.com/brian-lau/MatlabStan/releases. BAYMAG has been tested with v2.15.1.0.  Alternatively, you can run the Stan models from Python or R, but in that case you will need to write your own wrapper like baymag_predict.

To model SST from Mg/Ca, S, pH, Omega, and cleaning values and include errors in S, pH, and Omega:
Use baymag_predict_err.m. This version allows the user to input priors for T, S, pH, and Omega and therefore include uncertainties in S, pH, and Omega. It iteratively solves for posteriors for all four variables, thus be advised that it is quite a bit slower than baymag_predict.m. It also has an option for a seawater correction.

In addition to Stan, baymag_predict requires the params.mat files and the dependent functions ChainConvergence.m and EarthChordDistances_2.m. ChainConvergence is used to calculate the Rhat and Neff statistics (c.f. Gelman et al., "Bayesian Data Analysis") to assess convergence.

The omegaph and tsget folders contain the accessory functions omgph.m, get_sea.m, and get_ann.m. These can be used to grab modern values of omega, pH, temperature, and salinity for any location. They require EarthChordDistances_2.m to calculate the nearest gridded location based on chordal distance, and the netcdfs and .mat files in the folders which contain GLODAPv2 and WOA13 data.



