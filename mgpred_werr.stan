data {
    int<lower=0> N; // number of data items
    int<lower=0> M; // number of posterior samples
    int<lower=0> id; // id flag for species
    vector[N] mg; // mg predictor
    vector[N] omega_m; // omega prior mean
    vector[N] omega_s; // omega prior sigma
    vector[N] clean; // clean predictor
    vector[N] s_m; // salinity prior mean
    vector[N] s_s; // salinity prior sigma
    vector[N] ph_m; // ph prior mean
    vector[N] ph_s; // ph prior sigma
    vector[M] betaT; // betaT
    vector[M] betaO; // betaO
    vector[M] betaC; // betaC
    vector[M] betaS; // betaS
    vector[M] betaP; // betaSP
    vector[M] alpha; // alpha
    vector[M] sigma; // sigma
    vector[N] prior_mu;
    real prior_sig;
}
parameters {
    matrix<lower=-2.5>[N,M] t; // temperature to estimate
    matrix<lower=0>[N,M] s; // salinity to estimate
    matrix<lower=0,upper=14>[N,M] ph; // ph to estimate
    matrix<lower=0>[N,M] omega; // omega to estimate
}
model {
   vector[N] mu; //mean value
for (m in 1:M) {
    //set priors
    t[:,m] ~ normal(prior_mu,prior_sig);
    s[:,m] ~ normal(s_m,s_s);
    ph[:,m] ~ normal(ph_m,ph_s);
    omega[:,m] ~ normal(omega_m,omega_s);
   if (id < 3) {
      mu = alpha[m] + t[:,m] * betaT[m] + s[:,m] * betaS[m] + ph[:,m] * betaP[m] + omega[:,m] * betaO[m] + (1 - clean * betaC[m]);
   } else {
      mu = alpha[m] + t[:,m] * betaT[m] + s[:,m] * betaS[m] + omega[:,m] * betaO[m] + (1 - clean * betaC[m]);
    }
    mg ~ normal(mu,sigma[m]);
 }
}