data {
    int<lower=0> N; // number of data items
    int<lower=0> M; // number of posterior samples
    int<lower=0> id; // id flag for species
    vector[N] mg; // mg predictor
    vector[N] omega; // omega predictor
    vector[N] clean; // clean predictor
    vector[N] s; // salinity predictor
    vector[N] ph; // ph predictor
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
}
model {
   vector[N] mu; //mean value
for (m in 1:M) {
    t[:,m] ~ normal(prior_mu,prior_sig);
   if (id < 3) {
      mu = alpha[m] + t[:,m] * betaT[m] + s * betaS[m] + ph * betaP[m] + omega * betaO[m] + (1 - clean * betaC[m]);
   } else {
      mu = alpha[m] + t[:,m] * betaT[m] + s * betaS[m] + omega * betaO[m] + (1 - clean * betaC[m]);
    }
    mg ~ normal(mu,sigma[m]);
 }
}