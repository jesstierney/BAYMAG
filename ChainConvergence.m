function [Rhat,Neff]=ChainConvergence(chains,M)
% Calculates the R-hat and Neff statistics from Gelman et al., "Bayesian
% Data Analysis" for monitoring  convergence of multiple parallel MCMC
% runs.
%
% INPUTS:
% chains: a matrix of MCMC chains, each of the same length. 
% M: the number of different chains - must be one of the dimensions of
% chains. 
%
% Written by Martin P. Tingley (2012)
%
if ismember(M, size(chains))==0
    disp('Error: Second input must be the number of chains, so one of the dimensions of the first input.')
    return
elseif length(chains(:,1))==M
    chains=chains';
end

% each column is a chain. 
N=length(chains(:,1));

%quantities from Gelman:

psi_bar_dot_j = mean(chains);
psi_bar_dot_dot=mean(psi_bar_dot_j);

B= (N/(M-1))*sum((psi_bar_dot_j-psi_bar_dot_dot).^2);

s2_j= (1/(N-1))*sum((chains-kron(psi_bar_dot_j, ones(N,1))).^2);

W=(1/M)*sum(s2_j);

var_hat_pos=((N-1)/N)*W + (1/N)*B;

Rhat=sqrt(var_hat_pos/W);
Neff=M*N*min(var_hat_pos/B, 1);