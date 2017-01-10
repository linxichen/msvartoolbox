%% housekeeping
clear;
clc;
close all
addpath(genpath('../'))

load sim_for_mmu

%% specify the empirical model
model.p = 3;
model.K = model.N*model.p+1;
model.T = length(regimes_mat);
model.breakdate = ceil(0.5*model.T);
condition.regimes = regimes_mat;
condition.y_table = y;

%% specify priors for hyperpara, very flat
p = model.p;
N = size(y,2);
M1 = max(regimes_mat(:,1));
M2 = max(regimes_mat(:,2));
priors.mmu_mean = zeros(N*2*M1,1);
priors.mmu_cov = 1e5*eye(N*2*M1);
priors.pphi_mean = zeros(N*N*p*M1,1);
priors.pphi_cov = 1e5*eye(N*N*p*M1);
priors.Ssigma_SSR = 0;
priors.Ssigma_nnu = 0;

%% specify options for the Gibbs sampler
options.burnin = 2e4;
options.R = 5e4;
options.pphi_init = zeros(p*N*N*M1,1);
options.mmu_init = zeros(N*2*M1,1);
for i_regime2 = 1:M2
	options.Ssigma_array_init(:,:,i_regime2) = eye(N);
end

tic
draws = dummyVAR_Gibbs(y,regimes_mat,model,priors,options);
toc
%% look at result
median(draws.Ssigma_array,4)
mean(draws.pphi,2)
mean(draws.mmu,2)
