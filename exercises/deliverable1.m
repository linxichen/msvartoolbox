%% housekeeping
clear;
clc;
close all
addpath(genpath('../'))

% load data
load manual_select_dummy_FRED
start = 0; % first start date
shift = 0; % shift some of the variable forward or backward
% forward_resid_rCIPI = data(start+2+shift:end,2);
% share_rCIPI_pot = (rCIPI(2:end) - rCIPI(1:end-1))./GDPPOT(1:end-1);
share_rCIPI_pot = share_rCIPI_potential(2:end);
% forward_stock = data(start+1+shift:end,3);
% match_GDP = ln_rGDP;
match_GDP = (rGDP(2:end) - rGDP(1:end-1))./rGDP(1:end-1);
match_sales = (rSalesGoods(2:end) - rSalesGoods(1:end-1))./rSalesGoods(1:end-1);
data = 100*[match_sales match_GDP share_rCIPI_pot];

% data label
datelabel = (1949.25:0.25:2016.50)'; % because first differenced
yearnum = floor(datelabel);
monthnum = 12*(datelabel - yearnum)+2;
date_serial = datenum(yearnum,monthnum,ones(size(yearnum)));
break_dummy = (datelabel >= 1984.00);

%% load regimes from MLE
load regimes_MLE regimes_mat

%% specify the empirical model
model.p = 4;
model.T = length(regimes_mat);
model.breakdate = find(datelabel == 1984.00,1,'first');
condition.regimes = regimes_mat;
condition.y_table = data;

%% specify priors for hyperpara, very flat
p = model.p;
N = size(data,2);
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
draws = dummyVAR_Gibbs(data,regimes_mat,model,priors,options);
toc
%% look at result
median(draws.Ssigma_array,4)
mean(draws.pphi,2)
mean(draws.mmu,2)
