%% housekeeping
clear;
clc;
close all
addpath(genpath('../'))

%% set hyperparameters
var1 = 0.5;
var2 = 2.2;
cov12 = -0.5;
true_Ssigma1 = 0.1*[var1, cov12; cov12, var2];
true_Ssigma2 = 0.1*[0.37, 0.11; 0.11, 0.12];
% true_Ssigma1 = eye(2);
% true_Ssigma2 = eye(2);
true_pphi1 = 0.8*eye(2);
true_pphi2 = -0.6*eye(2);
true_pphi = [true_pphi1;true_pphi2];
true_mmu1 = [0.1;0.2];
true_mmu2 = [0.3;0.4];

%% Simulate data, painstakingly sim a VAR(2) with switching mean and covariance
desired_size = 500; 
burnin = ceil(2*desired_size);
T = burnin+desired_size;
date = (1:T)';
regimes = ones(T,1);
regime2_marker = sort(datasample(date,ceil(0.3*T),'Replace',false));
regimes(regime2_marker) = 2;
dum_regime1 = regimes == 1;
dum_regime2 = regimes == 2;
Ssigma_series = kron(dum_regime1,true_Ssigma1)+kron(dum_regime2,true_Ssigma2);
mmu_series = kron(dum_regime1,true_mmu1)+kron(dum_regime2,true_mmu2);
y_zero = zeros(2,1);
y_minusone = zeros(2,1);

t = 1;
mmu_now = dum_regime1(t)*true_mmu1+dum_regime2(t)*true_mmu2;
mmu_lag1 = true_mmu1;
mmu_lag2 = true_mmu2;
ylag1 = y_zero;
ylag2 = y_minusone;
Ssigma_now = dum_regime1(t)*true_Ssigma1+dum_regime2(t)*true_Ssigma2;
y(:,t) = mmu_now + [true_pphi1 true_pphi2]*([ylag1; ylag2]-[mmu_lag1; mmu_lag2]) + mvnrnd([0;0],Ssigma_now)';

t = 2;
mmu_now = dum_regime1(t)*true_mmu1+dum_regime2(t)*true_mmu2;
mmu_lag1 = dum_regime1(t-1)*true_mmu1+dum_regime2(t-1)*true_mmu2;
mmu_lag2 = true_mmu1;
ylag2 = y_zero;
ylag1 = y(:,1);
Ssigma_now = dum_regime1(t)*true_Ssigma1+dum_regime2(t)*true_Ssigma2;
y(:,t) = mmu_now + [true_pphi1 true_pphi2]*([ylag1; ylag2]-[mmu_lag1; mmu_lag2]) + mvnrnd([0;0],Ssigma_now)';

for t = 3:T
	mmu_now = dum_regime1(t)*true_mmu1+dum_regime2(t)*true_mmu2;
	mmu_lag1 = dum_regime1(t-1)*true_mmu1+dum_regime2(t-1)*true_mmu2;
	mmu_lag2 = dum_regime1(t-2)*true_mmu1+dum_regime2(t-2)*true_mmu2;
	ylag2 = y(:,t-2);
	ylag1 = y(:,t-1);
	Ssigma_now = dum_regime1(t)*true_Ssigma1+dum_regime2(t)*true_Ssigma2;
	y(:,t) = mmu_now + [true_pphi1 true_pphi2]*([ylag1; ylag2]-[mmu_lag1; mmu_lag2]) + mvnrnd([0;0],Ssigma_now)';
end

% remove the burnin sections
y = y(:,burnin+1:end);
regimes = regimes(burnin+1:end,1);

save('sim_VAR','y','regimes')
%% one pass of drawing mmu
N = 2;
p = 2;
M = 2;
model.N = 2;
model.p = 2;
model.K = model.N*model.p+1;
model.T = T;
model.M = 2;

condition.regimes = regimes;
condition.dum_mat = [dum_regime1 dum_regime2];
temp_Pphi = [true_pphi1,true_pphi2];
condition.pphi = temp_Pphi(:);

Ssigma_array = zeros(N,N,p);
Ssigma_array(:,:,1) = true_Ssigma1;
Ssigma_array(:,:,2) = true_Ssigma2;
condition.Ssigma_array = Ssigma_array;
condition.y_table = y';

S2star = zeros(T*N,M*N);
y2star = zeros(T*N,1);
for t = 1:T
    S2star(1+N*(t-1):N+N*(t-1),:) = transform_St2star_for_mmu(t,model,condition);
    y2star(1+N*(t-1):N+N*(t-1),:) = transform_yt2star_for_mmu(t,model,condition);
end

% more elaborate Gibbs sampling, now it's simply a univariate regression
mmu_model.T = N*(T-p);
mmu_prior.bbeta_mean = zeros(N*M,1); % how confident about the fake SSR in prior
mmu_prior.bbeta_cov = 1e5*eye(N*M); % the S in the notes
mmu_condition.y = y2star(N*p+1:end,:); % the dependent TN-by-1, T # of obs, N # of variables
mmu_condition.Z = S2star(N*p+1:end,:); % the TN-by-K matrix, K = sum k_n
mmu_condition.Ssigma = 1; % the drown Ssigma matrix from previous step

for r = 1:100
	mmu_draws(:,r) = post_draw_bbeta_indie(mmu_model,mmu_condition,mmu_prior);
end
mean(mmu_draws,2)
prctile(mmu_draws,10,2)
prctile(mmu_draws,100-10,2)

%% Draw pphi
% first construct needed transformed data
condition.mmu = mean(mmu_draws,2);
[pphi_y2star,pphi_Z] = transform_for_pphi(model,condition);

% quick OLS test
(pphi_Z(N*p+1:end,:)'*pphi_Z(N*p+1:end,:))\(pphi_Z(N*p+1:end,:)'*pphi_y2star(N*p+1:end,:))

% more elaborate Gibbs sampling, now it's simply a univariate regression
pphi_model.T = N*(T-p);
pphi_prior.bbeta_mean = zeros(N*N*p,1); % how confident about the fake SSR in prior
pphi_prior.bbeta_cov = 1e5*eye(N*N*p); % the S in the notes
pphi_condition.y = pphi_y2star(N*p+1:end,:); % the dependent TN-by-1, T # of obs, N # of variables
pphi_condition.Z = pphi_Z(N*p+1:end,:); % the TN-by-K matrix, K = sum k_n
pphi_condition.Ssigma = 1; % the drown Ssigma matrix from previous step

R = 1000;
pphi_draws = zeros(N*N*p,R); 
for r = 1:100
	pphi_draws(:,r) = post_draw_bbeta_indie(pphi_model,pphi_condition,pphi_prior);
end
mean(pphi_draws,2)
prctile(pphi_draws,10,2)
prctile(pphi_draws,100-10,2)

%% get Ssigma per regime
condition.mmu = [true_mmu1;true_mmu2]; % change this
pphi_mat = [true_pphi1 true_pphi2];
condition.pphi = pphi_mat(:);
[E_table_forSsigma] = transform_for_Ssigma(model,condition);
E_table_forSsigma = E_table_forSsigma(p+1:end,:);
regimes_forSsigma = regimes(p+1:end,:);

% more elaborate Ssigma model/condition
Ssigma_prior.SSR = 0;
Ssigma_prior.nnu = 0;
Ssigma_condition.bbeta = zeros(N,1);
Ssigma_model.N = 2; 
rows_forSsigma = 1:T-p;
for m = 1:M
	Em = E_table_forSsigma(regimes_forSsigma==m,:);
	Ssigma_model.T = size(Em,1);
	temp = Em';
	Ssigma_condition.y = temp(:);
	Ssigma_condition.Z = ones(Ssigma_model.T*Ssigma_model.N,Ssigma_model.N);
	for r = 1:1000
		% pack things in
		Ssigma_draw(:,:,m,r) = post_draw_Ssigma_indie(Ssigma_model,Ssigma_condition,Ssigma_prior);
	end
end

prctile(Ssigma_draw(2,2,2,:),10)
prctile(Ssigma_draw(2,2,2,:),100-10)
mean(Ssigma_draw,4)