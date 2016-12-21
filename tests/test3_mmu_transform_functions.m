%% housekeeping
clear;
clc;
close all
addpath(genpath('../helperfns'))
rng(1234)

%% set hyperparameters
var1 = 0.5;
var2 = 2.2;
cov12 = -0.5;
true_Ssigma1 = [var1, cov12; cov12, var2];
true_Ssigma2 = rand(2,2).*true_Ssigma1;
true_Ssigma2 = true_Ssigma2*true_Ssigma2';
true_pphi1 = 0.5*eye(2);
true_pphi2 = 0.4*eye(2);
true_pphi = [true_pphi1;true_pphi2];
true_mmu1 = rand(2,1);
true_mmu2 = rand(2,1);

%% Simulate data
T = 300;
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

pphi_mat = reshape(true_pphi,2,4);
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

%% Test St2star

N = 2;
p = 2;
model.N = 2;
model.p = 2;
model.K = model.N*model.p+1;
model.T = T;
model.M = 2;

condition.regimes = regimes;
condition.dum_mat = [dum_regime1 dum_regime2];
condition.Pphi_mat = [true_pphi1 true_pphi2];

Ssigma_array = zeros(N,N,p);
Ssigma_array(:,:,1) = true_Ssigma1;
Ssigma_array(:,:,2) = true_Ssigma2;
condition.Ssigma_array = Ssigma_array;

t = 100;
St2star = transform_St2star_for_mmu(t,model,condition)
L = chol(dum_regime1(t)*true_Ssigma1+dum_regime2(t)*true_Ssigma2,'lower');

true_St2star_1 = L\(eye(N)*dum_regime1(t) ...
    - (dum_regime1(t-1)*true_pphi1) ...
    - (dum_regime1(t-2)*true_pphi2));

true_St2star_2 = L\(eye(N)*dum_regime2(t) ...
    - (dum_regime2(t-1)*true_pphi1) ...
    - (dum_regime2(t-2)*true_pphi2));
disp([true_St2star_1 true_St2star_2]);

%% Test y2star
condition.y_table = y';
t = 50;
L = chol(dum_regime1(t)*true_Ssigma1+dum_regime2(t)*true_Ssigma2,'lower');
transform_yt2star_for_mmu(t,model,condition)
true_yt2star = L\(eye(N)*condition.y_table(t,:)' ...
    -true_pphi1*condition.y_table(t-1,:)' ...
    -true_pphi2*condition.y_table(t-2,:)' )