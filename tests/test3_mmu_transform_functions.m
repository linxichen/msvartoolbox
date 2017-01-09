%% housekeeping
clear;
clc;
close all
addpath(genpath('../helperfns'))
rng(1234)

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
