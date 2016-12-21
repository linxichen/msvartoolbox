%% housekeeping
clear;
clc;
close all
addpath(genpath('../'))

%% set hyperparameters
var1 = 0.5;
var2 = 2.2;
cov12 = -0.5;
true_Ssigma = [var1, cov12; cov12, var2];
bbeta1 = 3.5;
bbeta2 = 10;
bbeta3 = -6;

%% Simulate data
T = 300;
eps = mvnrnd(zeros(2,1),true_Ssigma,T);

x1 = rand(T,1);
x2 = rand(T,1);
x3 = rand(T,1);

y1 = bbeta1*x1 + bbeta2*x2 + eps(:,1);
y2 = bbeta3*x3 + eps(:,2);

%% create data input 
Y = [y1,y2]; Y = Y'; y = Y(:);
N = 2;
K = 3;
for t = 1:T
	Z(1,t,:) = [x1(t) x2(t) 0]; %#ok<*SAGROW>
	Z(2,t,:) = [ 0     0    x3(t)];
end
Z = reshape(Z,[T*2,3]);

%% specify my priors and initial things
model.T = T;
model.N = 2;
model.K = 3;
options.burnin = 5e3;
options.ndraws = 5e3;
options.bbeta_init = zeros(K,1);
options.Ssigma_init = eye(N,N);
prior.bbeta_mean = zeros(K,1);
prior.bbeta_cov = 10000*eye(K,K);
prior.SSR = zeros(N,N);
prior.nnu = 0;
draws = MVReg_Gibbs(y,Z,model,prior,options);

%% Get point estimates
post_med_Ssigma = median(draws.Ssigma,3)
true_Ssigma
post_med_bbeta = median(draws.bbeta,3)
prctile(draws.bbeta,100*(0.05/2),3)
prctile(draws.bbeta,100*(1-0.05/2),3)