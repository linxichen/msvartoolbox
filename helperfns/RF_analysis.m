%% housekeeping
clear;
clc;
close all
addpath(genpath('../helperfns'))

%% loading data
load deliverable1.mat

%% Apply the IRF thing to everything
R = length(draws.pphi);
horizon = 40;
FEVD_array = zeros(N,N,N,M1,M2,R);
cov_array = zeros(N,N,M1,M2,R);
corr_array = zeros(N,N,M1,M2,R);
OIRF_array = zeros(N,N,horizon+1,M1,M2,R);
parfor i_draw = 1:R
	pphi = draws.pphi(:,i_draw);
	Ssigma_array = draws.Ssigma_array(:,:,:,i_draw);
	[FEVD_array(:,:,:,:,:,i_draw),...
		cov_array(:,:,:,:,i_draw),...
		OIRF_array(:,:,:,:,:,i_draw),...
		corr_array(:,:,:,:,i_draw),...
		] = ...
		RFVAR_var_analysis(horizon,pphi,Ssigma_array,model);
end
save deliverable1_analysis.mat
