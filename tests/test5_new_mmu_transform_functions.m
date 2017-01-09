%% housekeeping
clear;
clc;
close all
addpath(genpath('../helperfns'))
rng(1234)

load sim_for_mmu

%% Test everything
model.M1 = 2;
model.M2 = 2;
model.T = length(y);

condition.regimes = regimes_mat;
condition.y_table = y;
condition.AR_dum_mat  = [regimes_mat(:,1)==1 regimes_mat(:,1)==2];
condition.cov_dum_mat = [regimes_mat(:,2)==1 regimes_mat(:,2)==2];
condition.pphi = pphi;
breakdate = ceil(0.5*T);
condition.breakdate = ceil(0.5*T);
condition.break_dum_mat = [[ones(breakdate,1);zeros(T-breakdate,1)],...
	[zeros(breakdate,1);ones(T-breakdate,1)]];

Ssigma_array = tmp_Ssigma_array;
condition.Ssigma_array = Ssigma_array;

N = model.N;
for t = model.p+1:T
	S2star(1+N*(t-1):N+N*(t-1),:) = transform_St2star_for_mmu(t,model,condition);
	y2star(1+N*(t-1):N+N*(t-1),:) = transform_yt2star_for_mmu(t,model,condition);
end

Y = y2star; % the dependent TN-by-1, T # of obs, N # of variables
Z = S2star; % the TN-by-K matrix, K = sum k_n
est_mmu = (Z'*Z)\(Z'*Y);
diff = est_mmu - mmu