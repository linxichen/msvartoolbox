%% housekeeping
close all;
clc;
clear;
addpath('../helperfns');
rng(12345);

%% generate all inputs
tmp_N = 2;
tmp_T = 1e4;
tmp_p = 3;
tmp_num_AR_regimes = 2;
tmp_P1 = [0.94 0.06; 1-0.75 0.75];
tmp_P2 = [0.97 1-0.97; 1-0.97 0.97];
tmp_mmu = [];
for i_regime = 1:tmp_num_AR_regimes
	tmp_mmu = [tmp_mmu; rand([tmp_N,1]); 20*rand([tmp_N,1])];
end

for i_lag = 1:tmp_p
	for i_regime = 1:tmp_num_AR_regimes
		tmp_Pphi(:,:,i_lag,i_regime) = diag(1/tmp_p*rand([tmp_N,1]));
		if i_regime == 2
			tmp_mat = 5*rand(tmp_N,tmp_N);
		else
			tmp_mat = 0.5*rand(tmp_N,tmp_N);
		end
		tmp_Ssigma_array(:,:,i_regime) = tmp_mat'*tmp_mat;
	end
end

%% input loading into input
model.N          = tmp_N;
model.p          = tmp_p;
T                = tmp_T;
mmu              = tmp_mmu;
model.breakdate  = ceil(0.5*T);
pphi             = tmp_Pphi(:);
Ssigma_array     = tmp_Ssigma_array;
transprob{1} = tmp_P1;
transprob{2} = tmp_P2;

%% callaing the simulation
% clear tmp*;
[y,regimes_mat] = simulate_MSVAR(T,mmu,pphi,Ssigma_array,transprob,model);
save sim_for_mmu
plot(1:T,y(:,1)',1:T,regimes_mat(:,2)');

%% a bit of testing
y_pre = y(1:model.breakdate,:);
regime1_pre = regimes_mat(1:model.breakdate,1);
y_pre_regime1 = y_pre(regime1_pre==1,:);
y_pre_regime2 = y_pre(regime1_pre==2,:);

y_post = y(model.breakdate+1:end,:);
regime1_post = regimes_mat(model.breakdate+1:end,1);

y_post_regime1 = y_post(regime1_post==1,:);
y_post_regime2 = y_post(regime1_post==2,:);

mmu_array = reshape(mmu,model.N,2,2);

mean(y_pre_regime1)'
mean(y_post_regime1)'
mean(y_pre_regime2)'
mean(y_post_regime2)'
