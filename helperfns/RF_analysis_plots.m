%% housekeeping
clear;
clc;
close all
addpath(genpath('../helperfns'))
load('deliverable1_analysis.mat')

%% Create tables and graphs after analysis
% Get the steady state numbers
mmu_array = reshape(draws.mmu,[N 2 M1 R]);
%%
ar_regime = 1;
sb_period = 1;
for i_var = 1:N
	i_var
	prctile(squeeze(mmu_array(i_var,sb_period,ar_regime,:)),16)
	prctile(squeeze(mmu_array(i_var,sb_period,ar_regime,:)),50)
	prctile(squeeze(mmu_array(i_var,sb_period,ar_regime,:)),84)
	disp('==========')
end

% Get the steady state numbers
ar_regime = 2;
sb_period = 1;
for i_var = 1:N
	i_var
	prctile(squeeze(mmu_array(i_var,sb_period,ar_regime,:)),16)
	prctile(squeeze(mmu_array(i_var,sb_period,ar_regime,:)),50)
	prctile(squeeze(mmu_array(i_var,sb_period,ar_regime,:)),84)
	disp('==========')
end

%% Get FEVD tables first, low AR and low COV regimes
ar_regime = 1;
cov_regime = 1;

table_16 = zeros(N,N);
table_50 = zeros(N,N);
table_84 = zeros(N,N);
for i = 1:N
	for j = 1:N
		table_16(i,j) = prctile(squeeze(corr_array(i,j,ar_regime,cov_regime,:)),16);
		table_50(i,j) = prctile(squeeze(corr_array(i,j,ar_regime,cov_regime,:)),50);
		table_84(i,j) = prctile(squeeze(corr_array(i,j,ar_regime,cov_regime,:)),84);
	end
end

%% Get FEVD tables first, low AR and low COV regimes
ar_regime = 2;
cov_regime = 1;

table_16 = zeros(N,N);
table_50 = zeros(N,N);
table_84 = zeros(N,N);
for i = 1:N
	for j = 1:N
		table_16(i,j) = prctile(squeeze(corr_array(i,j,ar_regime,cov_regime,:)),16);
		table_50(i,j) = prctile(squeeze(corr_array(i,j,ar_regime,cov_regime,:)),50);
		table_84(i,j) = prctile(squeeze(corr_array(i,j,ar_regime,cov_regime,:)),84);
	end
end

%% Get FEVD tables first, low AR and low COV regimes
ar_regime = 1;
cov_regime = 2;

table_16 = zeros(N,N);
table_50 = zeros(N,N);
table_84 = zeros(N,N);
for i = 1:N
	for j = 1:N
		table_16(i,j) = prctile(squeeze(corr_array(i,j,ar_regime,cov_regime,:)),16);
		table_50(i,j) = prctile(squeeze(corr_array(i,j,ar_regime,cov_regime,:)),50);
		table_84(i,j) = prctile(squeeze(corr_array(i,j,ar_regime,cov_regime,:)),84);
	end
end

%% Get FEVD tables first, low AR and low COV regimes
ar_regime = 2;
cov_regime = 2;

table_16 = zeros(N,N);
table_50 = zeros(N,N);
table_84 = zeros(N,N);
for i = 1:N
	for j = 1:N
		table_16(i,j) = prctile(squeeze(corr_array(i,j,ar_regime,cov_regime,:)),16);
		table_50(i,j) = prctile(squeeze(corr_array(i,j,ar_regime,cov_regime,:)),50);
		table_84(i,j) = prctile(squeeze(corr_array(i,j,ar_regime,cov_regime,:)),84);
	end
end

%% serious IRF plotting
width = 8;     % Width in inches
height = 5;    % Height in inches
alw = 1.1;    % AxesLineWidth
fsz = 11;      % Fontsize
lw = 1.5;      % LineWidth
msz = 8;       % MarkerSize

h = figure('Units','inches',...
	'Position',[0 0 width height]);

plot_horizon = 15;
subplot(1,2,1)
plot(0:plot_horizon,prctile(squeeze(OIRF_array(1,3,1:plot_horizon+1,1,1,:)),50,2),'b-',...
	0:plot_horizon,prctile(squeeze(OIRF_array(1,3,1:plot_horizon+1,2,1,:)),50,2),'r-.','LineWidth',lw,'MarkerSize',msz ...
	)
legend('Expansion','Recession')
xlabel('Periods From Imparct','FontSize',10)
ylabel('Percentage')
ylabh = get(gca,'yLabel');
set(ylabh,'Position',get(ylabh,'Position') - [2 0 0])
title('Inventory Investment')

subplot(1,2,2)
plot(0:plot_horizon,prctile(squeeze(OIRF_array(1,2,1:plot_horizon+1,1,1,:)),50,2),'b-',...
	0:plot_horizon,prctile(squeeze(OIRF_array(1,2,1:plot_horizon+1,2,1,:)),50,2),'r-.','LineWidth',lw,'MarkerSize',msz...
	)
xlabel('Periods From Imparct')
title('GDP Growth')
legend('Expansion','Recession')
print('OIRF_sales','-depsc2','-r300');


%% plotting
figure
subplot(2,1,1)
plot(0:horizon,prctile(squeeze(OIRF_array(1,3,:,1,1,:)),16,2),...
	0:horizon,prctile(squeeze(OIRF_array(1,3,:,1,1,:)),50,2),...
	0:horizon,prctile(squeeze(OIRF_array(1,3,:,1,1,:)),84,2)...
	)
subplot(2,1,2)
plot(0:horizon,prctile(squeeze(OIRF_array(1,3,:,2,1,:)),16,2),...
	0:horizon,prctile(squeeze(OIRF_array(1,3,:,2,1,:)),50,2),...
	0:horizon,prctile(squeeze(OIRF_array(1,3,:,2,1,:)),84,2)...
	)



figure
subplot(2,1,1)
plot(0:horizon,prctile(squeeze(OIRF_array(3,1,:,1,1,:)),16,2),...
	0:horizon,prctile(squeeze(OIRF_array(3,1,:,1,1,:)),50,2),...
	0:horizon,prctile(squeeze(OIRF_array(3,1,:,1,1,:)),84,2)...
	)
subplot(2,1,2)
plot(0:horizon,prctile(squeeze(OIRF_array(3,1,:,2,1,:)),16,2),...
	0:horizon,prctile(squeeze(OIRF_array(3,1,:,2,1,:)),50,2),...
	0:horizon,prctile(squeeze(OIRF_array(3,1,:,2,1,:)),84,2)...
	)

figure
subplot(2,1,1)
plot(0:horizon-1,prctile(squeeze(IRF1(1,1,:,:)),16,2),...
	0:horizon-1,prctile(squeeze(IRF1(1,1,:,:)),50,2),...
	0:horizon-1,prctile(squeeze(IRF1(1,1,:,:)),84,2)...
	)
subplot(2,1,2)
plot(0:horizon-1,prctile(squeeze(IRF2(1,1,:,:)),16,2),...
	0:horizon-1,prctile(squeeze(IRF2(1,1,:,:)),50,2),...
	0:horizon-1,prctile(squeeze(IRF2(1,1,:,:)),84,2)...
	)

figure
subplot(2,1,1)
plot(0:horizon-1,prctile(squeeze(IRF1(1,2,:,:)),16,2),...
	0:horizon-1,prctile(squeeze(IRF1(1,2,:,:)),50,2),...
	0:horizon-1,prctile(squeeze(IRF1(1,2,:,:)),84,2)...
	)
subplot(2,1,2)
plot(0:horizon-1,prctile(squeeze(IRF2(1,2,:,:)),16,2),...
	0:horizon-1,prctile(squeeze(IRF2(1,2,:,:)),50,2),...
	0:horizon-1,prctile(squeeze(IRF2(1,2,:,:)),84,2)...
	)

figure
subplot(2,1,1)
plot(0:horizon-1,prctile(squeeze(IRF1(3,2,:,:)),16,2),...
	0:horizon-1,prctile(squeeze(IRF1(3,2,:,:)),50,2),...
	0:horizon-1,prctile(squeeze(IRF1(3,2,:,:)),84,2)...
	)
subplot(2,1,2)
plot(0:horizon-1,prctile(squeeze(IRF2(3,2,:,:)),16,2),...
	0:horizon-1,prctile(squeeze(IRF2(3,2,:,:)),50,2),...
	0:horizon-1,prctile(squeeze(IRF2(3,2,:,:)),84,2)...
	)