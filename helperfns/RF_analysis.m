%% loading data
load deliverable1.mat

%% Apply the IRF thing to everything
R = length(draws.pphi);
horizon = 40;
IRF1 = zeros(N,N,horizon,R);
IRF2 = IRF1;
for i_draw = 1:R
	tmp_Pphi_array = reshape(draws.pphi(:,i_draw),N,N,p,M1);
	Pphi1_array = reshape(tmp_Pphi_array(:,:,:,1),N,N,p);
	Pphi2_array = reshape(tmp_Pphi_array(:,:,:,2),N,N,p);
	IRF1(:,:,:,i_draw) = OIRF_RFVAR(horizon,Pphi1_array,draws.Ssigma_array(:,:,1,i_draw),N,p);
	IRF2(:,:,:,i_draw) = OIRF_RFVAR(horizon,Pphi2_array,draws.Ssigma_array(:,:,1,i_draw),N,p);
end

%% plotting
subplot(2,1,1)
plot(0:horizon-1,prctile(squeeze(IRF1(2,1,:,:)),16,2),...
	0:horizon-1,prctile(squeeze(IRF1(2,1,:,:)),50,2),...
	0:horizon-1,prctile(squeeze(IRF1(2,1,:,:)),84,2)...
	)
subplot(2,1,2)
plot(0:horizon-1,prctile(squeeze(IRF2(2,1,:,:)),16,2),...
	0:horizon-1,prctile(squeeze(IRF2(2,1,:,:)),50,2),...
	0:horizon-1,prctile(squeeze(IRF2(2,1,:,:)),84,2)...
	)