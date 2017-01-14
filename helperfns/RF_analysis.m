%% loading data
load deliverable1.mat

%% Apply the IRF thing to everything
R = length(draws.pphi);
for i_draw = 1:R
	tmp_Pphi_array = reshape(draws.pphi(:,i_draw),N,N,p,M1);
	Pphi1_array = reshape(tmp_Pphi_array(:,:,:,1),N,N,p);
	Pphi2_array = reshape(tmp_Pphi_array(:,:,:,2),N,N,p);
	IRF1(:,:,:) = OIRF_RFVAR(20,Pphi1_array,draws.Ssigma_array(:,:,1,i_draw),N,p);
	IRF2(:,:,:) = OIRF_RFVAR(20,Pphi2_array,draws.Ssigma_array(:,:,1,i_draw),N,p);
end
