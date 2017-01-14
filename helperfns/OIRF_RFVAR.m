function OIRF = OIRF_RFVAR(horizon,Pphi_array,Ssigma,N,p)
% The usual Orthogonalized IRF for a demeaned VAR
% The initial shock size is normalized to be unity
% Pphi_array = reshape(pphi,N,N*p);
N = size(Pphi_array,1);
p = size(Pphi_array,3);
L = chol(Ssigma,'lower');
% [V,D] = eig(Ssigma);
structual_shocks = zeros(N,horizon);

impulsevec = zeros(N,N,horizon);
controlvec = zeros(N,N,horizon);

% on impact, apply shock
% if whichvar > 0
% 	structual_shocks(whichvar,1) = 1/L(whichvar,whichvar);
% else
% 	structual_shocks = zeros(N,horizon);
% end

% convert to reduced form shocks
pastresponse = zeros(N,N,p);
RF_shocks = L/diag(diag(L));
impulsevec(:,:,1) = RF_shocks;

% first p lags involve pre-shock states
for t = 2:horizon
	pastresponse(:,:,2:p)  = pastresponse(:,:,1:p-1);
	pastresponse(:,:,1) = impulsevec(:,:,t-1);
	tmp = zeros(N,N);
	for j = 1:p
		tmp = tmp + Pphi_array(:,:,j)*pastresponse(:,:,j);
	end
	impulsevec(:,:,t) = tmp;
end

OIRF = impulsevec - controlvec;
% for i = 1:N
% 	OIRF(:,i,:) = tmp_OIRF(:,i,:)/L(i,i);
% end

end
