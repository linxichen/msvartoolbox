function OIRF = OIRF_RFVAR(whichvar,horizon,pphi,Ssigma,N,p)
% The usual Orthogonalized IRF for a demeaned VAR
% The initial shock size is normalized to be unity
Pphi_mat = reshape(pphi,N,N*p);
L = chol(Ssigma,'lower');
structual_shocks = zeros(N,horizon);

impulsevec = zeros(N,horizon);
controlvec = zeros(N,horizon);

% on impact, apply shock
if whichvar > 0
	structual_shocks(whichvar,1) = 1/L(whichvar,whichvar);
else
	structual_shocks = zeros(N,horizon);
end

% convert to reduced form shocks
RF_shocks = L*structual_shocks;
impulsevec(:,1) = RF_shocks(:,1);

% first p lags involve pre-shock states
for t = 2:horizon
	pastresponse  = impulsevec(:,max(1,t-p):t-1);
	alllags = [pastresponse,zeros(N,p-size(pastresponse,2))];
	impulsevec(:,t) = Pphi_mat*alllags(:) + RF_shocks(:,t);
end

OIRF = impulsevec - controlvec;

end