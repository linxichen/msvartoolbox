function sim_data = simulate_RFVAR(horizon,pphi,Ssigma,N,p)
% The usual Orthogonalized IRF for a demeaned VAR
% The initial shock size is normalized to be unity
Pphi_mat = reshape(pphi,N,N*p);

% convert to reduced form shocks
sim_data(:,1) = mvnrnd(zeros(N,1),Ssigma)';

% first p lags involve pre-shock states
for t = 2:horizon
	pastresponse  = sim_data(:,max(1,t-p):t-1);
	alllags = [pastresponse,zeros(N,p-size(pastresponse,2))];
	sim_data(:,t) = Pphi_mat*alllags(:) + mvnrnd(zeros(N,1),Ssigma)';
end


end