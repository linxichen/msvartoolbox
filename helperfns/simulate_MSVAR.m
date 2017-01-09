function [obs,reg_mat] = simulate_MSVAR(T,mmu,pphi,Ssigma_array,transprob,model)
% Simulate data from a markov switching VAR in the demeaned form
%   1. model specifies the dimensions, number of vars, number of markov chains,
%      number of possible regimes
%   2. mmu,pphi,Ssigma_array are the true values used to simulate
%   3. transporb contains the true transition probability for all markov
%      chains
%   4. T is the length of simulation
%===============================================================================
% The function body documentations goes here
%===============================================================================
%% load dimension
N          = model.N; % number of variables
p          = model.p; % number of lags
breakdate  = model.breakdate;
M1         = length(transprob{1});
num_chains = length(transprob);
Pphi_array = reshape(pphi,N,N,p,M1);
mmu_array  = reshape(mmu,N,2,M1);
burnin     = 50*T;
bigT       = T+burnin;

%% simulate markov chains
regimes = zeros(bigT,num_chains);
for i_chain = 1:num_chains
	PPP = transprob{i_chain};
	MMM = length(PPP);
	regimes(:,i_chain) = markovsim(1:MMM,PPP,1,bigT)';
end
break_sections = 1+[zeros(burnin+breakdate,1);ones(T-breakdate,1)];

%% Initialize the simulation
y_lag = zeros(N,p);
mmu_lag = kron(ones(1,p),mmu_array(:,1,1));
demeaned_lag = y_lag - mmu_lag;

% convert to reduced form shocks
y = zeros(N,T+burnin);

% main simulation body
zeros_N = zeros(N,1);
for t = 1:T+burnin
	% find current things
	mmu_now = mmu_array(:,break_sections(t),regimes(t,1));
	Pphi_mat_now = reshape(Pphi_array(:,:,:,regimes(t,1)),N,N*p);
	Ssigma_now = Ssigma_array(:,:,regimes(t,2));

	% simulation step
	y(:,t) = mmu_now + Pphi_mat_now*demeaned_lag(:) + mvnrnd(zeros_N,Ssigma_now)';

	% prepare for tomorrow
	mmu_lag = [mmu_now mmu_lag(:,1:p-1)];
	y_lag = [y(:,t) y_lag(:,1:p-1)];
	demeaned_lag = y_lag - mmu_lag;
end

% leave only the last T samples
obs = y(:,end-(T-1):end)';
reg_mat = regimes(end-(T-1):end,:);

end