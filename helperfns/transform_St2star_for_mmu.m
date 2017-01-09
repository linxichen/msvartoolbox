function St2star = transform_St2star_for_mmu(t,model,condition)
% compute the S** matrix in generating for mmu
% assumes:
%   pphi = vec(Pphi) where Pphi is N-by-N-by-p-by-M1 that stores VAR
%   coeffs. Pphi(:,:,j,m_1) = j-th lag VAR coeff in AR_regime m_1
%   dum_mat is the T-by-M matrix where the m-th column is AR_regimes == m
%   now regimes contains two chains, one telling me the intercept, AR
%   coefficients and the other the covariance matrix

%% load stuff
p = model.p; % number of lags
M1 = model.M1; % number AR regimes
M2 = model.M2; % number covariacne regimes
N = model.N; % number of variable/equations
T = model.T; % sample size in time dimension
AR_dum_mat = condition.AR_dum_mat;

% transform pphi to pphi_mat
% Pphi_mat = reshape(condition.pphi,[N,p*N*M1]);
Pphi_array = reshape(condition.pphi,[N,N,p,M1]);
Pphi_lag = [];
for j = 1:p
	Pphi_lag = [Pphi_lag,squeeze(Pphi_array(:,:,j,condition.regimes(t,1)))];
end

current_m1 = condition.regimes(t,1);
current_m2 = condition.regimes(t,2);
break_dum_mat = condition.break_dum_mat;
break_periods = size(break_dum_mat,2);

%% Main Construction Body Done by the Smart Me in the Past
%==============================================================================%
% for details refer to the accompanying notes dummy_VAR.pdf e.g. what's St2star
%==============================================================================%
if t <= p
	St2star = NaN(N,M1*N*break_periods);
elseif t > T
	error('t out of sample size T');
else
	% getting dummy varialbes and Ssigma now
	Ssigma_current = condition.Ssigma_array(:,:,current_m2);
	L = chol(Ssigma_current,'lower');

	% Matrix algebra gymnastics that's probably not worth it
	% but it's fun and cool to do anyway, and somehow it's the only thing
	% that I can think of
	AR_dum_sec = AR_dum_mat(t-p:t,:);
	AR_dum_flip = flipud(AR_dum_sec);
	break_dum_sec = break_dum_mat(t-p:t,:);
	break_dum_flip = flipud(break_dum_sec);
	bigmess = kron(AR_dum_flip,ones(1,break_periods)) ...
		.* kron(ones(1,M1),break_dum_flip);

	St2star = L\([eye(N),-Pphi_lag]*kron(bigmess,eye(N)));
end
