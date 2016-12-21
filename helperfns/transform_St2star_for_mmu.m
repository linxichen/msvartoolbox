function St2star = transform_St2star_for_mmu(t,model,condition)
% compute the S** matrix in generating for mmu
% assumes:
%   pphi = vec([Pphi_1 ... Pphi_p]) where Pphi_j is the lagged VAR coeff
%   dum_mat is the T-by-M matrix where the m-th column is regimes == m
p = model.p; % number of lags
M = model.M; % number regimes
N = model.N; % number of variable/equations
T = model.T; % sample size in time dimension

% transform pphi to pphi_mat
Pphi_mat = reshape(condition.pphi,[N,p*N]);

if t <= p
	St2star = NaN(N,M*N);
elseif t > T
	error('t out of sample size T');
else
	% getting dummy varialbes and Ssigma now
	dum_mat = condition.dum_mat;
	m = condition.regimes(t);
	Ssigma_current = condition.Ssigma_array(:,:,m);
	L = chol(Ssigma_current,'lower');
	
	% actually computing the N-by-M*N St matrix
	dum_sec = dum_mat(t-p:t,:);
	dum_flip = flipud(dum_sec);
	St2star = L\([eye(N),-Pphi_mat]*kron(dum_flip,eye(N)));
end