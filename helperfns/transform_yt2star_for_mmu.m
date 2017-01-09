function yt2star = transform_yt2star_for_mmu(t,model,condition)
% compute the S** matrix in generating for mmu
p = model.p;
M1 = model.M1;
N = model.N;
T = model.T;

% transform pphi to pphi_mat
Pphi_array = reshape(condition.pphi,[N,N,p,M1]);
m1 = condition.regimes(t,1);

if t <= p
	yt2star = NaN(N,1);
elseif t > T
	error('t out of sample size T');
else
	% getting dummy varialbes and Ssigma now
	m2 = condition.regimes(t,2);
	Ssigma_current = condition.Ssigma_array(:,:,m2);
	L = chol(Ssigma_current,'lower');

	% actually computing the N-by-1 yt vector
	ytstar = condition.y_table(t,:)';
	for j = 1:p
		ytstar = ytstar - Pphi_array(:,:,j,m1)*condition.y_table(t-j,:)';
	end
	yt2star = L\ytstar;
end

end
