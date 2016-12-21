function yt2star = transform_yt2star_for_mmu(t,model,condition)
% compute the S** matrix in generating for mmu
p = model.p;
M = model.M;
N = model.N;
T = model.T;

% transform pphi to pphi_mat
Pphi_mat = reshape(condition.pphi,[N,p*N]);

if t <= p
	yt2star = NaN(N,1);
elseif t > T
	error('t out of sample size T');
else
	% getting dummy varialbes and Ssigma now
	m = condition.regimes(t);
	Ssigma_current = condition.Ssigma_array(:,:,m);
	L = chol(Ssigma_current,'lower');
    
	% actually computing the N-by-1 yt vector
    ytstar = condition.y_table(t,:)';
    for j = 1:p
        ytstar = ytstar - Pphi_mat(:,1+N*(j-1):N+N*(j-1))*condition.y_table(t-j,:)';
    end
	yt2star = L\ytstar;
end

end