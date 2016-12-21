function [y2star,Z] = transform_for_pphi(model,condition)
%% unpack things
T = model.T;
N = model.N;
M = model.M;
p = model.p;

%% unload conditional info
regimes = condition.regimes;
mmu_array = reshape(condition.mmu,[N,M]);
y2star = zeros(T*N,1);
Z = zeros(T*N,p*N^2);
invL_array = condition.Ssigma_array;
for m = 1:M
    L = chol(condition.Ssigma_array(:,:,m),'lower');
    invL_array(:,:,m) = L^(-1);
end

%% first find de-mean 
ystar_table = zeros(T,N);
for t = 1:T
	ystar_table(t,:) = condition.y_table(t,:) - mmu_array(:,regimes(t))'; 
end

xstar_table = lagmatrix(ystar_table,1:p);

for t = 1:T
	% we all need the current Ssigma
    invL = invL_array(:,:,regimes(t));
    
	% the easy one is the lhs variable
	y2star(1+N*(t-1):N+N*(t-1),:) = invL*ystar_table(t,:)';
	
	% fill the Z, which is T*N-by-p*N*N	
	Z(1+N*(t-1):N+N*(t-1),:) = kron(xstar_table(t,:),invL);
end

end