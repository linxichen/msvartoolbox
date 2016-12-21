function E = transform_for_Ssigma(model,condition)
%% unpack things
T = model.T;
N = model.N;
M = model.M;
p = model.p;

%% unload conditional info
regimes = condition.regimes;
mmu_array = reshape(condition.mmu,[N,M]);

%% first find de-mean 
ystar_table = zeros(T,N);
for t = 1:T
	ystar_table(t,:) = condition.y_table(t,:) - mmu_array(:,regimes(t))'; 
end

xstar_table = lagmatrix(ystar_table,1:p);
E = zeros(T*N,1);
for t = 1:T
	E(1+N*(t-1):N+N*(t-1)) = ystar_table(t,:)'-(kron(xstar_table(t,:),eye(N)))*condition.pphi;
end
E = reshape(E,N,T);
E = E';
end