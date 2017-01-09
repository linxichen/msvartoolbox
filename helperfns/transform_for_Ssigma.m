function E = transform_for_Ssigma(model,condition)
%% unpack things
T = model.T;
N = model.N;
M1 = model.M2;
M2 = model.M2;
p = model.p;

%% unload conditional info
regimes = condition.regimes;
mmu_array = reshape(condition.mmu,[N,2,M1]);
pphilong_array = reshape(condition.pphi,[N,N*p,M1]);

%% first find de-mean
ystar_table = zeros(T,N);
for t = 1:T
	if t <= model.breakdate
		ystar_table(t,:) = condition.y_table(t,:) - mmu_array(:,1,regimes(t,1))';
	else
		ystar_table(t,:) = condition.y_table(t,:) - mmu_array(:,2,regimes(t,1))';
	end
end

xstar_table = lagmatrix(ystar_table,1:p);
E = zeros(T*N,1);
for t = 1:T
	pphilong_now = pphilong_array(:,:,regimes(t,1));
	E(1+N*(t-1):N+N*(t-1)) = ystar_table(t,:)' ...
		-(kron(xstar_table(t,:),eye(N)))*pphilong_now(:);
end
E = reshape(E,N,T);
E = E';
end
