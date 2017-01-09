function [y2star,Z] = transform_for_pphi(model,condition)
%% unpack things
T = model.T;
N = model.N;
M1 = model.M1;
M2 = model.M2;
p = model.p;

%% unload conditional info
regimes = condition.regimes;
mmu_array = reshape(condition.mmu,[N,2,M1]);
y2star = zeros(T*N,1);
Z = zeros(T*N,p*N^2*M1);
invL_array = condition.Ssigma_array;
for m = 1:M2
    L = chol(condition.Ssigma_array(:,:,m),'lower');
    invL_array(:,:,m) = L^(-1);
end

%% first find de-mean
ystar_table = zeros(T,N);
for t = 1:T
	if t <= model.breakdate
		ystar_table(t,:) = condition.y_table(t,:) - mmu_array(:,1,regimes(t,1))';
	else
		ystar_table(t,:) = condition.y_table(t,:) - mmu_array(:,2,regimes(t,1))';
	end
end

for t = p+1:T
	% we all need the current Ssigma
	invL = invL_array(:,:,regimes(t,2));

	% the easy one is the lhs variable
	y2star(1+N*(t-1):N+N*(t-1),:) = invL*ystar_table(t,:)';

	% fill the Z, which is T*N-by-p*N*N
	z2star_mt = zeros(N,N^2,p,M1);
	for m1 = 1:M1
		for j = 1:p
			z2star_mt(:,:,j,m1) = (regimes(t,1)==m1)*kron(ystar_table(t-j,:),eye(N));
		end
	end
	Z(1+N*(t-1):N+N*(t-1),:) = invL*reshape(z2star_mt,[N,p*M1*N^2]);
end

end
