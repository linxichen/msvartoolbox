function est = RFVAR_OLS(y_table,X_table,N)
T = size(y_table,1);
K = size(X_table,2);
X = zeros(N*T,N*K);
for t = 1:T
	X(1+N*(t-1):N+N*(t-1),:) = kron(X_table(t,:),eye(N));
end
yprime = y_table';
Y = yprime(:);
B = (X'*X)\(X'*Y);
E = Y-X*B;
E = reshape(E,[N T]);

est.coeff = (X'*X)\(X'*Y);
est.Ssigma = E*E'/(T-K);

end