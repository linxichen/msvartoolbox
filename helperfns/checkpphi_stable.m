function stable = checkpphi_stable(pphi,N,p)
% check whether pphi = vec([Pphi_1 ... Pphi_p]) is a stable VAR
Pphi_mat = reshape(pphi,N,N*p);
% construct the companion matrix as you write VAR(p) into VAR(1)
F = [Pphi_mat; kron(eye(p-1),eye(N)), zeros(N*(p-1),N)];
eigenvalues = eig(F,'vector');
stable = max(abs(eigenvalues)) < 1;
end
