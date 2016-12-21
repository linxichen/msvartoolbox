function OLSest = VAROLS(data, nlag)

[fullT, N] = size(data);
T = fullT - nlag;
K = N*nlag + 1; % number of parameters in one equation only
Y = data(nlag+1:end,:); % LHS matrix in OLS regression
X = [];
for lag = 1:nlag
	X = [X data(nlag+1-lag:nlag+T-lag,:)];
end
X = [X ones(T,1)];
invXprimeX = inv(X'*X);

OLSest.nobs = T; % usable observations (in sample)
OLSest.nlag = nlag;
OLSest.K = K; % number of parameters in one equation only
OLSest.Pphi = invXprimeX*(X'*Y); %#ok<*MINV> % OLS for AR matrix
OLSest.SSR = (Y-X*OLSest.Pphi)'*(Y-X*OLSest.Pphi); % sum of squared residuals from OLS
OLSest.invXprimeX = inv(X'*X);
OLSest.Ssigma = OLSest.SSR/(T-K);

end
