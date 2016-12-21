function [Pphi, Ssigma] = RFdraw(OLSest,prior)
% Draw from posterior distribution of parameters (Pphi, Ssigma)

% Preliminaries
S_hat = OLSest.SSR;
T = OLSest.nobs;
K = OLSest.K;
invXprimeX = OLSest.invXprimeX;
Pphi_hat = OLSest.Pphi;

% By default, draw from Zeller's inproper prior which has no parameters
Ssigma = iwishrnd(S_hat,T-K);
Pphi = mvnrnd(Pphi_hat(:),kron(Ssigma,invXprimeX));
Pphi = reshape(Pphi,size(Pphi_hat));

end