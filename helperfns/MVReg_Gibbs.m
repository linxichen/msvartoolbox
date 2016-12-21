function draws = MVReg_Gibbs(y,Z,model,prior,options)
% Estimate a multivariate model y = Z*bbeta + E
% For model with N depedents and potentially different regressors for each
% equations in the following form: For the n-th variable at observation t
% y_n,t = Z_t * beta + e_t.
% where y_n,t is 1-by-1 scalar, Z_t is N-by-K regressors, K = sum k_n, k_n
% the number of regressors for n-th equation.
% Input: y is N*T-by-1 vector [y1 y2  ... yT]' where y_t is N-by-1
% Z is N*T-by-K  stack data matrix

%% Unpack model settings
T = model.T;
N = model.N;
K = model.K;

%% Basic error checking
if isequal(size(y),[T*N,1])
else
	error('dependent vector is not expected size');
end

if isequal(size(Z),[T*N,K])
else
	error('regressor matrix is not expected size');
end

%% pack condition with initials;
condition.Ssigma = options.Ssigma_init;
condition.bbeta = options.bbeta_init;

%% put data into the conditionals
condition.y = y;
condition.Z = Z;

% draw bbeta, Ssigma
burnin = options.burnin;
ndraws = options.ndraws; % post burnin draws
count = 0;

% Burnin draws
while count < burnin
	bbeta_draw = post_draw_bbeta_indie(model,condition,prior);
	condition.bbeta = bbeta_draw;
	Ssigma_draw = post_draw_Ssigma_indie(model,condition,prior);
	condition.bbeta = Ssigma_draw;
	count = count + 1;
	if mod(count,100) == 0
		fprintf('Drawed %d for burning.\n',count);
	end
	if mod(count,666) == 0
		fprintf('Drawed %d for burning.\n',count);
	end
end

% Draws post burnin
draws.bbeta = zeros(K,1,ndraws);
draws.Ssigma = zeros(N,N,ndraws);
for i_draw = 1:ndraws
	bbeta_draw = post_draw_bbeta_indie(model,condition,prior);
	condition.bbeta = bbeta_draw;
	draws.bbeta(:,i_draw) = bbeta_draw;
	
	Ssigma_draw = post_draw_Ssigma_indie(model,condition,prior);
	condition.bbeta = Ssigma_draw;
	draws.Ssigma(:,:,i_draw) = Ssigma_draw;
	
	%% show something I don't panic
	if mod(i_draw,100) == 0
		fprintf('Drawed %d for storage.\n',i_draw);
	end
end

end