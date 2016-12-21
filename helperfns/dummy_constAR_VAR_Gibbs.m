function draws = dummy_constAR_VAR_Gibbs(Y_table,S,model,prior,options)
% Estimate a multivariate model y = Z*bbeta + E
% For model with N depedents and potentially different regressors for each
% equations in the following form: For the n-th variable at observation t
% y_n,t = z'_n,t * beta_n + e_n,t.
% where y_n,t is 1-by-1 scalar, z'_n,t is 1-by-k_n regressors, K = sum k_n, k_n
% 
% 
% the number of regressors for n-th equation.
% Input: y_table is T-by-N matrix where the columns are variables
% Z_table is T-by-K  regressor table where K = N*p+1 are number of regressor for
% each equation (assume each equation has the same set of regressors)
% 
%% Unpack model settings
T = model.T; % this is the all observations, first p rows won't be usable
N = model.N;
M = max(S);
p = model.p;
K = N*p+1;

%% Basic error checking
if isequal(size(y),[T*N,1])
else
	error('dependent vector is not expected size');
end

if isequal(size(Z),[T*N,K])
else
	error('regressor matrix is not expected size');
end

if isequal(size(options.pphi_init),[M*(N*p+1),1])
else
	error('regressor matrix is not expected size');
end

%% create dummies and lagged table of regressors
dum_mat = S == 1;
for m = 2:M
	dum_mat = [dum_mat,S==m]; %#ok<AGROW>
end
ylag_table = lagmatrix(Y_table,1:p);
X_table = [ones(T-p,1) ylag_table(p+1:end,:)];
y_table = Y_table(p+1:end,:);

[y_vec,X_vec] = table2vec(y_table,X_table,model);

%% pack condition with initials;
condition.pphi = options.pphi_init;
condition.mmu = options.mmu_init;
condition.Ssgima = options.Ssigma_init;

%% generate mmu
% first you construct the transformed data as in the note
ystar = y;
Pphi = reshape(condition.pphi,[N N*p]);
for lag = 1:p
	ystar = ystar - condition.pphiPphi*X_table(:,2:end)';
end
for

% %% put data into the conditionals
% condition.y = y;
% condition.Z = Z;
% 
% % draw bbeta, Ssigma
% burnin = options.burnin;
% ndraws = options.ndraws; % post burnin draws
% count = 0;
% 
% % Burnin draws
% while count < burnin
% 	bbeta_draw = post_draw_bbeta_indie(model,condition,prior);
% 	condition.bbeta = bbeta_draw;
% 	Ssigma_draw = post_draw_Ssigma_indie(model,condition,prior);
% 	condition.bbeta = Ssigma_draw;
% 	count = count + 1;
% 	if mod(count,100) == 0
% 		str = sprintf('Drawed %d for burning.\b',count);
% 		disp(str);
% 	end
% 	if mod(count,666) == 0
% 		str = sprintf('Drawed %d for burning.\b',count);
% 		disp(str);
% 	end
% end
% 
% % Draws post burnin
% draws.bbeta = zeros(K,1,ndraws);
% draws.Ssigma = zeros(N,N,ndraws);
% for i_draw = 1:ndraws
% 	bbeta_draw = post_draw_bbeta_indie(model,condition,prior);
% 	condition.bbeta = bbeta_draw;
% 	draws.bbeta(:,i_draw) = bbeta_draw;
% 	
% 	Ssigma_draw = post_draw_Ssigma_indie(model,condition,prior);
% 	condition.bbeta = Ssigma_draw;
% 	draws.Ssigma(:,:,i_draw) = Ssigma_draw;
% 	
% 	%% show something I don't panic
% 	if mod(i_draw,100) == 0
% 		str = sprintf('Drawed %d for storage.\b',i_draw);
% 		disp(str);
% 	end
% end

end