function bbeta_draw = post_draw_bbeta_indie(model,condition,prior)
%% unpack things
T = model.T; 
prior_bbeta_mean = prior.bbeta_mean; % how confident about the fake SSR in prior
prior_bbeta_cov = prior.bbeta_cov; % the S in the notes
y = condition.y; % the dependent TN-by-1, T # of obs, N # of variables
Z = condition.Z; % the TN-by-K matrix, K = sum k_n
Ssigma = condition.Ssigma; % the drown Ssigma matrix from previous step

%% compute the mean and variance of post distribution
% see companion note
post_bbeta_cov = (prior.bbeta_cov\eye(length(prior_bbeta_mean)) + Z'*kron(eye(T),Ssigma^-1)*Z)\eye(length(prior_bbeta_mean));
post_bbeta_mean =  post_bbeta_cov*(prior_bbeta_cov\prior_bbeta_mean + ...
	               Z'*kron(eye(T),inv(Ssigma))*y);
			   
%% actually draw the bbeta
[L,check] = chol(post_bbeta_cov,'lower');
if  ~issymmetric(post_bbeta_cov)
	% warning('posteriors covariance not numerically symmetric');
	post_bbeta_cov = L*L';
end
if check > 0 
	warning('posteriors covariance not positive definite');
	post_bbeta_cov = zeros(size(post_bbeta_cov));
end
bbeta_draw = mvnrnd(post_bbeta_mean,post_bbeta_cov);
bbeta_draw = bbeta_draw';

end