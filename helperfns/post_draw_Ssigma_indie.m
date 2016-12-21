function Ssigma_draw = post_draw_Ssigma_indie(model,condition,prior)
%% unpack things
T = model.T; 
N = model.N;
prior_SSR = prior.SSR; % how confident about the fake SSR in prior
prior_nnu = prior.nnu; % the S in the notes

%% unpack data and conditoned draws from previous steps
E = condition.y-condition.Z*condition.bbeta; % the observed? error terms
E = reshape(E,[N,T]);

%% compute the mean and variance of post distribution
% see companion note
post_SSR = prior_SSR + E*E';
post_nnu =  T + prior_nnu;
			   
%% actually draw the bbeta
[~,check] = chol(post_SSR,'lower');
if check ~=0
	check
end
Ssigma_draw = iwishrnd(post_SSR,post_nnu);
end