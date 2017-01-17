function draws = dummyVAR_Gibbs(y_table,regimes,model,priors,options)
%% Unpack stuff
N       = size(y_table,2);
model.N = N;
p       = model.p;
T       = size(y_table,1);
model.T       = T;
breakdate = model.breakdate;
condition.regimes = regimes;
condition.y_table = y_table;

%% derivative things/ some constructions from inputs
M1 = max(regimes(:,1)); % assumes all regimes are present
M2 = max(regimes(:,2)); % assumes all regimes are present
model.M1 = M1;
model.M2 = M2;
condition.Ssigma_array = zeros(N,N,M2);
for m1 = 1:M1
	condition.AR_dum_mat(:,m1) = regimes(:,1) == m1;
end
for m2 = 1:M2
	condition.cov_dum_mat(:,m2) = regimes(:,2) == m2;
end
pphi_init = zeros(M1*p*N^2,1);
condition.pphi = options.pphi_init;
condition.Ssigma_array = options.Ssigma_array_init;
condition.break_dum_mat = zeros(T,2); % first col pre, second post break 84
condition.break_dum_mat(1:breakdate,1) = ones(breakdate,1);
condition.break_dum_mat(breakdate+1:end,2) = ones(T-breakdate,1);

%% initializae containers
draws.pphi         = options.pphi_init;
draws.mmu          = options.mmu_init;
draws.Ssigma_array = options.Ssigma_array_init;

%%% main body of Gibbs sampler
% specification for transformed model for drawing mmu
mmu_model.T          = N*(T-p);
mmu_prior.bbeta_mean = priors.mmu_mean; % prior mean for intercepts
mmu_prior.bbeta_cov  = priors.mmu_cov; % prior covariance wrt intercepts
mmu_condition.Ssigma = 1; % homogenized by the transformation

% spec for drawing pphi
pphi_model.T = N*(T-p);
pphi_prior.bbeta_mean = priors.pphi_mean; % confident re fake SSR in prior
pphi_prior.bbeta_cov = priors.pphi_cov; % the S in the notes
pphi_condition.Ssigma = 1; % the drown Ssigma matrix from previous step

% spec for drawing Ssigma
Ssigma_prior.SSR = priors.Ssigma_SSR;
Ssigma_prior.nnu = priors.Ssigma_nnu;
Ssigma_condition.bbeta = zeros(N,1);
Ssigma_model.N = N; % because there's N shocks
Ssigma_array = zeros(N,N,M2);

count = 0;
i_draw = 1;
while i_draw <= options.R
	%% draw mmu
	S2star = zeros(T*N,2*M1*N);
	y2star = zeros(T*N,1);
	for t = 1:T
		S2star(1+N*(t-1):N+N*(t-1),:) = transform_St2star_for_mmu(t,model,condition);
		y2star(1+N*(t-1):N+N*(t-1),:) = transform_yt2star_for_mmu(t,model,condition);
	end
	mmu_condition.y = y2star(N*p+1:end,:); % the dependent TN-by-1, T # of obs, N # of variables
	mmu_condition.Z = S2star(N*p+1:end,:); % the TN-by-K matrix, K = sum k_n
	mmu = post_draw_bbeta_indie(mmu_model,mmu_condition,mmu_prior);
	condition.mmu = mmu; % key step of updating cond from draw
	if count >= options.burnin
		draws.mmu(:,i_draw) = mmu;
	end

	%% draw pphi
	[pphi_y2star,pphi_Z] = transform_for_pphi(model,condition);
	pphi_condition.y = pphi_y2star(N*p+1:end,:); % the dependent TN-by-1, T # of obs, N # of variables
	pphi_condition.Z = pphi_Z(N*p+1:end,:); % the TN-by-K matrix, K = sum k_n
	stable = false;
	stable_vec = false(M1,1);
	stable_count = 0;
	while stable ~= 1 % only accept draws that are stable
		pphi = post_draw_bbeta_indie(pphi_model,pphi_condition,pphi_prior);
		pphi_array = reshape(pphi,N*N*p,M1);
		for m1 = 1:M1
			stable_vec(m1) = checkpphi_stable(pphi_array(:,m1),model.N,model.p);
		end
		stable = min(stable_vec);
		stable_count = stable_count + 1;
		stable = true;%!!!!!!!!!!!!!!! bypass stability check for once
	end
	disp(strcat('stable found after:| ',num2str(stable_count),'| draws'));

	condition.pphi = pphi; % key step of updating cond from draw

	if count >= options.burnin
		draws.pphi(:,i_draw) = pphi;
	end

	%% draw Ssigma_array
	[E_table_forSsigma] = transform_for_Ssigma(model,condition);
	E_table_forSsigma = E_table_forSsigma(p+1:end,:);
	regimes_forSsigma = regimes(p+1:end,2);
	for m2 = 1:M2
		Em = E_table_forSsigma(regimes_forSsigma==m2,:);
		Ssigma_model.T = size(Em,1);
		temp = Em';
		Ssigma_condition.y = temp(:);
		Ssigma_condition.Z = ones(Ssigma_model.T*Ssigma_model.N,Ssigma_model.N);
		% pack things in
		Ssigma_array(:,:,m2) = post_draw_Ssigma_indie(Ssigma_model,Ssigma_condition,Ssigma_prior);
	end
	condition.Ssigma_array = Ssigma_array; % key step of updating cond from draw
	if count >= options.burnin
		draws.Ssigma_array(:,:,:,i_draw) = Ssigma_array;
	end

	%% Display for the burnin part
	if count < options.burnin
		if mod(count,100) == 0
			fprintf('Drawed %d for burning.\n',count);
		end
		if mod(count,666) == 0
			fprintf('Drawed %d for burning.\n',count);
		end
		count = count + 1;
	else
		%% Display
		if mod(i_draw,100) == 0
			fprintf('Drawed %d for storage.\n',i_draw);
		end
		i_draw = i_draw + 1;
	end

end

end
