function [sim_val,sim_ind] = markovsim(states,transition,init_ind,periods)
% inputs: 
%   states         | n-by-1 vector of possible states
%   transition     | n-by-n Markov matrix
%   init_state_ind | index (1,...,n) 
%   periods        | length of simulation
% outpus:
%   sim_val        | simulated 1-by-periods vector of state values
%   sim_ind        | simulated 1-by-periods vector of state indexes
T = periods;
n = length(states);
sim_ind = zeros(1,T);
sim_val = zeros(1,T);
u = rand(1,periods);
CDF_mat = cumsum(transition,2);

sim_ind(1,1) = init_ind;
sim_val(1,1) = states(init_ind);
for t = 2:T
	j = find(u(t)<=CDF_mat(sim_ind(t-1),:),1,'first');
	sim_val(1,t) = states(j);
	sim_ind(1,t) = j;
end

end