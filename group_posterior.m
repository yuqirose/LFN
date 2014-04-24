function [ prob ] = group_posterior(p, q, Dpqm, Dqpm, Gqpm, Fpq, T_feature, group_cnt_p, Gq_weight, params, hyper)
% this function infers the posterior probability distribution
% that the pair of user (p, q) at time slide m have dialog on a specific topic

% key components: 
% c1. Dpqm, Dqpm: both are included because change of Gpqm effects the value
%       of both of them
% c2. (Gqpm, T_feature --> F(p,q)) 
%       T_feature = Tp_weight*Tq_weight
%       is used to calculate the link probability F(p,q) and F(q,p)
% c3. Gq_weight: the average weight vector of Gqpm, used to calculate the grouping feature


K = hyper.K;
Theta_prime= params.Theta_prime;
B = params.B;
Tau = params.Tau;

prob = zeros(1,K);

for k = 1:K
        
    % 1. the middle part, followship from topics and groups information
    group_cnt = group_cnt_p;
    group_cnt(k) = group_cnt(k)+1;  %???
    Gp_weight= group_cnt/sum(group_cnt);

    % followship between user p and q
	F_prob= 1/(1+exp(-Tau(1)*T_feature - Tau(2)*Gq_weight*Gp_weight-Tau(3)));

    % 2. the right part, followship from param and to dialogue
    link_prob = (B(k,Gqpm)^Dpqm)*((1-B(k,Gqpm))^(1-Dpqm)); % bernoulli prob
    link_prob = link_prob*(B(Gqpm,k)^Dqpm)*((1-B(Gqpm,k))^(1-Dqpm));

    prob(k) = (Theta_prime(p,k))*link_prob*(F_prob^Fpq)*(1-F_prob)^(1-Fpq);
   
end

prob = prob /sum(prob); 

end

