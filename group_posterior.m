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
Theta_prime(p,:) = Theta_prime(p,:)+1e-32;
Theta_prime(p,:) = Theta_prime(p,:)/sum(Theta_prime(p,:));
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
    link_prob = Dpqm*log(B(k,Gqpm)+1e-32)+ (1-Dpqm)*log(1-B(k,Gqpm)+1e-32); % bernoulli prob
    link_prob = link_prob + Dqpm*log(B(Gqpm,k)+1e-32)+ (1-Dqpm)*log(1-B(Gqpm,k)+1e-32); % bernoulli prob
    %link_prob = link_prob*(B(Gqpm,k)^Dqpm)*((1-B(Gqpm,k))^(1-Dqpm));

    prob(k) = log(Theta_prime(p,k)+1e-32) + link_prob + Fpq*log(F_prob+1e-32) + (1-Fpq)*log(1-F_prob+1e-32);
%    prob(k) = log(Theta_prime(p,k)) + link_prob ;
end

prob = exp(prob-max(prob));
prob = prob/sum(prob); 

    if(isnan(sum(prob)))
        %keyboard
    end
end

