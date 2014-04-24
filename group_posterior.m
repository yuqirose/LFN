function [ prob ] = group_posterior(p,q, Dpqm, Dqpm, Gqpm, Fpq, topic_cnt_p,group_cnt_p, params,hyper)
%GROUP_POSTERIOR Summary of this function goes here
%   Detailed explanation goes here

K = hyper.K;
Theta_prime= params.Theta_prime;
Phi = params.Phi;
B = params.B;

prob = zeros(1,K);

for k = 1:K
        
   % 1. the middle part, followship from topics and groups information
    group_cnt = group_cnt_p;
    group_cnt(k) = group_cnt(k)+1;  %???
    group_prob= group_cnt/sum(group_cnt);
    topic_prob= topic_cnt_p/sum(topic_cnt_p);

   % followship between user p and q
	F_prob= 1/(1+exp(topic_prob*Phi*group_prob));

   % 2. the right part, followship from param and to dialogue
%    link_prob  = B(k,Gqpm)^Dpqm+(1-B(k,Gqpm))^(1-Dpqm);
    link_prob = (B(k,Gqpm)^Dpqm)*((1-B(k,Gqpm))^(1-Dpqm)); % bernoulli prob
    link_prob = link_prob*(B(Gqpm,k)^Dqpm)*((1-B(Gqpm,k))^(1-Dqpm));
   % ???
   % should we also consider the effect from B(q, Gpqm) where q ~= p?
   % ???

    prob(k) = (Theta_prime(p,k))*link_prob*(F_prob^Fpq)*(1-F_prob)^(1-Fpq);
   
end

prob = prob /sum(prob); 

end

