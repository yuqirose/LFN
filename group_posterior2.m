function [ prob ] = group_posterior( p,q, Dpqm, Gqpm, topic_cnt_p,group_cnt_p, params,hyper)
%GROUP_POSTERIOR Summary of this function goes here
%   Detailed explanation goes here

K = hyper.K;
Theta= params.Theta;
Phi = params.Phi;
B = params.B;

prob = zeros(1,K);

for k = 1:K
        
   % 1. the middle part, followship from topics and groups information
   group_cnt = group_cnt_p;
   group_cnt(k) = group_cnt(k)+1;  %???
   group_prob= group_cnt/sum(group_cnt);
   topic_prob= topic_cnt_p/sum(topic_cnt_p);
   tmp = 1/(1+exp(topic_prob*Phi*group_prob));

   % 2. the right part, followship from param and to dialogue
%    link_prob  = B(k,Gqpm)^Dpqm+(1-B(k,Gqpm))^(1-Dpqm);
   link_prob  = (B(k,Gqpm)^Dpqm)*((1-B(k,Gqpm))^(1-Dpqm)); % bernoulli prob
   prob(k) = (Theta(p,k)+ Theta(q,k))/2*link_prob*tmp;
   
end

prob = prob /sum(prob); 

end
