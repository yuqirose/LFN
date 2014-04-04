function [ prob ] = topic_posterior( p, Wpc, topic_cnt_p, group_cnt_p, params,hyper )
%TOPIC_POSTERIOR Summary of this function goes here
%   Detailed explanation goes here
K = hyper.K;

Theta= params.Theta;
Phi = params.Phi;
Beta = params.Beta;

prob = zeros(1,K);

for k = 1: K
    
   topic_cnt = topic_cnt_p;
   topic_cnt(k) = topic_cnt(k)+1;
   topic_prob= topic_cnt/sum(topic_cnt);
 
   group_cnt_p= group_cnt_p/sum(group_cnt_p);

   
   tmp = 1/(1+exp(topic_prob*Phi*group_prob'));
   prob(k) =  Theta(p,k)*Beta(k,Wpc) * tmp;  
end

prob(k)/sum(prob);
end

