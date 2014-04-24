function prob = topic_posterior(p, Wpc, topic_cnt_p, Tp_weights, Gp_feature, F, params, hyper)
% function [ prob ] = topic_posterior( p, Wpc, topic_cnt_p, group_cnt_p, F, params,hyper )
%TOPIC_POSTERIOR Summary of this function goes here
%   Detailed explanation goes here

    K = hyper.K;

    Theta= params.Theta;
    Phi = params.Phi;
    Beta = params.Beta;

    prob = zeros(1,K);

    group_prob= bsxfun(@rdivide, group_cnt_p, sum(group_cnt_p,2));

    for k = 1: K

        topic_cnt = topic_cnt_p;
        topic_cnt(k) = topic_cnt(k)+1;
        topic_prob= topic_cnt/sum(topic_cnt);

        % when the topic of one use is updated, 
        % it should effect the followship relationship with all other users
        %
        F_prob = 0;
        for q = 1:size(group_prob,1)
            if(q~=p)
                followprob = 1/(1+exp(-Tau(1)*topic_prob*Tp_weights{q}-Tau(2)*Gp_feature(p,q)-Tau(3)));
                F_prob = F_prob + F(p,q)*log(followprob) + (1-F(p,q))*log(1-followprob);
            end
        end

        prob(k) = log(Theta(p,k)) + log(Beta(k,Wpc)) + F_prob;% * tmp;  
    end

    prob = exp(prob-max(prob));
    prob = prob./sum(prob);
end
