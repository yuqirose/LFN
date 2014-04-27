function prob = topic_posterior(p, Wpc, topic_cnt_p, Tp_weights, Gp_feature, F, params, hyper)
% function [ prob ] = topic_posterior( p, Wpc, topic_cnt_p, group_cnt_p, F, params,hyper )
%TOPIC_POSTERIOR Summary of this function goes here
%   Detailed explanation goes here

% Tp_weights{p,1}: the avearge weight of other users with connection
% Tp_weights{p,2}: the avearge weight of other users without connection

    K = hyper.K;

    Theta= params.Theta;
    Tau = params.Tau;
    Beta = params.Beta;
    Beta = Beta+1e-32;
    Beta = bsxfun(@rdivide, Beta, sum(Beta,2));
    prob = zeros(1,K);

    qorder = randperm(hyper.N);
    F1 = sum(F(:));
    F0 = length(F(:))-F1;
    for k = 1: K

        topic_cnt = topic_cnt_p;
        topic_cnt(k) = topic_cnt(k)+1;
        topic_prob= topic_cnt/sum(topic_cnt);

        % when the topic of one use is updated, 
        % it should effect the followship relationship with all other users
        %
        F_prob = 0;
        Fcount = 0;
        nc = sum(F(p,:));
        nd = hyper.N-nc-1;
        if(nc>0)
            followprob = 1/(1+exp(-Tau(1)*topic_prob*Tp_weights{p,1}'-Tau(2)*Gp_feature(p,1)-Tau(3)));
            F_prob = F_prob + nc*log(followprob+1e-32);
        end
        followprob = 1/(1+exp(-Tau(1)*topic_prob*Tp_weights{p,2}'-Tau(2)*Gp_feature(p,2)-Tau(3)));
        F_prob = F_prob + nd*log(1-followprob+1e-32);
        
        prob(k) = log(Theta(p,k)) + log(Beta(k,Wpc)) + F_prob;% * tmp;  
    end

    prob = exp(prob-max(prob));
    prob = prob./sum(prob);
    if(isnan(sum(prob)))
        keyboard 
    end
end
