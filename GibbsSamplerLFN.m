function [T,G, params, LogLike_List] = GibbsSamplerLFN( WS , DS ,F, nTopic , BURNIN  )
%GIBBSSAMPLERLFN Summary of this function goes here
%   Wrapper for LFN Gibbs sampling methods: WP word probility
nDoc = max(DS);
nTerm = max(WS);

W = cell(nDoc,1);
for i = 1:length(DS)
    W{DS(i)} = [W{DS(i)} WS(i)];
end

data.F =  F;

data.D  = D;

hyper.K = nTopic;
hyper.V = nTerm;
hyper.M = 1;
[ T,G, params, LogLike_List] = Gibbs_sampling(data,  hyper );
    
end

