function [T,G, params, LogLike_List] = GibbsSamplerLFN( WS , DS , D, F, nTopic  )
%GIBBSSAMPLERLFN Summary of this function goes here
%   Wrapper for LFN Gibbs sampling methods: WP word probility
nDoc = max(DS);
nTerm = max(WS);

W = cell(1,nDoc);
for i = 1:length(DS)
    W{DS(i)} = [W{DS(i)} WS(i)];
end

D = num2cell(full(D));

data.W = W;
data.F =  F;
data.D  = D;

hyper.K = nTopic;
hyper.V = nTerm;
hyper.M = 1;
hyper.N= nDoc;

[ T,G, params, LogLike_List] = Gibbs_sampling(data,  hyper );
    
end

