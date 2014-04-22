%% test on the toy real dataset
clear; clc;
addpath(genpath('.'));

load('toy_N=102.mat');

%%
data.W = W;
data.F = F;
data.D = D;



hyper.K = 5; % topic/group number
N = size(F,1); % graph size
hyper.N = N;
V = 0.0;
for p = 1:hyper.N
    V_tmp = max(max(W{p}));
    if(find(W{p}==0))
        fprintf('find zeros\n')
    end
    if(V_tmp > V)
        V = V_tmp;
    end
end

hyper.V = double(V); % dictionary length
hyper.M = 1;

% threshold D
D = double(D>0);
data.D = mat2cell(D,ones(N,1),ones(N,1));

%%
[ LFN ] = Gibbs_sampling(data,  hyper );
save('toy_rslt.mat','LFN','hyper');