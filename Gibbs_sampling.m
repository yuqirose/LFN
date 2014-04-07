% function [ LFN ] = Gibbs_sampling(data,  hyper )
%GIBBS_SAMPLING Summary of this function goes here
%   Detailed explanation goes her
clear; clc;

load('fake_data.mat');

thres = 1e-3;
%% intialization
K = hyper.K; % topic/group number
N = hyper.N; % graph size
V = hyper.V; % dictionary length
M = hyper.M;

W = data.W; % word counts
F = data.F;
D = data.D;


params_org = params;

%% params
Theta = ones(N, K)*(1/K);
Beta = ones(K,V)*(1/V);
Phi =  diag(ones(K,1));
B = diag(ones(K,1)*0.9) +0.1*ones(K,K);

params.Theta = Theta;
params.Beta = Beta;
params.Phi = Phi;
params.B = B;

%% latent variables
% T = cell(1,N);
% G = cell(N,N);
% T = data.T;
% G = data.G;



%% randomly assign labels

[ T, G] = Rnd_generateLFN(W, D, K,M);
% [T,G] = generateLFN(Theta,Beta, Phi, B);

%% topic/link count
       

T_count = zeros(N,K);
G_count = zeros(N,N,K);

for p = 1:N
    T_p = T{p};
    C = numel(T_p); % total words
    for c = 1:C
        T_count(p,T_p(c)) = T_count(p,T_p(c))+1;
    end
    for q = 1:N
        G_pq = G{p,q};
        for m = 1:M
            G_count (p,q,G_pq(m)) =G_count (p,q,G_pq(m))+1;
        end
    end
end
%%

MaxIter = 100;
MaxSubIter = 50;
LogLike = loglike_LFN (W,F,D,T,G,params,hyper);
for iter = 1:MaxIter
    for subiter = 1:MaxSubIter
        % Sample T 
        for p = 1:N
            C = numel(W{p}); % total words
            for c = 1:C
                 Tpc_old = T{p}(c);
                 T_count(p,Tpc_old) = T_count(p,Tpc_old) -1;
                 % TBD: process with bows
                 Wpc = W{p}(c);
                 Gp_count = squeeze(sum(G_count(p,:,:),2));
                 prob_T = topic_posterior(p, Wpc, T_count(p,:), Gp_count, params,hyper );
                 Tpc_new = find(mnrnd(1, prob_T)==1);
                 T{p}(c) = Tpc_new;
                 T_count(p,Tpc_new) = T_count(p,Tpc_new)+1;

            end
            
            % Sample G
            for q = 1:N
                for m = 1:M
                     Gpqm_old = G{p,q}(m);
                     
                     G_count(p,q,Gpqm_old) = G_count(p,q,Gpqm_old) -1;
                     
                     % TBD:process with networks
                     Dpqm = D{p,q}(m);
                     Gpq_count = squeeze(G_count(p,q,:));
                     
                     prob_G = group_posterior( p,q, Dpqm, Gpqm_old, T_count(p,:),Gpq_count, params,hyper);
                     Gpqm_new = find(mnrnd(1,prob_G)==1);
                     G{p,q}(m) = Gpqm_new;
                     G_count(p,q, Gpqm_new) = G_count(p,q, Gpqm_new)+1;

                end
            end

        fprintf('%d ',p);
        end
        
                    
        % Check convergence
        LogLike_new = loglike_LFN(W,F,D,T,G,params,hyper);
        if(abs(LogLike_new-LogLike) < thres )
           break;
        end
        disp(LogLike_new);
    end
    
    % Update parameters
      params = update_params_LFN(params,T_count, G_count);
    
end

% end


