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
B = diag(ones(K,1)*0.8) +0.1*ones(K,K);

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
[LogLike, LogLike1, LogLkie2, LogLike3] = loglike_LFN3 (W,F,D,T,G,params,hyper);
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
                 % G_count: NxNxK tensor
                 Gp_count = squeeze(G_count(p,:,:));
                 prob_T = topic_posterior(p, Wpc, T_count(p,:), Gp_count, F, params, hyper);
                 Tpc_new = find(mnrnd(1, prob_T)==1);
                 T{p}(c) = Tpc_new;
                 T_count(p,Tpc_new) = T_count(p,Tpc_new)+1;
            end
            if(rem(p,5)==0)
                fprintf('%d ',p);
            end
        end
        
        for p=1:N
            % Sample G
            for q = p+1:N
                for m = 1:M
                    Dpqm = D{p,q}(m);
                    Dqpm = D{q,p}(m);
                    
                    % p => q
                    Gpqm_old = G{p,q}(m);
                    Gqpm_old = G{q,p}(m);
                    G_count(p,q,Gpqm_old) = G_count(p,q,Gpqm_old) -1;
                     
                    % TBD:process with networks
                    Gpq_count = squeeze(G_count(p,q,:));
                    
                    prob_G = group_posterior( p,q, Dpqm, Dqpm, Gqpm_old, F(p,q), T_count(p,:), Gpq_count, params,hyper);
%                     [a,b] = max(prob_G);
%                     Gpqm_new = b;
                    Gpqm_new = find(mnrnd(1,prob_G)==1);
                    G{p,q}(m) = Gpqm_new;
                    G_count(p,q, Gpqm_new) = G_count(p,q, Gpqm_new)+1;

                    % q => p
                    Gqpm_old = G{q,p}(m);
                    Gpqm_old = G{p,q}(m);
                    G_count(q,p,Gqpm_old) = G_count(q,p,Gqpm_old) -1;
                     
                    % TBD:process with networks
                    
                    Gqp_count = squeeze(G_count(q,p,:));
                     
                    prob_G = group_posterior( q,p, Dqpm, Dpqm, Gpqm_old, F(q,p), T_count(q,:), Gqp_count, params,hyper);
%                     [a,b] = max(prob_G);
%                     Gqpm_new = b;
                    Gqpm_new = find(mnrnd(1,prob_G)==1);
                    G{q,p}(m) = Gqpm_new;
                    G_count(q,p, Gqpm_new) = G_count(q, p, Gqpm_new)+1;
                    
                end
            end

            if(rem(p,5)==0)
                fprintf('%d ',p);
            end
        end
        
                    
        % Check convergence
        [LogLike_new, loglike_new1, loglike_new2, loglike_new3] = loglike_LFN(W,F,D,T,G,params,hyper);
        if(abs(LogLike_new-LogLike) < thres )
           break;
        end
        disp(num2str([LogLike_new loglike_new1 loglike_new2 loglike_new3]));
    end

    

    % Update parameters
      params = update_params_LFN(params,T_count, G_count);
end

% end


