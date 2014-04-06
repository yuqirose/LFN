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


%% params
% Theta = ones(N, K)*(1/K);
% Beta = ones(K,V)*(1/V);
% Phi =  diag(ones(K,1));
% B = diag(ones(K,1));

%% latent variables
T = cell(1,N);
G = cell(N,N);

%% topic/link count
T_count = zeros(N,K);
G_count = zeros(N,K);

%% randomly assign labels
% [T,G] = generateLFN(Theta,Beta, Phi, B);


MaxIter = 100;
MaxSubIter = 50;
LogLike = loglike_LFN ();
for iter = 1:MaxIter
    for subiter = 1:MaxSubIter
        
    % Sample T 
        for p = 1:N
            C = numel(W{p}); % total words
            for c = 1:C
                 Tpc_old = T{p}(c);
                 T_count(n,Tpc_old) = T_count(n,Tpc_old) -1;
                 % TBD: process with bows
                 Wpc = W{p}(c);
                 prob_T = topic_posterior(p, Wpc, topic_cnt_p, group_cnt_p, params,hyper );
                 Tpc_new = find(mnrnd(prob_T,1)==1);
                 T{p}(c) = Tpc_new;
                 T_count(n,Tpc_new) = T_count(n,Tpc_new)+1;

            end

            % Sample G
            for q = 1:N
                for m = 1:M
                     Gpm_old = G{p}(m);
                     G_count(n,Gpm_old) = G_count(n,Gpm_old) -1;
                     % TBD:process with networks
                     prob_G = group_posterior( p,q, Dpqm, Gqpm, topic_cnt_p,group_cnt_p, params,hyper);
                     Gpm_new = find(mnrnd(prob_G,1)==1);
                     G{p}(m) = Gpm_new;
                     G_count(n,Gpm_new) = G_count(n,Gpm_new)+1;

                end
            end
            
            % Check convergence
             LogLike_new = loglike_LFN();
             if(abs(LogLike_new-LogLike) < thres )
                 break;
             end

        end
           
    
    end
    
    
    
end

% end


