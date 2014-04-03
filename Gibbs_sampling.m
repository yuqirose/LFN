% function [ LFN ] = Gibbs_sampling(data,  hyper )
%GIBBS_SAMPLING Summary of this function goes here
%   Detailed explanation goes her

load('fake_data.mat');
thres = 1e-3;
%% intialization
K = hyper.K; % topic/group number
N = hyper.N; % graph size
V = hyper.V; % dictionary length

W = data.W; % word counts
F = data.F;
D = data.D;


%% params
Theta = ones(N, K)*(1/K);
Beta = ones(K,V)*(1/V);
Phi =  diagonal(ones(K,1));
B = diagonal(ones(K,1));

%% latent variables
T = cell(N,1);
G = cell(N,N);

%% topic/link count
T_count = zeros(N,K);
G_count = zeros(N,K);

%% randomly assign labels
[T,G] = generateLFN(Theta,Beta, Phi, B);


MaxIter = 100;
MaxSubIter = 50;
LogLike = loglike_LFN ();
for iter = 1:MaxIter
    for subiter = 1:MaxSubIter
        
    % Sample T 
        for p = 1:N
            C = numel(W{p}); % total words
            for c = 1:C
                 Tnc_old = T{p}(c);
                 T_count(n,Tnc_old) = T_count(n,Tnc_old) -1;
                 % TBD: process with bows
                 prob_T = topic_posterior(p, Wpc, topic_cnt_p, group_cnt_p, params,hyper );
                 Tnc_new = find(mnrnd(prob_T,1)==1);
                 T{p}(c) = Tnc_new;
                 T_count(n,Tnc_new) = T_count(n,Tnc_new)+1;

            end

            % Sample G
            for q = 1:N
                for m = 1:M
                     Gnm_old = G{p}(m);
                     G_count(n,Gnm_old) = G_count(n,Gnm_old) -1;
                     % TBD:process with networks
                     prob_G = group_posterior( p,q, Dpqm, Gqpm, topic_cnt_p,group_cnt_p, params,hyper);
                     Gnm_new = find(mnrnd(prob_G,1)==1);
                     G{p}(m) = Gnm_new;
                     G_count(n,Gnm_new) = G_count(n,Gnm_new)+1;

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


