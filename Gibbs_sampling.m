function [ T,G, params, LogLike_List] = Gibbs_sampling(data,  hyper )
% GIBBS_SAMPLING : Gibbs sampling method for LFN model
    % Direct implementation of Gibb sampling
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

    Tp_count = zeros(N,K);
    Tw_count = zeros(K,V);
    G_count = zeros(N,N,K);

    for p = 1:N
        T_p = T{p};
        C = numel(T_p); % total words
        for c = 1:C
            Tp_count(p,T_p(c)) = Tp_count(p,T_p(c))+1;
            Wpc = W{p}(c);
            Tw_count(T_p(c),Wpc) = Tw_count(T_p(c),Wpc) +1;
        end
        for q = 1:N
            G_pq = G{p,q};
            for m = 1:M
                G_count (p,q,G_pq(m)) =G_count (p,q,G_pq(m))+1;
            end
        end
    end
    %%

    MaxIter = 5;
    MaxSubIter = 20;
    LogLike_List = [];
    [LogLike, LogLike1, LogLkie2, LogLike3] = loglike_LFN(W,F,D,T,G,params,hyper);
    
    for iter = 1:MaxIter
        for subiter = 1:MaxSubIter
            % Sample T             
            Tpcount = cell(N,1);
            Gpcount = cell(N,1);
            for p=1:N
                Tpcount{p} = Tp_count(p,:);
                Gpcount{p} = squeeze(G_count(p,:,:));
            end
            parfor p=1:N

                %[T{p}, Tp_count(p,:)] = sampleUser(p, W{p}, Tp_count(p,:), T{p}, squeeze(G_count(p,:,:)), F);
                [T{p}, Tpcount{p}] = sampleUser(p, W{p}, Tpcount{p}, T{p}, Gpcount{p}, F, params, hyper);
%                 fprintf('%d ',p);
            end
            for p=1:N
                Tp_count(p,:)=Tpcount{p};
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

                        prob_G = group_posterior( p,q, Dpqm, Dqpm, Gqpm_old, F(p,q), Tp_count(p,:), Gpq_count, params,hyper);
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

                        prob_G = group_posterior( q,p, Dqpm, Dpqm, Gpqm_old, F(q,p), Tp_count(q,:), Gqp_count, params,hyper);
    %                     [a,b] = max(prob_G);
    %                     Gqpm_new = b;
                        Gqpm_new = find(mnrnd(1,prob_G)==1);
                        G{q,p}(m) = Gqpm_new;
                        G_count(q,p, Gqpm_new) = G_count(q, p, Gqpm_new)+1;

                    end
                end
% 
%                 if(rem(p,5)==0)
%                     fprintf('%d ',p);
%                 end
            end


        % Check convergence
        [LogLike_new, loglike_new1, loglike_new2, loglike_new3] = loglike_LFN(W,F,D,T,G,params,hyper);
        if(abs(LogLike_new-LogLike) < thres )
           break;
        end
%         disp(num2str([LogLike_new loglike_new1 loglike_new2 loglike_new3]));
        fprintf('SubIter # %d\n',subiter);
        end



    % Update parameters
    
      for p = 1:N
        T_p = T{p};
        C = numel(T_p); % total words
        for c = 1:C
            Wpc = W{p}(c);
            Tw_count(T_p(c),Wpc) = Tw_count(T_p(c),Wpc) +1;
        end
      end
      params = update_params_LFN(F,D,params,Tp_count,Tw_count, G_count, G);
      fprintf('Iter # %d\n',iter );
      disp(num2str([LogLike_new loglike_new1 loglike_new2 loglike_new3]));
      LogLike_List = [LogLike_List , LogLike_new ];

    end
end




function [Tp, Tp_count] = sampleUser(p, Wp, Tp_count, Tp, Gp_count, F, params, hyper)
% Wp: the content of user p
% Tp: the topic of user p
% Tp_count: the count of topics of user p
% Gp_count: 
% function [T, Tp_count] = sampleUser(W, p, Tp_count)

    C = numel(Wp); % total words
    for c = 1:C
         Tp_count(Tp(c)) = Tp_count(Tp(c)) -1;

         % G_count: NxNxK tensor
%          Gp_count = squeeze(G_count(p,:,:));
         
         prob_T = topic_posterior(p, Wp(c), Tp_count, Gp_count, F, params, hyper);
         Tpc_new = find(mnrnd(1, prob_T)==1);
         Tp(c) = Tpc_new;
         Tp_count(Tpc_new) = Tp_count(Tpc_new)+1;
    end

end
