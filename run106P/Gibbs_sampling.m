function [ T,G, params, LogLike_List] = Gibbs_sampling(data,  hyper,params, prefix)
% GIBBS_SAMPLING : Gibbs sampling method for LFN model
% Direct implementation of Gibb sampling
    thres = 1e-3;
    %% intialization
    % preparison
    W = data.W; 
    F = data.F;
    D = data.D;

    K = hyper.K; % topic/group number
    N = hyper.N; % graph size
    V = hyper.V; % dictionary length
    M = hyper.M;

    % randomly initialize the parameters
    Theta = ones(N, K)*(1/K);
    Theta_prime = ones(N, K)*(1/K);
    Beta = ones(K,V)*(1/V);
    Tau = rand(3,1);
    B = diag(ones(K,1)*0.8) +0.1*ones(K,K);
    
    params.Theta = Theta;
    params.Theta_prime = Theta_prime;
    params.Beta = Beta;
    params.Tau = Tau;
    params.B = B;


    % randomly sample the initial latent variables
    % T: NxN cell array
    % G: NxN cell array
    [ T, G] = Rnd_generateLFN(W, D, K,M);
    
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
    
    %% EM algorithm using Gibbs Sampling for Inference at E stage
    MaxIter = 20;
    MaxSubIter = 5;%10;%20;
    LogLike_List = [];
    [LogLike, LogLike1, LogLike2, LogLike3] = loglike_LFN(W,F,D,T,G,params,hyper);
    disp(['Initialization-> L: ' num2str(floor(LogLike)) ' L1: ' num2str(floor(LogLike1)) ' L2: ' num2str(floor(LogLike2)) ' L3: ' num2str(LogLike3)]);
    
    process.InitParam = params;
    process.InitL = [LogLike, LogLike1, LogLike2, LogLike3];
    Tpcount = cell(N,1);
    Tp_weights = cell(N,1); % assistant variables, used to calc P(F|T,G)
    Gp_weights = cell(N,1); % assistant variables, used to calc P(F|T,G)
    for p=1:N
        Tpcount{p} = Tp_count(p,:);
        Tp_weights{p} = Tpcount{p}'/sum(Tpcount{p});
        Gp_weights{p} = bsxfun(@rdivide, squeeze(G_count(p,:,:)), sum(squeeze(G_count(p,:,:)),2));
    end

    Gp_feature = zeros(N,N);
    for iter = 1:MaxIter
        % E stage
        for subiter = 1:MaxSubIter
            % Left part (T): Sampling Users' content
            
            % note 1. during sampling of Left part, the right part, G is
            % unchanged, use Gp_weight as input to the User-sampler
    
            % before sampling latent topics
            % calculate the Grouping feature 
            % which is used a lot but does not change during T sampling
            for p=1:N
                for q=p+1:N
                    Gp_feature(p,q) = Gp_weights{p}(q,:)*Gp_weights{q}(p,:)';
                    Gp_feature(q,p) = Gp_feature(p,q);
                end
            end

            % G_weight_Average
            G_weights = zeros(N,2);
            TpcountOther = cell(N,2);
            for p=1:N
                TpcountOther{p,1} = zeros(1,K);
                TpcountOther{p,2} = zeros(1,K);
            end
            F = F-diag(diag(F));
            for p = 1:N
                qc = find(F(p,:)>0);
                qd = find(F(p,:)==0);
                qd = qd(qd~=p);
                if(~isempty(qc))
                    G_weights(p,1) = mean(Gp_feature(p,qc));
                    TpcountOther{p,1} = sum(Tp_count(qc,:),1)/length(qc);
                    if(length(TpcountOther{p,1})<5)
                        keyboard
                    end
                end
                if(~isempty(qd))
                    G_weights(p,2) = mean(Gp_feature(p,qd));
                    TpcountOther{p,2} = sum(Tp_count(qc,:),1)/length(qd);
                    if(length(TpcountOther{p,2})<5)
                        keyboard
                    end
                end
            end


            parfor p=1:N
                [T{p}, Tpcount{p}] = ...
                    sampleUser(p, W{p}, Tpcount{p}, T{p}, TpcountOther, G_weights, F, params, hyper);
                if(rem(p,10)==0)
                    fprintf('%d ',p);
                end
            end
            fprintf('\n');

            % after sampling each user, update the Tp_weight
            for p=1:N
                Tp_count(p,:)=Tpcount{p};
                Tp_weights{p} = Tpcount{p}'/sum(Tpcount{p});
            end
            
            
            % before sampling grouping nodes G
            % calculate the topic feature 
            % which is used a lot but does not change during G sampling 
            Tp_features = zeros(N,N);
            for p=1:N
                for q=p+1:N
                    Tp_features(p,q) = Tp_weights{p}'*Tp_weights{q};
                    Tp_features(q,p) = Tp_features(p,q);
                end
            end
            
            % Right part (G) : Sampling User-User communication
            for p=1:N
                for q = p+1:N
                    for m = 1:M
                        % p => q
                        G_count(p,q,G{p,q}(m)) = G_count(p,q,G{p,q}(m))-1;

                        Gpq_count = squeeze(G_count(p,q,:));
                        prob_G = ...
                            group_posterior( p, q, D{p,q}(m), D{q,p}(m), G{q,p}(m), ...
                                F(p,q), Tp_features(p,q), Gpq_count, Gp_weights{q}(p,:), params,hyper);
                        Gpqm_new = find(mnrnd(1,prob_G)==1);
                        G{p,q}(m) = Gpqm_new;
                        G_count(p,q, Gpqm_new) = G_count(p,q, Gpqm_new)+1;
                    end
                    Gp_weights{p} = bsxfun(@rdivide, squeeze(G_count(p,:,:)), sum(squeeze(G_count(p,:,:)),2));
                
                    for m = 1:M
                        % q => p
                        G_count(q,p,G{q,p}(m)) = G_count(q,p,G{q,p}(m))-1;

                        Gqp_count = squeeze(G_count(q,p,:));
                        prob_G = ...
                            group_posterior( q, p, D{q,p}(m), D{p,q}(m), G{p,q}(m), ...
                                F(q,p), Tp_features(q,p), Gqp_count, Gp_weights{p}(q,:), params,hyper);
                        Gpqm_new = find(mnrnd(1,prob_G)==1);
                        G{q,p}(m) = Gpqm_new;
                        G_count(q,p, Gpqm_new) = G_count(q,p, Gpqm_new)+1;
                    end
                    Gp_weights{q} = bsxfun(@rdivide, squeeze(G_count(q,:,:)), sum(squeeze(G_count(q,:,:)),2));

                end
            end

            % Check convergence
            [LogLike_new, loglike_new1, loglike_new2, loglike_new3] = loglike_LFN(W,F,D,T,G,params,hyper);
            if(abs(LogLike_new-LogLike) < thres )
                break;
            end
            fprintf('SubIter # %d\n',subiter);
            disp(['SubIter # ' num2str(subiter) '-> L: ' num2str(floor(LogLike_new)) ' L1: ' num2str(floor(loglike_new1)) ' L2: ' num2str(floor(loglike_new2)) ' L3: ' num2str(loglike_new3)]);
            LogLike = LogLike_new;
        end


    
        % Update parameters
        Tw_count = Tw_count*0;
        for p = 1:N
            T_p = T{p};
            C = numel(T_p); % total words
            for c = 1:C
                Wpc = W{p}(c);
                Tw_count(T_p(c),Wpc) = Tw_count(T_p(c),Wpc) +1;
            end
        end
        Tw_count(:,1:4)
        params = update_params_LFN(F,D,params,Tp_count,Tw_count, G_count, G);
        fprintf('Iter # %d\n',iter );
        disp(num2str([LogLike_new loglike_new1 loglike_new2 loglike_new3]));
        LogLike_List = [LogLike_List , LogLike_new ];
        process.params{iter} = params;
        process.L(iter,1) = LogLike_new;
        process.L(iter,2) = loglike_new1;
        process.L(iter,3) = loglike_new2;
        process.L(iter,4) = loglike_new3;

        matname = sprintf('N30_iter%d.mat',iter);
        matname = [prefix matname];
        save(matname, 'params','hyper','process');
    end
    save([prefix 'LearnProcess.mat'],'process');
end



function [Tp, Tp_count] = sampleUser(p, Wp, Tp_count, Tp, Tp_weightOther, Gp_feature, F, params, hyper);
% function [Tp, Tp_count] = sampleUser(p, Wp, Tp_count, Tp, Gp_count, F, params, hyper)
% Wp: the content of user p
% Tp: the topic of user p
% Tp_count: the count of topics of user p
% Gp_count: 
% function [T, Tp_count] = sampleUser(W, p, Tp_count)
    wunic = unique(Wp);
    C = length(wunic);
    % C = numel(Wp); % total words
    for c = 1:C
        idx = find(Wp==wunic(c));
        idxc = idx(1);

        Tp_count(Tp(idxc)) = Tp_count(Tp(idxc)) - 1;
        prob_T = topic_posterior(p, Wp(idxc), Tp_count, Tp_weightOther, Gp_feature, F, params, hyper);
        Tp_count(Tp(idxc)) = Tp_count(Tp(idxc)) + 1;
        
        for i=1:length(idx)
            Tp_count(Tp(idx(i))) = Tp_count(Tp(idx(i))) - 1;
            Tpc_new = find(mnrnd(1, prob_T)==1);
            Tp(idx(i)) = Tpc_new;
            Tp_count(Tpc_new) = Tp_count(Tpc_new)+1;
        end
    end
end
