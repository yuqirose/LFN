function [ T,G, params, LogLike_List] = Gibbs_sampling(data,  hyper, params, nameprefix )
% GIBBS_SAMPLING : Gibbs sampling method for LFN model
% Direct implementation of Gibb sampling
    
    

    thres = 5e-5;
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
    %{
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
    %}
    Theta = params.Theta;
    Theta_prime = params.Theta_prime;
    Beta = params.Beta;
    Tau = params.Tau;
    B = params.B;

    % randomly sample the initial latent variables
    % T: NxN cell array
    % G: NxN cell array
    [ T, G] = Rnd_generateLFN(W, D, K,M);
    maxC = 0;
    for i=1:N
        if(maxC<length(W{i}))
            maxC = length(W{i});
        end
    end
    Ws = zeros(maxC, N);
    for i=1:N
        Ws(1:length(W{i}), i) = W{i}';
    end
    Ts = Ws*0;
    for i=1:N
        Ts(1:length(T{i}),i) = T{i};
    end
    numToken = sum(Ws>0);

    Tp_count = DataCount(Ts,N,K); % K x N matrix
    Tw_count = zeros(K,V);
    Tw_count = WordTopicCount(Ws, Ts, Tw_count);
    
    G_count = zeros(N,N,K);
    
    for p = 1:N
        for q = 1:N
            G_pq = G{p,q};
            for m = 1:M
                G_count (p,q,G_pq(m)) =G_count (p,q,G_pq(m))+1;
            end
        end
    end
    
    %% EM algorithm using Gibbs Sampling for Inference at E stage
    MaxIter = 25;
    histLogLike = cell(MaxIter,1);

    MaxSubIter = 10;%20;
    LogLike_List = [];
    [LogLike, LogLike1, LogLike2, LogLike3] = loglike_LFN(W,F,D,Ts,G,params,hyper);
    disp(['Initialization-> L: ' num2str(floor(LogLike)) ' L1: ' num2str(floor(LogLike1)) ' L2: ' num2str(floor(LogLike2)) ' L3: ' num2str(LogLike3)]);
    
    process.InitParam = params;
    process.InitL = [LogLike, LogLike1, LogLike2, LogLike3];
    Gp_weights = cell(N,1); % assistant variables, used to calc P(F|T,G)
    for p=1:N
        Gp_weights{p} = bsxfun(@rdivide, squeeze(G_count(p,:,:)), sum(squeeze(G_count(p,:,:)),2));
    end

    Gp_feature = zeros(N,N);
    for iter = 1:MaxIter
        histLogLike{iter} = [];

        % E stage
        for subiter = 1:MaxSubIter
            % Left part (T):  Topic Model
            for p=1:N
                for q=p+1:N
                    Gp_feature(p,q) = Gp_weights{p}(q,:)*Gp_weights{q}(p,:)';
                    Gp_feature(q,p) = Gp_feature(p,q);
                end
            end
            [tp, tpcount] = GibbsTopicPosterior(Ws, Ts, numToken, Tp_count, Gp_feature, F, Tau, Theta', Beta);
            % attention: the input Theta should be K x N matrix

            tp = tp+1; % from 0~(K-1) in c++ to 1~K in MATLAB
            Ts = reshape(tp, [size(Ts,1), size(Ts,2)]);
            Tp_count = reshape(tpcount, [size(Tp_count,1), size(Tp_count,2)]);
            %fprintf('change in Ws: %d\n', sum(abs(sum(Ws>0)-numToken)));
            %fprintf('change in Ts: %d\n', sum(abs(sum(Ts>0)-numToken)));
            %fprintf('change in Tp: %d\n', sum(abs(sum(Tp_count)-numToken)));
            for p=1:N
                for q=p+1:N
                    Tp_features(p,q) = Tp_count(:,p)'*Tp_count(:,q)/numToken(p)/numToken(q);
                    Tp_features(q,p) = Tp_features(p,q);
                end
            end
            % Right part (G) : Group Model
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
            [LogLike_new, loglike_new1, loglike_new2, loglike_new3] = loglike_LFN(W,F,D,Ts,G,params,hyper);
            histLogLike{iter} = [histLogLike{iter}; LogLike_new, loglike_new1, loglike_new2, loglike_new3];
            fprintf('SubIter # %d\n',subiter);
            disp(['SubIter # ' num2str(subiter) '-> L: ' num2str(floor(LogLike_new)) ' L1: ' num2str(floor(loglike_new1)) ' L2: ' num2str(floor(loglike_new2)) ' L3: ' num2str(loglike_new3)]);
            if(abs(LogLike_new-LogLike)/abs(LogLike_new) < thres )
                break;
            end
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
        params = update_params_LFN(F,D,params,Tp_count',Tw_count, G_count, G);
        Theta = params.Theta;
        Beta = params.Beta;
        % previously the Tp_count is a NxK matrix 
        fprintf('Iter # %d\n',iter );
        disp(num2str([LogLike_new loglike_new1 loglike_new2 loglike_new3]));
        LogLike_List = [LogLike_List , LogLike_new ];
        process.params{iter} = params;
        process.L(iter,1) = LogLike_new;
        process.L(iter,2) = loglike_new1;
        process.L(iter,3) = loglike_new2;
        process.L(iter,4) = loglike_new3;

        matname = [nameprefix num2str(iter) '.mat'];
%         matname = sprintf('N30_iter%d.mat',iter);
        save(matname, 'params','hyper','process', 'histLogLike');
    end
    save('LearnProcess.mat','process');
end



function [Tp, Tp_count, time1, time2, time3] = sampleUser(p, Wp, Tp_count, Tp, Tp_weights, Gp_feature, F, params, hyper);
% function [Tp, Tp_count] = sampleUser(p, Wp, Tp_count, Tp, Gp_count, F, params, hyper)
% Wp: the content of user p
% Tp: the topic of user p
% Tp_count: the count of topics of user p
% Gp_count: 
% function [T, Tp_count] = sampleUser(W, p, Tp_count)
    wunic = unique(Wp);
    C = length(wunic);
    time1 = 0;
    time2 = 0;
    time3 = 0;
    % C = numel(Wp); % total words
    for c = 1:C
        idx = find(Wp==wunic(c));
        idxc = idx(1);

        Tp_count(Tp(idxc)) = Tp_count(Tp(idxc)) - 1;
        [prob_T, time1unit, time2unit] = topic_posterior(p, Wp(idxc), Tp_count, Tp_weights, Gp_feature, F, params, hyper);
        Tp_count(Tp(idxc)) = Tp_count(Tp(idxc)) + 1;
        time1 = time1+time1unit;
        time2 = time2+time2unit;

        start = tic;
        for i=1:length(idx)
            Tp_count(Tp(idx(i))) = Tp_count(Tp(idx(i))) - 1;
            Tpc_new = find(mnrnd(1, prob_T)==1);
            Tp(idx(i)) = Tpc_new;
            Tp_count(Tpc_new) = Tp_count(Tpc_new)+1;
        end
        time3 = time3+toc(start);
    end
end

function Count = DataCount(data, N, C)
    % data is N column matrix
    Count = zeros(C,N);
    for c=1:C
        Count(c,:) = sum(data==c);
    end
end

function Count = WordTopicCount(Ws, Ts, Count)
    % Count: K X V matrix
    Count = Count*0;
    Ws = Ws(:);
    Ts = Ts(:);
    Ts = Ts(Ws>0);
    Ws = Ws(Ws>0);
    for i=1:length(Ts)
        Count(Ts(i), Ws(i)) = Count(Ts(i), Ws(i))+1;
    end
end
