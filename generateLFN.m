function [ T, G, F ] = generateLFN()
%GENERATELFN Summary of this function goes here
%   Detailed explanation goes here
% step 1. specify the model parameters: 
%   number of users, N
%   number of words, V
%   number of topics = number of groups, K
%   Theta: KxN 
%   beta: KxV
%   Phi: KxK
%   Bern: KxK
    N = 100;
    V = 500;
    K = 5;

    Theta = rand(K,N); % can consider using Dirichlet r.v. instead
    Theta = bsxfun(@rdivide, Theta, sum(Theta));
    Beta = rand(K,V);
    Beta = bsxfun(@rdivide, Beta, sum(Beta,2));
    Phi = randn(K,K);
    B = rand(K,K); % each element is one parameter of Bernoulli dist.

% setp 2. sample the 1) topics, 2) words, 3) group connection, 4) dialog and 5) following-relations
    
    % 2.1 sample Left part
    C = floor(rand(N,1)*1000+1);
    % C(n): number of words from user n, from 1 to 1000
    fprintf(2,'generating Left topic model part\n');
    for n=1:N
        % sample C(n) topics for user n
        T{n} = sampleCat(Theta(:,n), C(n),1);
        % sample C(n) words for user n
        for k=1:K
            idx = find(T{n}==k);
            if(~isempty(idx))
                W{n}(idx) = sampleCat(Beta(k,:)', length(idx), 1);
            end
        end
    end

    fprintf(2,'generating Right component model part\n');
    % 2.2 sample Right part
    M = 50;%floor(rand(N,1)*10+1);
    % M(n): number of dialog involving user n, from 1 to 10
    G = cell(N,N);
    D = cell(N,N);
    for p=1:N
        for q=1:N
            theta = 0.5*(Theta(:,p)+Theta(:,q));
            % sample M(p) diaglog (latent) variables
            G{p,q} = sampleCat(theta, M,1);
        end
    end
    
    for p=1:N
        for q=1:N
            D{p,q} = zeros(M,1);
            for m=1:M
                D{p,q}(m) = rand(1,1)<= B(G{p,q}(m), G{p,q}(m)); % bernoulli parameter, M(p)x1
            end
        end
    end

    % 2.3 sample the Following relation
    fprintf(2,'generating followship part\n');
    F = zeros(N,N);
    % average topic weight of one user
    mT = zeros(N,K);
    for n=1:N
        for k=1:K
            mT(n,k) = sum(T{n}==k);
        end
    end
    mT = bsxfun(@rdivide, mT, sum(mT,2));

    mG = cell(N,N);
    for p=1:N
        for q=1:N
            mG{p,q} = zeros(K,1);
            for k=1:K
                mG{p,q}(k) = sum(G{p,q}==k);
            end
            mG{p,q} = mG{p,q}/sum(mG{p,q});
        end
    end

    for p=1:N
        for q=1:N
            F(p,q) = rand(1,1)<1/(1+exp(mT(p,:)*Phi*mG{p,q}));
            %F(p,q) = mT(p,:)*Phi*mG{p,q};
        end
    end
%    keyboard
    % average grouping weight of one pair
    
    Theta = Theta';
%    save('fake_data.mat',);
end  


