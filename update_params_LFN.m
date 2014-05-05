function [ params_new ] = update_params_LFN (F,D, params, Tp_count, Tw_count, G_count, G )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    % 1. update Theta
    Theta = bsxfun(@rdivide, Tp_count, sum(Tp_count,2));
    [N,K] = size(Theta);
   

    % 2. update B (Right part)
    numer_B = zeros(K,K);
    denom_B = zeros(K,K);
    M = length(G{1,1});
    % number of time windows

    for p = 1:N
        for q = p+1:N
            for m = 1:M
                i = G{p,q}(m);
                j = G{q,p}(m);
                % {i, j} are topic id
                numer_B(i,j) = numer_B(i,j)+ 1*D{p,q}(m);
                denom_B (i,j) = denom_B(i,j) +1; 
                numer_B(j,i) = numer_B(j,i)+ 1*D{q,p}(m);
                denom_B (j,i) = denom_B(j,i) +1; 
            end
        end
    end
    denom_B = denom_B+1e-32;
    B = numer_B./ denom_B;
    if(isnan(sum(B(:))))
        %keyboard
    end
    % 3. update Beta (left part)
    % Tw_count: KxV matrix
    % normalize each row respectively
    Beta = bsxfun(@rdivide, Tw_count, sum(Tw_count,2));

    % 4. update Tau (bottom part)
    % no analytical solution
    % use numerical method
    Tp_weight = bsxfun(@rdivide, Tp_count, sum(Tp_count,2));
    G_weight = bsxfun(@rdivide, G_count, sum(G_count,3));

    feature = cell(N,N);
    for p=1:N
        for q=p+1:N
            feature{p,q} = [Tp_weight(p,:)*Tp_weight(q,:)'; squeeze(G_weight(p,q,:))'*squeeze(G_weight(q,p,:)); 1];
            feature{q,p} = feature{p,q};
        end
    end
    

    Tau = params.Tau;
    lambda = 0.001;
    [Tau, negL] = minimize(Tau, @calcFTGTau, 100, feature, F, lambda);
    Tau
     
    % 5. update Theta_prime
    for p = 1:N
        tmp = squeeze(G_count(p,:,:));
        Theta_prime(p,:) = sum(tmp)./ sum(sum(tmp));
    end
    
    % verify the bsxfun calculation -->
    verif = 0;
    if(verif)
        err = 0;
        for p=1:N
            for q=p+1:N
                tmp1 = squeeze(G_count(p,q,:)/sum( G_count(p,q,:)));
                err = err+sum(abs(tmp1-squeeze(G_weight(p,q,:))));
                tmp2 = squeeze(G_count(q,p,:)/sum( G_count(q,p,:)));
                err = err+sum(abs(tmp2-squeeze(G_weight(q,p,:))));
            end
        end
        if(err>1e-3)
           error('the bsxfun used incorrectly'); 
        end
    end
    % <-- verify the bsxfun calculation 

    params_new.Theta =Theta;
    params_new.Theta_prime = Theta_prime;
    params_new.Tau = Tau;
    params_new.Beta = Beta;
    params_new.B = B;


end

function [value, grad] = calcFTGTau(Tau, feature, F, lambda)
    N = size(F,1);
    K = length(Tau);
    
    % 4.1 do inference
    pred = zeros(N,N);
    for p=1:N
        for q=p+1:N
            pred(p,q) = 1/(1+exp(-Tau'*feature{p,q}));
%            pred(q,p) = pred(q,p);
        end
    end
    
    % 4.2 calculate the gradient
    grad = zeros(K,1);
    for p=1:N
        for q=p+1:N
            grad = grad + (F(p,q)-pred(p,q))*feature{p,q};
%            grad = grad + (F(q,p)-pred(q,p))*feature{q,p};
        end
    end
    
    % 4.3 calculate the value
    value = 0;
    for p=1:N
        for q=p+1:N
            value = value + F(p,q)*log(pred(p,q)+1e-32) + (1-F(p,q))*log(1-pred(p,q)+1e-32);
%            value = value + F(q,p)*log(pred(q,p)+1e-32) + (1-F(q,p))*log(1-pred(q,p)+1e-32);
        end
    end

    grad = -grad(:);
    value = -value;
    value = value + 0.5*lambda*Tau'*Tau;
    grad = grad + lambda*Tau;
end

