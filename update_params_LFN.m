function [ params_new ] = update_params_LFN (W,F,D, params, Tp_count, Tw_count, G_count, G )
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
    B = numer_B./ denom_B;

    % 3. update Beta (left part)
    % Tw_count: KxV matrix
    % normalize each row respectively
    Beta = bsxfun(@rdivide, Tw_count, sum(Tw_count,2));

    % 4. update Phi (bottom part)
    % Tp_count: NxK matrix
    % Tp_count(i,j): for user i, the number of words in topic j
    numer_Phi = zeros(K,K);
    denom_Phi = zeros(K,K);
    Tp_weight = bsxfun(@rdivide, Tp_count, sum(Tp_count,2));
    G_weight = bsxfun(@rdivide, G_count, sum(G_count,3));

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

    for p=1:N
        for q=p+1:N
            inc = squeeze(G_weight(p,q,:))* Tp_weight(p,:);
            numer_Phi = numer_Phi + F(p,q).*inc;
            denom_Phi = denom_Phi + inc;
            inc = squeeze(G_weight(q,p,:))* Tp_weight(q,:);
            numer_Phi = numer_Phi + F(q,p).* inc;
            denom_Phi = denom_Phi + inc;
        end
    end
    Phi = numer_Phi./denom_Phi;

    params_new.Theta =Theta;
    params_new.Phi = Phi;
    params_new.Beta = Beta;
    params_new.B = B;


end