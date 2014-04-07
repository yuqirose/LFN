function X = sampleCat(theta, m, n)
% sample multinomial (categorical) random variables of size [mxn]
    N = m*n;
    K = length(theta);
    count = mnrnd(N,theta,1); 
    % a rown vector indicating number of examples from each group
    count = cumsum([0 count]);
    X = zeros(N,1);
    for k=1:K
        X(count(k)+1:count(k+1)) = X(count(k)+1:count(k+1))+k;
    end
    X = X(randperm(length(X)));
    X = reshape(X,[m,n]);
end

