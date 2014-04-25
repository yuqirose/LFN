function perplexity = lfn_perplexity(params, WS_test, DS_test, hyper)

    theta = params.Theta;
    beta= params.Beta;
    beta = beta+1e-32;
    beta = bsxfun(@rdivide, beta, sum(beta,2));
    reVal = 0;

    for i = 1:length(DS_test)
        tmp = 0;
        for k = 1:hyper.K
            tmp  =  tmp + theta(DS_test(i),k) * beta(k,WS_test(i));
        end
        reVal = reVal  - log(tmp);
        if(isinf(reVal))
            keyboard
        end
    end
    perplexity = exp( reVal / length(DS_test));
end 
