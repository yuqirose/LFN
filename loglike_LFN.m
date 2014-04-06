function [ loglike ] = loglike_LFN( W,F,D,T,G,params,hyper )
%LOGLIKE_LFN Summary of this function goes here
%   Detailed explanation goes here
K = hyper.K; % topic/group number
N = hyper.N; % graph size
V = hyper.V; % dictionary length
M = hyper.M;
%%
Theta = params.Theta;
Phi = params.Phi;
Beta = params.Beta;
B = params.B;

loglike = 0;

for p = 1:N
    % topics
    C = numel(W{p});
    for c= 1:C
        Tpc = T{p};
        Wpc = W{p}(c); 
        loglike = loglike + log(Theta(p,Tpc))+ log(Beta(Tpc,Wpc));
    end
    % groups
    for m = 1:M
         for q = 1:N
        
            Gpqm = G{p,q}(m);
            Gqpm = G{q,p}(m);
            f1 = (Theta(p,Gpqm)+Theta(q,Gpqm))/2;
            loglike = loglike + Gpqm*log(f1) + (1-Gpqm) *log(1-f1);
            loglike = loglike + log(B(Gpqm, Gqpm));
        end
    end
    % follow , calculate average
    Tavg = mean(T{p});
    for q = 1:N
        Fpq = F(p,q);
        Gavg = mean(G{p,q});
        f2 = 1/(1+exp(Tavg'*Phi*Gavg));
        loglike = loglike + Fpq*log(f2) +(1-Fpq)*log(1-f2);
    end
end

