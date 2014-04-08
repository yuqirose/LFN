function [ loglike, loglike1, loglike2, loglike3 ] = loglike_LFN( W,F,D,T,G,params,hyper )
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
loglike1 = 0;
loglike2 = 0;
loglike3 = 0;

for p = 1:N
    % topics
    C = numel(W{p});
    for c= 1:C
        Tpc = T{p}(c);
        Wpc = W{p}(c); 
        loglike = loglike + log(Theta(p,Tpc))+ log(Beta(Tpc,Wpc));
        loglike1 = loglike1 + log(Theta(p,Tpc))+ log(Beta(Tpc,Wpc));
    end
    % groups
end
Tavg = cell(N,1);
for p=1:N
    Tavg{p} = histc(T{p},[1:K])/sum(histc(T{p},[1:K]));
end
Gavg = cell(N,N);
for p=1:N
    for q=p+1:N
        Gavg{p,q} = histc(G{p,q},[1:K])/sum(histc(G{p,q},[1:K]));
        Gavg{q,p} = histc(G{q,p},[1:K])/sum(histc(G{q,p},[1:K]));
    end
end


for p=1:N
    for q = p+1:N
        for m = 1:M
            Dpqm = D{p,q}(m);
            Dqpm = D{q,p}(m);
            Gpqm = G{p,q}(m);
            Gqpm = G{q,p}(m);
        
            % part right: (G|\theta) and (D|G,B)
            f1 = (Theta(p,:)+Theta(q,:))/2;
            f2a = B(Gpqm, Gqpm);
            f2b = B(Gqpm, Gpqm);

            
            loglike = loglike + log(f1(Gpqm)+1e-32) + log(f1(Gqpm)+1e-32);
            loglike = loglike + Dpqm*log(f2a+1e-32) + (1-Dpqm)*log(1-f2a+1e-32);
            loglike = loglike + Dqpm*log(f2b+1e-32) + (1-Dqpm)*log(1-f2b+1e-32);

            loglike2 = loglike2 + log(f1(Gpqm)+1e-32) + log(f1(Gqpm)+1e-32);
            loglike2 = loglike2 + Dpqm*log(f2a+1e-32) + (1-Dpqm)*log(1-f2a+1e-32);
            loglike2 = loglike2 + Dqpm*log(f2b+1e-32) + (1-Dqpm)*log(1-f2b+1e-32);
            if(isnan(loglike2))
                keyboard
            end
        end
    end
    
    % follow , calculate average
    for q = p+1:N
        Fpq = F(p,q);
        Fqp = F(q,p);
        %Gavg = histc(G{p,q},[1:K])/sum(histc(G{p,q},[1:K]));
        f3a = 1/(1+exp(Tavg{p}'*Phi*Gavg{p,q}));
        f3b = 1/(1+exp(Tavg{q}'*Phi*Gavg{q,p}));
   
        loglike = loglike + Fpq*log(f3a+1e-32) +(1-Fpq)*log(1-f3a+1e-32);
        loglike = loglike + Fqp*log(f3b+1e-32) +(1-Fqp)*log(1-f3b+1e-32);
        
        loglike3 = loglike3 + Fpq*log(f3a+1e-32) +(1-Fpq)*log(1-f3a+1e-32);
        loglike3 = loglike3 + Fqp*log(f3b+1e-32) +(1-Fqp)*log(1-f3b+1e-32);
    end
end

