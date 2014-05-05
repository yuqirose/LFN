function [ loglike, loglike1, loglike2, loglike3 ] = loglike_LFN( W,F,D,T,G,params,hyper )
%LOGLIKE_LFN Summary of this function goes here
%   Detailed explanation goes here
K = hyper.K; % topic/group number
N = hyper.N; % graph size
V = hyper.V; % dictionary length
M = hyper.M;
%%
Theta = params.Theta;
Theta_prime = params.Theta_prime;
Beta = params.Beta;
B = params.B;
Tau = params.Tau;

loglike = 0;
loglike1 = 0;
loglike2 = 0;
loglike3 = 0;

for p = 1:N
    % topics
    C = numel(W{p});
    for c= 1:C
        Tpc = T(c,p);%T{p}(c);
        Wpc = W{p}(c); 
        loglike = loglike + log(Theta(p,Tpc))+ log(Beta(Tpc,Wpc));
        loglike1 = loglike1 + log(Theta(p,Tpc))+ log(Beta(Tpc,Wpc));
        if(isinf(loglike1))
            %keyboard
        end
    end
    % groups
end

Tavg = cell(N,1);
for p=1:N
    Tavg{p} = histc(T(:,p),[1:K])/sum(histc(T(:,p),[1:K]));
    %Tavg{p} = histc(T{p},[1:K])/sum(histc(T{p},[1:K]));
    Tavg{p} = reshape(Tavg{p},[K,1]);
end
Gavg = cell(N,N);
for p=1:N
    for q=p+1:N
        Gavg{p,q} = histc(G{p,q},[1:K])/sum(histc(G{p,q},[1:K]));
        Gavg{p,q} = reshape(Gavg{p,q},[K,1]);
        Gavg{q,p} = histc(G{q,p},[1:K])/sum(histc(G{q,p},[1:K]));
        Gavg{q,p} = reshape(Gavg{q,p},[K,1]);
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
%             f1 = (Theta(p,:)+Theta(q,:))/2;
            f1a = Theta_prime(p,:);
            f1b = Theta_prime(q,:);
            f2a = B(Gpqm, Gqpm);
            f2b = B(Gqpm, Gpqm);

            
            loglike = loglike + log(f1a (Gpqm)+1e-32) + log(f1b(Gqpm)+1e-32);
            loglike = loglike + Dpqm*log(f2a+1e-32) + (1-Dpqm)*log(1-f2a+1e-32);
            loglike = loglike + Dqpm*log(f2b+1e-32) + (1-Dqpm)*log(1-f2b+1e-32);

            loglike2 = loglike2 + log(f1a(Gpqm)+1e-32) + log(f1b(Gqpm)+1e-32);
            loglike2 = loglike2 + Dpqm*log(f2a+1e-32) + (1-Dpqm)*log(1-f2a+1e-32);
            loglike2 = loglike2 + Dqpm*log(f2b+1e-32) + (1-Dqpm)*log(1-f2b+1e-32);
            if(isnan(loglike2))
                %keyboard
            end
        end
    end
    
    % follow , calculate average
    for q = p+1:N
        Fpq = F(p,q);
        Fqp = F(q,p);
        %Gavg = histc(G{p,q},[1:K])/sum(histc(G{p,q},[1:K]));
        f3 = 1/(1+exp(-Tau(1)*Tavg{p}'*Tavg{q} - Tau(2)*Gavg{p,q}'*Gavg{q,p}));% TBD:transpose
   
        loglike = loglike + 2*Fpq*log(f3+1e-32) + 2*(1-Fpq)*log(1-f3+1e-32);
        
        loglike3 = loglike3 + 2*Fpq*log(f3+1e-32) + 2*(1-Fpq)*log(1-f3+1e-32);
    end
%     if(isinf(loglike))
%         keyboard
%     end
end
