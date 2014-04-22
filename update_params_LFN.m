function [ params_new ] = update_params_LFN (W,F,D, params, Tp_count, Tw_count, G_count, G )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here




% update Theta
[N,K] = size(Tp_count);
M = length(G{1,1});
for p = 1:N
    Theta(p,: ) = Tp_count(p,:)/sum(Tp_count(p,:));
end

% update Phi

numer_B = zeros(K,K);
denom_B = zeros(K,K);

for p = 1:N
    for q = 1:N
        for m = 1:M
            if q~=p
            i = G{p,q}(m);
            j = G{q,p}(m);
            numer_B(i,j) = numer_B(i,j)+ 1*D{p,q}(m);
            denom_B (i,j) = denom_B(i,j) +1; 
            end
        end
    end
end

B = numer_B./ denom_B;

% update Beta
for k = 1:K
    Beta(k,:) = Tw_count (k,:) /sum(Tw_count (k,:));
end


% update Phi
numer_Phi = zeros(K,K);
denom_Phi = zeros(K,K);

for p = 1:N
    for q = 1:N
           if q~=p
            tmp1 = Tp_count(p,:)/sum(Tp_count(p,:));
            tmp2 = squeeze(G_count(p,q,:)/sum( G_count(p,q,:)));
            numer_Phi = numer_Phi + F(p,q).* tmp1*tmp2;
            denom_Phi = denom_Phi + tmp1*tmp2; 
           end
    end
end

   
Phi = numer_Phi ./ denom_Phi;


params_new.Theta =Theta;
params_new.Phi = Phi;
params_new.Beta = Beta;
params_new.B = B;

end

