function [ T, G] = Rnd_generateLFN(W,D, K, M)
N = length(W);

T = cell(1,N);
G = cell(N,N);

for p = 1:N
    C = length(W{p});
    T{p}= zeros(C,1);
    for c = 1:C
     T{p}(c) = randi([1,K],1);
    end
    for q = 1:N
        G{p,q} =  zeros(M,1);
        for m = 1:M
            G{p,q}(m) = randi([1,K],1);
        end
    end
end


