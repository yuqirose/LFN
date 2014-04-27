Theta = ones(N, K)*(1/K);
Theta_prime = ones(N, K)*(1/K);
Beta = ones(K,V)*(1/V);
Tau = rand(3,1);
B = diag(ones(K,1)*0.8) +0.1*ones(K,K);

params.Theta = Theta;
params.Theta_prime = Theta_prime;
params.Beta = Beta;
params.Tau = Tau;
params.B = B;
