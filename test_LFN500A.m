% load the training data set:

names = {'Data/N=1225/N-2010-09--2010-12.mat'};
nFile = length(names);

N = 500;
data.W = cell(1,N);
data.D = cell(N,N);
for p=1:N
    for q=1:N
        data.D{p,q} = zeros(nFile,1);
    end
    data.W{p} = [];
end

for iFile = 1:nFile
    d = load(names{iFile});
    % read the dialog records
    AD = full(d.AD);
    disp(['time slide ' num2str(iFile) ' has ' num2str(sum(AD(:))) ' dialogs']);
    for p=1:N
        for q=1:N
            data.D{p,q}(iFile) = double(AD(p,q)>0);
        end
        data.W{p} = [data.W{p} d.WS(d.DS==p)];
    end
end

Fname = 'Data/N=1225/F_N=1225.txt';
data.F = load(Fname);
data.F = data.F(1:N,1:N);

V = 0;
for p=1:N
    V = max(V, max(data.W{p}));
end
hyper.K = 5;
hyper.N = N;
hyper.V = V;
hyper.M = nFile;

K = 5;

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

%params = load('init.mat');
%params = params.params;


[ T,G,params, LogLike_List ] = Gibbs_sampling(data, hyper, params,'fullLFNA');
save('results1.mat','T','G','params','LogLike_List');