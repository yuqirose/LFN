% load the training data set:

names = {'Data/small_N=30/N=30-2010-09--2010-12.mat';
'Data/small_N=30/N=30-2011-05--2011-08.mat';
'Data/small_N=30/N=30-2012-01--2012-04.mat';
'Data/small_N=30/N=30-2012-09--2012-12.mat';};

%{
'Data/small_N=30/N=30-2011-01--2011-04.mat';
'Data/small_N=30/N=30-2011-05--2011-08.mat';
'Data/small_N=30/N=30-2011-09--2011-12.mat';
'Data/small_N=30/N=30-2012-01--2012-04.mat';
'Data/small_N=30/N=30-2012-05--2012-08.mat';
'Data/small_N=30/N=30-2012-09--2012-12.mat';
'Data/small_N=30/N=30-2013-01--2013-04.mat'}
%}
nFile = length(names);

N = 30;
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

Fname = 'Data/small_N=30/F_N=30.txt';
data.F = load(Fname);
keyboard
    
V = 0;
for p=1:N
    V = max(V, max(data.W{p}));
end
keyboard
hyper.K = 5;
hyper.N = N;
hyper.V = V;
hyper.M = nFile;


[ T,G,paramsL, LogLike_List ] = Gibbs_sampling(data, hyper);
