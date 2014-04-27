names = {'Data/small_N=30/N=30-2010-09--2010-12.mat';
'Data/small_N=30/N=30-2011-01--2011-04.mat';
'Data/small_N=30/N=30-2011-05--2011-08.mat';
'Data/small_N=30/N=30-2011-09--2011-12.mat';
'Data/small_N=30/N=30-2012-01--2012-04.mat';
'Data/small_N=30/N=30-2012-05--2012-08.mat';
'Data/small_N=30/N=30-2012-09--2012-12.mat';
'Data/small_N=30/N=30-2013-01--2013-04.mat'};
nFile = length(names);

prefix = 'N30_iter';

perplexity = zeros(nFile, 11);
for iFile = 1:nFile
    load(names{iFile});  %load the data
    
    for iter = 1:25
        name = [prefix num2str(iter) '.mat'];
        load(name); % load the model 
        perplexity(iFile, iter) = lfn_perplexity(params, WS, DS, hyper);
        fprintf('on data %d, iteration %d, perplexity is %d\n', iFile,  iter, floor(perplexity(iFile, iter))); 
    end
end
