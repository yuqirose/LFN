names = {'Data/small_N=106/N=106-2010-09--2010-12.mat';
'Data/small_N=106/N=106-2011-01--2011-04.mat';
'Data/small_N=106/N=106-2011-05--2011-08.mat';
'Data/small_N=106/N=106-2011-09--2011-12.mat';
'Data/small_N=106/N=106-2012-01--2012-04.mat';
'Data/small_N=106/N=106-2012-05--2012-08.mat';
'Data/small_N=106/N=106-2012-09--2012-12.mat';
'Data/small_N=106/N=106-2013-01--2013-04.mat'};
nFile = length(names);

prefix = 'results';
hyper.N = 106;
hyper.K = 5;

perplexity = zeros(nFile, 7);
for iFile = 1:nFile
    load(names{iFile});  %load the data
    
    for iter = 1:6
        name = [prefix num2str(iter) '.mat'];
        load(name); % load the model 
        perplexity(iFile, iter) = lfn_perplexity(params, WS, DS, hyper);
        fprintf('on data %d, iteration %d, perplexity is %d\n', iFile,  iter, floor(perplexity(iFile, iter))); 
    end
end
