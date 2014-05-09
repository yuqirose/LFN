names = {'Data/N=1225/N-2010-09--2010-12.mat';
'Data/N=1225/N-2011-01--2011-04.mat';
'Data/N=1225/N-2011-05--2011-08.mat';
'Data/N=1225/N-2011-09--2011-12.mat';
'Data/N=1225/N-2012-05--2012-08.mat';
'Data/N=1225/N-2013-01--2013-04.mat'};
nFile = length(names);

modelnames = {'results1.mat';
'results2.mat';
'results3.mat'};
nModel = length(modelnames);


hyper.N = 30;
hyper.K = 5;

perplexity = zeros(nFile, 4);
for iFile = 1:nFile
    load(names{iFile});  %load the data
    
    WS = WS(DS<=500);
    DS = DS(DS<=500);

    
    for modelid = 1:3
        load(modelnames{modelid}); % load the model 
        perplexity(iFile, modelid) = lfn_perplexity(params, WS, DS, hyper);
        fprintf('on data %d, modelid %d, perplexity is %d\n', iFile,  modelid, floor(perplexity(iFile, modelid))); 
    end
end
