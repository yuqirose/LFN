%% test on the toy real dataset
clear; clc;
addpath(genpath('.'));

nTopic = 50;
N = int2str(252);

F = load(strcat('F_N=',N,'.txt'));

times = ['2010-09--2010-12';'2011-01--2011-04';'2011-05--2011-08';'2011-09--2011-12';'2012-01--2012-04'...
     ;'2012-05--2012-08';'2012-09--2012-12';'2013-01--2013-04'];
pers = [];

for t = 1: size(times,1)-1
    Train_time = times(t,:);
    Test_time =  times(t+1,:);
    disp(Train_time);


    datafname = strcat('Data/small_N=252/N=',N,'-',Test_time,'.mat');
    load(datafname); 
    DS_test = DS;
    WS_test = WS;

    %training set
    datafname = strcat('Data/small_N=252/N=',N,'-',Train_time,'.mat');
    load(datafname); 

    fprintf('Data loaded \n');
    %% take WS and DS as input, output LFN learned model
        
    
     [T,G, params, LogLike_List] = GibbsSamplerLFN( WS , DS , AD, F,nTopic );
    fprintf('Training finished \n');
     %%
     theta = params.Theta;
     beta= params.beta;
     
     reVal = 0;
     for i = 1:length(DS_test)
        tmp = 0;
        for k = 1:nTopic
            tmp  =  tmp + theta(DS_test(i),k) * beta(WS_test(i), k);
        end
        reVal = reVal  - log(tmp);
    end
    perplexity = exp( reVal / length(DS_test));
    
    pers = [pers perplexity];  
    save(strcat('./Result/',Train_time, '_N=',N,'_LFN.mat'),'T','G','params','LogLike_List','perplexity');     
    fprinft('Result saved\n');
end 

disp(pers);
