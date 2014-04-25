%% test on the toy real dataset
clear; clc;
addpath(genpath('.'));

nTopic = 5;

F = load('F_N=30.txt');

times = ['2010-09--2010-12';'2011-01--2011-04';'2011-05--2011-08';'2011-09--2011-12';'2012-01--2012-04'...
     ;'2012-05--2012-08';'2012-09--2012-12';'2013-01--2013-04'];
pers = [];

for t = 1: size(times,1)-1
    Train_time = times(t,:);
    Test_time =  times(t+1,:);
    disp(Train_time);


    datafname = strcat('N=30-',Test_time,'.mat');
    load(datafname); 
    DS_test = DS;
    WS_test = WS;

    %training set
    datafname = strcat('N=30-',Train_time,'.mat');
    load(datafname); 


    %% take WS and DS as input, output WP, Z
        

     [T,G, params, LogLike_List] = GibbsSamplerLFN( WS , DS , AD, F,nTopic );

     theta = params.Theta;
     beta= params.beta;
     for i = 1:length(DS_test)
        tmp = 0;
        for k = 1:nTopic
            tmp  =  tmp + theta(DS_test(i),k) * beta(WS_test(i), k);
        end
        reVal = reVal  - log(tmp);
    end
    perplexity = exp( reVal / length(DS_test));
    pers = [pers perplexity];
end 
