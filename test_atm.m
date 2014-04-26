clear; clc;


% The starting seed number
SEED = 1;
OUTPUT = 0;
N = int2str(252);
T = 50; 

BETA  = 0.01;
ALPHA = 1;
BURNIN   = 100;

%%
tic
%WS: word ID
%DS: doc ID

times = ['2010-09--2010-12';'2011-01--2011-04';'2011-05--2011-08';'2011-09--2011-12';'2012-01--2012-04'...
     ;'2012-05--2012-08';'2012-09--2012-12';'2013-01--2013-04'];
pers = [];
% for t = 1: size(times,1)-1      
%     Train_time = times(t,:);
%     Test_time =  times(t+1,:);
%     disp(Train_time);
%     [WS,DS,AD] = prepareATMinput(N , Train_time );
%     tic
%     datafname = strcat('N=',N,'-',Train_time,'.mat');
%     load(datafname); 
%     [ WP, AT , Z , X ] = GibbsSamplerAT( WS , DS , AD , T , BURNIN , ALPHA , BETA , SEED , OUTPUT );
%     toc
%     
%     [perplexity ]  =  eval_perplexity( N, Train_time, Test_time, T,  WP,Z);
%     pers = [pers, perplexity];
%     save(strcat('./Result/',Train_time, '_N=',N,'_ATM.mat'),'WP','AT','Z','X','perplexity');
% end


for t = 1:size(times,1)
    time = times(t,:);
    disp(time);
    prepareATMinput(N , time );
end