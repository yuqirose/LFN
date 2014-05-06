% The starting seed number
SEED = 1;
OUTPUT = 1;
N = int2str(252);
T = 50; 

BETA  = 1;
ALPHA = 1;
BURNIN   = 100;

%%
tic
%WS: word ID
%DS: doc ID

times = ['2010-09--2010-12';'2011-01--2011-04';'2011-05--2011-08';'2011-09--2011-12';'2012-01--2012-04'...
     ;'2012-05--2012-08';'2012-09--2012-12';'2013-01--2013-04'];
pers = [];
for t = 1: 4  
    Train_time = times(t,:);
    Test_time =  times(t,:);
    disp(Train_time);
    tic
    datafname = strcat('N=',N,'-',Train_time,'.mat');
    load(datafname); 
    [ WP, AT , Z , X ] = GibbsSamplerAT( WS , DS , AD , T , BURNIN , ALPHA , BETA , SEED , OUTPUT );
    toc
    
    [per]  =  eval_perplexity( N, Train_time, Test_time, T,  WP,Z);
    disp(per);
    pers = [pers, per];
%     save(strcat('./Result/',Train_time, '_N=',N,'_ATM.mat'),'WP','AT','Z','X','perplexity');
end
save(strcat('./Result/perplexityTrain_ATM','_N=',N,'.mat'),'pers');
