

% The starting seed number
SEED = 1;
OUTPUT = 0;
T = 5; 

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
for t = 1: size(times,1)-1
    Train_time = times(t,:);
    Test_time =  times(t+1,:);
    disp(Train_time);
    tic
    datafname = strcat('N=30-',Train_time,'.mat');
    load(datafname); 
    [ WP,Z ] = GibbsSamplerLDA( WS , DS , T , BURNIN , ALPHA , BETA , SEED , OUTPUT );
    toc
    
    [perplexity ]  =  eval_perplexity(Train_time, Test_time, WP,DP,Z);
    pers = [pers, perplexity];
end


