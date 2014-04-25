clear; clc;
time = '2010-09--2010-12';

datafname = strcat('N=30-',time,'.mat');
load(datafname); % V = 1000;


% The starting seed number
SEED = 1;
OUTPUT = 1;
T = 5; 

BETA  = 0.01;
ALPHA = 50/T;
BURNIN   = 100;

%%
tic
%WS: word ID
%DS: doc ID
[ WP, AT , Z , X ] = GibbsSamplerAT( WS , DS , AD , T , BURNIN , ALPHA , BETA , SEED , OUTPUT );
toc
