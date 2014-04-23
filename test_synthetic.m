% test script for LFN
%% test on the synthetic dataset
clear; clc;
addpath(genpath('.'));

load('fake_data01.mat');

%%

[ T,G,params,LogLike_List ] = Gibbs_sampling(data, hyper);
save('synthetic_rslt.mat','T','G','params','LogLike_List');