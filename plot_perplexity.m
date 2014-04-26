% times = ['2010-09--2010-12';'2011-01--2011-04';'2011-05--2011-08';'2011-09--2011-12';'2012-01--2012-04'...
%      ;'2012-05--2012-08';'2012-09--2012-12';'2013-01--2013-04'];
% 
% pers = [];
% for t = 1: size(times,1)-1
%     Train_time = times(t,:);
%     fname = strcat('./Result/',Train_time, '_N=',N,'_ATM.mat');
%     load(fname);
%     pers = [pers, perplexity];
% end




fname = 'perplexity_LDA_N=106.mat';
per_LDA = load(fname);
per_LDA = per_LDA.pers;
disp(mean(per_LDA));
fname = 'perplexity_ATM_N=106.mat';
per_ATM = load(fname);
per_ATM = per_ATM.pers;
disp(mean(per_ATM));

 
 
% hold all;
% plot(pers,'r');
% plot(pers_lda,'b');
% 
% legend('ATM','LDA');