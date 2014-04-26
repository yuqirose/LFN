%PREPAREATMINPUT Summary of this function goes here
%   Detailed explanation goes here
function [WS,DS,AD] = prepareATMinput(N , time )
DS = [];
WS = [];

% times = ['2010-09--2010-12';'2011-01--2011-04';'2011-05--2011-08';'2011-09--2011-12';'2012-01--2012-04'...
%      ;'2012-05--2012-08';'2012-09--2012-12';'2013-01--2013-04'];

fname = strcat('W-small-N=',N,'-',time,'.txt');
fid = fopen(fname);
p = 1;

while ~feof(fid)
    line = fgets(fid); %# read line by line
    words = textscan(line,'%s');
    words = words{1};
    Wp = [];
    Dp = [];
    for i = 1: numel(words);
        pair = textscan(words{i},'%d,%d');
        wordID = pair{1}+1;
        wordcount = pair{2};
%         if(wordID>V)
%             continue;
%         end
        for c = 1:wordcount;
            Wp = [Wp wordID];
            Dp = [Dp p];
        end
        
    end
    WS = [WS,Wp];
    DS = [DS,Dp];
%     fprintf('%d\n',p);

    p = p+1;
 
end
WS = double(WS);

AD = load(strcat('D-small-N=',N,'-',time,'.txt'));
% make sure there is at least author for each document
AD = AD+diag(ones(p-1,1));
AD = sparse(double(AD>0));
save(strcat('Data/small_N=252/N=',N,'-',time,'.mat'),'WS','DS','AD');
end

