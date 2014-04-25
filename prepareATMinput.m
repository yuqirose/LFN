%PREPAREATMINPUT Summary of this function goes here
%   Detailed explanation goes here
clear;clc;
DS = [];
WS = [];
time = '2010-09--2010-12';

fname = strcat('W-small-N=30-',time,'.txt');
fid = fopen(fname);
p = 1;
V = 3000;

while ~feof(fid)
    line = fgets(fid); %# read line by line
    words = textscan(line,'%s');
    words = words{1};
    Wp = [];
    Dp = [];
    for i = 1: numel(words);
        pair = textscan(words{i},'%d,%d');
        wordID = pair{1};
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
    fprintf('%d\n',p);

    p = p+1;
 
end
WS = double(WS);

AD = load(strcat('D-small-N=30-',time,'.txt'));
% make sure there is at least author for each document
AD = AD+diag(ones(p-1,1));
AD = sparse(double(AD>0));
save(strcat('N=30-',time,'.mat'),'WS','DS','AD');
