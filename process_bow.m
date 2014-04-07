% read W file and output the 
clear; clc;
N = 102;
W = cell(1,N );
fid = fopen('W_N=102.txt');

tline = fgets(fid);
p = 1;
while ~feof(fid)
    line = fgets(fid); %# read line by line
    words = textscan(line,'%s');
    words = words{1};
    Wp = [];
    for i = 1: numel(words);
        pair = textscan(words{i},'%d,%d');
        wordID = pair{1};
        wordcount = pair{2};
        for c = 1:wordcount;
            Wp = [Wp,wordID];
        end
        
    end
    W{p} = Wp;
    p = p+1;
    fprintf('%d\n',p);
end

fclose(fid);
