% read W file and output the bag of word
clear; clc;
N = 102;
W = cell(1,N );
fid = fopen('W_N=102.txt');

tline = fgets(fid);
p = 1;
V = 1000;
while ~feof(fid)
    line = fgets(fid); %# read line by line
    words = textscan(line,'%s');
    words = words{1};
    Wp = [];
    for i = 1: numel(words);
        pair = textscan(words{i},'%d,%d');
        wordID = pair{1};
        wordcount = pair{2};
        if(wordID>V)
            continue;
        end
        for c = 1:wordcount;
            Wp = [Wp,wordID];
        end
        
    end
    W{p} = Wp;
    p = p+1;
    fprintf('%d\n',p);
end

fclose(fid);

%%
F = load('F_N=102.txt');

%%
D = load('D_N=102.txt');
D = D(1:101,1:101);
save('toy_N=102.mat','W','F','D');


%%
N = 106;
D_all = zeros(N,N);
D_all = D0+D1+D2+D3+D4+D5+D6+D7;
total = sum(sum(D_all));
    