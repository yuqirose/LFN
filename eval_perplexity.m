% ATM
% WP (word, topic), AT (author topic ), Z (topic id array ), X (author id
% array )
function [perplexity, theta, phi ]  =  eval_perplexity(N, Train_time, Test_time, nTopic, WP,Z)

datafname = strcat('N=',N,'-',Test_time,'.mat');
load(datafname); 
DS_test = DS;
WS_test = WS;

%training set
datafname = strcat('N=',N,'-',Train_time,'.mat');
load(datafname); 

beta = 1;
alpha =1;
nTerm = max(WS);
nDoc = max(DS);
phi = zeros(nTerm, nTopic);
theta = zeros (nDoc, nTopic);

for j = 1:nTopic
    phi(:,j) = (WP(:,j)+ beta )/ ( sum(WP(:,j)) + nTerm  * beta);
end

DocTopicCount = zeros(nDoc,nTopic);
for i = 1: length(Z)
    DocTopicCount(DS(i),Z(i)) = DocTopicCount(DS(i),Z(i)) +1;
end
for d = 1:nDoc
    theta(d,:) = ( DocTopicCount(d,:)+alpha ) / ( sum(DocTopicCount(d,:)) + nTopic*alpha);
end
reVal = 0;

for i = 1: length(DS_test)
    tmp = 0;
    for k = 1:nTopic
        tmp  =  tmp + theta(DS_test(i),k) * phi (WS_test(i), k);
    end
    reVal = reVal  - log(tmp);
end
perplexity = exp( reVal / length(DS_test));

end
    

