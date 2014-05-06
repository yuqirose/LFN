load develinput.mat
N = 30;
Ws = input{1}(:,1:N);
Ts = input{2}(:,1:N);
numToken = full(input{3}(:,1:N));
Ws = full(Ws(1:max(numToken),:));
Ts = full(Ts(1:max(numToken),:));
unicW = zeros(2,N);
unicT = zeros(2,N);
for i=1:N
    unicW(1,i) = min(Ws(1:numToken(i),i));
    unicT(1,i) = min(Ts(1:numToken(i),i));
    unicW(2,i) = max(Ws(1:numToken(i),i));
    unicT(2,i) = max(Ts(1:numToken(i),i));
end

Tpcount = input{4}(:,1:N);
Gpfeature = input{5}(1:N,1:N);
F = input{6}(1:N,1:N);
Tau = input{7};
Theta = input{8}(:,1:N);
Beta = input{9};

Theta = bsxfun(@rdivide, Theta+1e-32, sum(Theta+1e-32));
Beta = bsxfun(@rdivide, Beta+1e-32, sum(Beta+1e-32,2));


startT = tic;
for round=1:5
[tp, tpcount] = GibbsTopicPosterior(Ws, Ts, numToken, Tpcount, Gpfeature, F, Tau, Theta, Beta);
end
stopT = toc(startT)

Ws = [Ws Ws Ws Ws Ws];
Ts = [Ts Ts Ts Ts Ts];
numToken = [numToken numToken numToken numToken numToken];
Tpcount = [Tpcount Tpcount Tpcount Tpcount Tpcount];
Gpfeature = [Gpfeature Gpfeature Gpfeature Gpfeature Gpfeature];
Gpfeature = [Gpfeature; Gpfeature; Gpfeature; Gpfeature; Gpfeature];
F = [F F F F F];
F = [F; F; F; F; F];
Theta = [Theta Theta Theta Theta Theta];

keyboard
startT = tic;
for round=1:5
[tp, tpcount] = GibbsTopicPosterior(Ws, Ts, numToken, Tpcount, Gpfeature, F, Tau, Theta, Beta);
end
stopT = toc(startT)
N= N*5;
tp = reshape(tp,[max(numToken),N]);
for i=1:N
    tp(1:numToken(i),i) = tp(1:numToken(i),i)+1;
end
