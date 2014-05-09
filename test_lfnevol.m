modelnames = {
'results1.mat';
'results2.mat';
'results3.mat'};
nModel = length(modelnames);

hyper.N = 500;
hyper.K = 5;

Thetas = cell(4,1);
ThetasPrime = cell(4,1);
Betas = cell(4,1);

for modelid = 1:3
    load(modelnames{modelid});
    Thetas{modelid} = params.Theta;
    ThetasPrime{modelid} = params.Theta_prime;
    Betas{modelid} = params.Beta;
end

for modelid = 1:2
    theta1 = Thetas{modelid}';
    theta2 = Thetas{modelid+1}';
    order = perms([1 2 3 4 5]);
    M = size(order,1);
    for i=1:M
        div(i) = -sum(sum(theta1.*log(theta2(order(i,:),:)+1e-32)));
    end
    [a,b] = min(div);
    theta2 = theta2(order(b,:),:);

    t1(modelid) = mean(sum((theta1.*log(theta2+1e-32))));
    
    theta1 = ThetasPrime{modelid}';
    theta2 = ThetasPrime{modelid+1}';
    order = perms([1 2 3 4 5]);
    M = size(order,1);
    for i=1:M
        div(i) = -sum(sum(theta1.*log(theta2(order(i,:),:)+1e-32)));
    end
    [a,b] = min(div);
    theta2 = theta2(order(b,:),:);

    t2(modelid) = mean(sum((theta1.*log(theta2+1e-32))));
    
%    t2(modelid) = mean(sum((ThetasPrime{modelid}-ThetasPrime{modelid+1}).^2,2));
end
