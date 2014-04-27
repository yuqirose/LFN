prefix = 'ModelN30D1wF/N30_iter';
nIter = 25; 
for iter = 1:nIter
    name = [prefix num2str(iter) '.mat'];
    load(name); % load the model 
    Thetas{iter} = params.Theta;
    ThetaPs{iter} = params.Theta_prime;
end
for iter = 1:nIter-1
    div(iter) = mean(sum((Thetas{iter}-Thetas{iter+1}).^2,2));
    divP(iter) = mean(sum((ThetaPs{iter}-ThetaPs{iter+1}).^2,2));
end
figure, plot(sqrt(divP),'-d','LineWidth',2,'MarkerSize',5)
hold on, plot(sqrt(div),'r-s','LineWidth',2,'MarkerSize',5)
