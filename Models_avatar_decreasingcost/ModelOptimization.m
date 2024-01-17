function [x, fval, xs, fvals] = ModelOptimization(trueSeqs, probeSeqs, behChoice, behCorrect, behRT, nIter, imodel, nreplay0)

% rng;

f = @(x)ModelCostFunction(x,trueSeqs, probeSeqs, behChoice, behCorrect, behRT, imodel, nreplay0);

fval = 0;
nTask = size(trueSeqs, 2);
for i = 1 : nIter
    disp(['iter...',num2str(i)])

    %initialize starting point
    x0 = zeros(1, 2);
    rng(i,'twister');
    x0(1) = rand();
    x0(2) = exp(9 * (rand() - 0.5));
    x0(3) = rand() * (nTask - 1)/nTask + 1 / nTask;
    
    %added by JJ
    x0(4) = rand() * log(1000);
    
    lb = [0 0.01 1/nTask 0];
    ub = [1 100 1 log(1000)];

    
    % edited by Yang
    options = optimoptions('fmincon','Display','off');

    [x1,fval1] = fmincon(f,x0, [], [], [], [], lb, ub, [], options);

    xs(i,:) = x1;
    fvals(i) = fval1;
    %
    if i == 1 || fval1 < fval
        fval = fval1;
        x = x1;
    end
end