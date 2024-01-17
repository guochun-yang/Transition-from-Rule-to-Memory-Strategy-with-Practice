function [AICs, betas, predictedRT,numFreePara,intreg,dm,regnames,LLs] = ModelComparison(x, trueSeqs, probeSeqs, behChoice, behCorrect, behRT, imodel)

[y,dm,noEnt,intreg] = ModelCostFunction(x, trueSeqs, probeSeqs, behChoice, behCorrect, behRT, imodel);
y = sqrt(y);

% minusV = [0 1 2 3 0];
% dm0 = dm(:, 3: end-minusV(imodel)); %trialwise+1+2+3+4+5(+6)
% dm0 = [[log(1:nTrial)]' resps ones(nTrial, 1)];
dm0 = dm(:,end-6:end);
% dm0 = dm0(behCorrect > 0,:);
behRT = behRT(behCorrect > 0);
betas = pinv(dm) * behRT;
betas0 = pinv(dm0) * behRT;
preRT0 = dm0 * betas0;
peRT0 = behRT - preRT0;
y0 = sqrt(sum(peRT0 .* peRT0) / length(behRT));

nTrialFit = length(behRT);

numFreePara = size(dm, 2) - size(dm0, 2);%5; %septrial*(1+sim) + septrial*(1+ent) + interaction
AICs(1) = (numFreePara + length(x)) * 2 + 2 * nTrialFit * log(y);
AICs(2) = 2 * nTrialFit * log(y0);

%log likelihood
k1 = numFreePara + length(x);
LLs(1) = -nTrialFit/2*log(2*pi)-nTrialFit/2*log(nTrialFit*y^2/(nTrialFit-k1))-(nTrialFit-k1)/2;
k2 = 0;
LLs(2) = -nTrialFit/2*log(2*pi)-nTrialFit/2*log(nTrialFit*y0^2/(nTrialFit-k2))-(nTrialFit-k2)/2;

predictedRT = behCorrect;
predictedRT(behCorrect > 0) = dm * betas;
assignin('base','betasx',betas);
if noEnt && imodel < 8
    betas = [betas(1),nan,betas(2:end)];
end

if imodel == 1 %base
    regnames = {'simCost','expEnt','trial','resp1','resp2','resp3','resp4','resp5','one'};
elseif imodel == 2 %interaction
    regnames = {'simCost','expEnt','interaction','trial','resp1','resp2','resp3','resp4','resp5','one'};
elseif imodel == 3 %rep+mem+pos1
    regnames = {'simCost','expEnt','simtrial','memtrial','trial','resp1','resp2','resp3','resp4','resp5','one'};
elseif imodel == 4 %rep+mem+pos1+interaction
    regnames = {'simCost','expEnt','simtrial','memtrial','interaction','trial','resp1','resp2','resp3','resp4','resp5','one'};
elseif imodel == 5 %addReverse %%%
%     regnames = {'simCost','expEnt','trial','resp1','resp2','resp3','resp4','resp5','one'};
    regnames = {'simCost','expEnt','simtrial','memtrial','interaction','trial','resp1','resp2','resp3','resp4','resp5','one'};
elseif imodel == 6 %rep+mem
    regnames = {'simCost','expEnt','simtrial','memtrial','trial','resp1','resp2','resp3','resp4','resp5','one'};
elseif imodel == 7 %rep+mem+interaction
    regnames = {'simCost','expEnt','simtrial','memtrial','interaction','trial','resp1','resp2','resp3','resp4','resp5','one'};
elseif imodel == 8 %no transition
    regnames = {'simCostFrom0','interaction','trial','resp1','resp2','resp3','resp4','resp5','one'};
end

if noEnt && imodel < 8
    idx_expEnt = find(ismember('expEnt',regnames));
    regnames(idx_expEnt) = [];
end