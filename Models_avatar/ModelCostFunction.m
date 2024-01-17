function [y,dm,noEnt,intreg, preRT] = ModelCostFunction(x, trueSeqs, probeSeqs, behChoice, behCorrect, behRT, imodel)

[posSelection, taskPrediction, ~, expectedEntropy, simPos, simSteps] = Model(trueSeqs, probeSeqs, x, imodel);

%compute average task error
nTrial = size(taskPrediction, 1);

taskP = taskPrediction;
taskP(behCorrect < 1) = 1 - taskP(behCorrect < 1);
% y = -2 * nTrial * log(mean(taskP)); 

nTask = size(posSelection, 2);
steps = posSelection * [1:nTask]';

if std(expectedEntropy) > 1e-6
    expectedEntropy = (expectedEntropy- mean(expectedEntropy)) / std(expectedEntropy);
    dm = [(probeSeqs(2, :)' - steps) expectedEntropy];
    noEnt = 0;
else
    dm = [(probeSeqs(2, :)' - steps)];
    noEnt = 1;
end

% added 20231001
if imodel == 8 %only simulation, no transition to memory
    dm = probeSeqs(2, :)';
    noEnt = 1;
end

pos1trial = find(probeSeqs(2,:) == 1);
% simulatedtrial = union(find(simPos ~= probeSeqs(2,:)'),find(probeSeqs(2,:) == 1)); % based on 240 trials
simulatedtrial = find(simPos ~= probeSeqs(2,:)');
memorytrial = setdiff(find(simPos == probeSeqs(2,:)'),pos1trial);

if ismember(imodel,[3 4 5])%%% %previously only 3 4. 5 was based on the all trial
    %%% separate the three trial types: position 1 will be treated as different
    [addcol1,addcol2,addcol3] = deal(zeros(240,1));
    addcol1(simulatedtrial) = 1;
    addcol2(memorytrial) = 1;
    addcol3(pos1trial) = 1;
    dm(setdiff(1:240,simulatedtrial),1) = 0;
    dm(setdiff(1:240,memorytrial),2) = 0;
%     dm = [dm addcol1 addcol2 addcol3];
    dm = [dm addcol1 addcol2]; % we have a column of one so we just need two of them
end

if ismember(imodel,[6 7])
    %%% separate the two trial types: position 1 will be treated as both replay and memory trials
    [addcol1,addcol2] = deal(zeros(240,1));
    addcol1([pos1trial,simulatedtrial']) = 1;
    addcol2([pos1trial,memorytrial']) = 1;
    dm(setdiff(1:240,simulatedtrial),1) = 0;
    dm(setdiff(1:240,memorytrial),2) = 0;
    dm = [dm addcol1 addcol2];
end


if ismember(imodel, [2 4 5 7 8]) %%%
    %%% add the interaction term of simulation-experience x simulation step
    simtrial_accum = [];
    for i = 1:size(simSteps,1)
        simtrial_accum(i,1) = length(find(simSteps(1:i)>1));
    end
    interaction = zscore(simtrial_accum) .* zscore(simSteps);
%     interaction = simtrial_accum .* simSteps;
    interaction(setdiff(1:240,simulatedtrial),1) = 0;
    intreg = interaction;
    dm = [dm,interaction];
else
    intreg = [];
end

% add the trial trend
dm = [dm [log(1:nTrial)]'];

% add the responses and the intercept column
resps = zeros(nTrial, nTask);
for i = 1 : nTrial
    if behChoice(i) > 0
        resps(i, behChoice(i)) = 1;
    end
end
dm = [dm resps ones(nTrial, 1)];

dm = dm(behCorrect > 0, :);

behRT = behRT(behCorrect > 0);
betas = pinv(dm) * behRT;
preRT = dm * betas;
peRT = behRT - preRT;
y = (sum(peRT .* peRT) / length(behRT));
