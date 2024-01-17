%Input:
%trueSeqs: trueSequences of tasks, 2D array of s (sequences) x t (tasks)
%probeSeqs: probes of trials, 2D array of 2 x n (trials), first row is
%sequence, second row is position
%paras: model parameters, see comments in code for how parameters are
%organized

%Output:
%posSelection: probability of choosing a postion to start simulation, used
%to estimate RT
%taskPrediction: probability of participant selecting the correct task
%expectedEntropy is the mathmatical expectation of entropy based on the
%probabilistic distribution of posSelection

function [posSelection, taskPrediction, curPs, expectedEntropy, simPos, simSteps, simSteps_ratio, curCostss, prederr,cuedEntropy,interaction, totalStepsReplayed] = Model(trueSeqs, probeSeqs, paras, imodel, nreplay0)
% usage: [posSelection, taskPrediction, curP, expectedEntropy, simPos, simSteps, simSteps_ratio, curCosts, prederr,cuedEntropy,interaction, totalStepsReplayed] = Model(trueSeqs, probeSeqs, paras, imodel, nreplay0)
% imodel: 
% 1: base model
% 2: interaction
% 3: sim_mem_pos1
% 4: sim_mem_pos1_interaction
% 5: addRevSim
% 6: sim_mem
% 7: sim_mem_interaction
% 8: notransition

nSeq = size(trueSeqs, 1);
nPos = size(trueSeqs, 2);
nTask = nPos;
nTrial = size(probeSeqs, 2);


learningRate = paras(1);
uncertaintyWeight = paras(2);
recencyEffect = paras(3);

%added by JJ
replaySpeedup = -exp(-paras(4));

pos0Bias = 0;
otherMemory = 1 / nTask;

% current probability of performing task k when starting simulation from 
% sequence i and position j.
% Starting at location 1 will have perfect knowleage, whereas other
% locations will start with random
curP = ones(nSeq, nPos, nTask) / nTask;
for i = 1 : nSeq
    curP(i, 1, :) = 0;
    curP(i, 1, trueSeqs(i, 1)) = 1;
    curP(i, nTask, :) = (1 - recencyEffect) / (nTask - 1);
    curP(i, nTask, trueSeqs(i, nTask)) = recencyEffect;

    for j = 2 : nTask - 1
        curP(i, j, :) = (1 - otherMemory) / (nTask - 1);
        curP(i, j, trueSeqs(i, j)) = otherMemory;
    end
end

posSelection = zeros(nTrial, nPos);
taskPrediction = zeros(nTrial, 1);
entropy = ones(nTrial, nTask);
expectedEntropy = zeros(nTrial, 1);

%added by JJ, should change the number of #replayed in training phase
totalStepsReplayed = nreplay0;

for i = 1 : nTrial
    %added by JJ
    replaySpeedupDiscount = exp(totalStepsReplayed * replaySpeedup);
    
    curSeq = probeSeqs(1, i);
    curPos = probeSeqs(2, i);

    %calculating cost, which is the sum of steps of simulation and
    %uncertainty of distribution
    if imodel == 5
        curCost_frd = (1:curPos) - curPos;
        clear entropy_frd entropy_bkd
        
        for j = 1 : curPos
            idx1 = curP(curSeq, j, :) > 0;
            entropy_frd(j) = sum(curP(curSeq, j, idx1) .* log(curP(curSeq, j, idx1)));
            curCost_frd(j) = curCost_frd(j) + uncertaintyWeight * entropy_frd(j);
        end
        % back-simulation
        backDifficulty = 1; %paras(4);
        curCost_bkd = curPos - (curPos:5);
        for j0 = 1: 5-curPos+1
            j = j0+curPos-1;
            idx1 = curP(curSeq, j, :) > 0;
            entropy_bkd(j0) = sum(curP(curSeq, j, idx1) .* log(curP(curSeq, j, idx1)));
            curCost_bkd(j0) = curCost_bkd(j0) / backDifficulty + uncertaintyWeight * entropy_bkd(j0);
        end
        curCost = [curCost_frd(1:end-1), curCost_bkd];
        entropy(i,:) = [entropy_frd(1:end-1), entropy_bkd];
        simCost(i) = max(curCost_frd(1),curCost_bkd(end)); % should always choose the smaller one (note curCost is negative)
        
    else
        curCost = (1:curPos) - curPos;
        
        %added by JJ
        curCost = curCost * replaySpeedupDiscount;

        simCost(i) = curPos-1;

        for j = 1 : curPos
            idx1 = curP(curSeq, j, :) > 0;
            entropy(i, j) = sum(curP(curSeq, j, idx1) .* log(curP(curSeq, j, idx1)));
            curCost(j) = curCost(j) + uncertaintyWeight * entropy(i, j);
        end
    end
    memCost(i) = -uncertaintyWeight * entropy(i, curPos); %memCost differs from expectedEntropy, since the latter can be in the middle, the memCost means memory of the cued position

    curPs{i} = squeeze(curP(curSeq,:,:));
    curCosts{i} = curCost;
    
    %added by JJ
    [~, maxIdx] = max(curCost);
    totalStepsReplayed = totalStepsReplayed + (curPos - maxIdx(1));

    curCost = curCost - max(curCost);

    %winner take all (%think about how to measure the difficulty of making decision here, and associate that to the brain activations.)
    curCost(abs(curCost) < 1e-6) = 1; %make the largest (above -max = 0) as 1
    curCost(curCost < -1e-6) = 0;
    preP = curCost / sum(curCost) * (1 - pos0Bias); %yang: preP: predicted position. e.g., [1 0] means A2 but choose to simulate from 1
    preP(1) = preP(1) + pos0Bias;

%     %use a softmax function to convert cost to decision
%     softMax = zeros(size(curCost));
%     idx1 = curCost > -20 / softMaxBeta;
%     softMax(idx1) = exp(curCost(idx1) * softMaxBeta);
%     preP = softMax / sum(softMax) * (1 - pos0Bias);
%     preP(1) = preP(1) + pos0Bias;
    
    posSelection(i, 1:length(preP)) = preP;

    expectedEntropy(i) = sum(entropy(i, 1:length(preP)) .* preP);
    % added 20231005: calculate the entropy for the cued position, so that
    % we can calculate the increasing association related brain activations
    % (note the expEntropy in memory trial would be for the cued position,
    % but the variance would be small).
    preP2 = zeros(size(preP)); preP2(end) = 1;
    cuedEntropy(i) = sum(entropy(i, 1:length(preP)) .* preP2);
    
    for j = 1 : curPos
        taskPrediction(i) = taskPrediction(i) + ...
            preP(j) * curP(curSeq, j, trueSeqs(curSeq, j));
    end
    %taskPrediction(i, :) = taskPrediction(i, :) / curPos;

    %learning process
    curP(curSeq, curPos, :) = curP(curSeq, curPos, :) * (1 - learningRate);
    curP(curSeq, curPos, trueSeqs(curSeq, curPos)) = ...
        curP(curSeq, curPos, trueSeqs(curSeq, curPos)) + learningRate;
    prederr(i) = 1 - curP(curSeq, curPos, trueSeqs(curSeq, curPos));
end
curCostss{1} = curCosts;
curCostss{2} = simCost;
curCostss{3} = memCost;

simPos = posSelection * [1:nTask]';
simSteps = (probeSeqs(2, :)' - simPos + 1);
if imodel == 8
    simSteps = probeSeqs(2, :)';
end
simSteps_ratio = (probeSeqs(2, :)' - simPos + 1)./probeSeqs(2, :)';

simtrial_accum = [];
for i = 1:size(simSteps,1)
    simtrial_accum(i,1) = length(find(simSteps(1:i)>1));
end
interaction = zscore(simtrial_accum) .* zscore(simSteps);