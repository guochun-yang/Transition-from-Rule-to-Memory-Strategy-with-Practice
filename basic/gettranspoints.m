function [transpoints,condtrials,subNos] = gettranspoints(datasdir,modelresultsdir,iftpfornotransition)
% !!!should be deprecated! because the current counted transpoint is based on
% the simtrials (without error) but the trials pre/post the transpoint are 
% based on the whole 24 trial for each condition (should be OK. see below
% test)
% iftpfornotransition: 1: define 0.5 but not 24.5 as transition points; 2: define 0.5 and 24.5 as transition points; 0: define 0.5 and 24.5 as nan;
% sim_dir = ['C:\Users\guoc_yang\OneDrive - University of Iowa\JiangLab\Experiments\fMRI\fMRI_program\data_all\analysis_results_allsubs'];

if ~exist("iftpfornotransition",'var')
    iftpfornotransition = 0;
end
load([modelresultsdir '\Modelresults.mat']);
load([datasdir '\datas.mat']);
nSub = size(datas,2);
for iSub_beh = 1:nSub
    subNos(iSub_beh) = datas{iSub_beh}.subNo;
    probeSeqs = datas{iSub_beh}.probeSeqs;
    simtrials = Modelresults.simtrials{iSub_beh};
    memtrials = Modelresults.memtrials{iSub_beh};
    for iseq = 1:2
        for ipos = 1:5
            condtrial = intersect(find(probeSeqs(1,:) == iseq),find(probeSeqs(2,:) == ipos));
            condtrials(iSub_beh,iseq,ipos,:) = condtrial;
            if ipos > 1
                simtrials_cond = intersect(simtrials,condtrial);
                memtrials_cond = intersect(memtrials,condtrial);
                transpoint = length(simtrials_cond) + 0.5;
%                 if transpoint > 24; transpoint = nan; end
                if iftpfornotransition == 0
                    if ismember(transpoint,[0.5,24.5]); transpoint = nan; end
                elseif iftpfornotransition == 1
                    if ismember(transpoint,24.5); transpoint = nan; end
                end
                transpoints(iSub_beh,iseq,ipos) = transpoint;
            end
        end
    end
end
if ~isfield(Modelresults,'transpoints')
    Modelresults.transpoints = transpoints;
    save([modelresultsdir '\Modelresults.mat'],'Modelresults')
end

% test if the simtrials above is calculated from all 24 trials or only
% correct trials
for iSub_beh = 1:nSub
    nsim(iSub_beh) = length(Modelresults.simtrials{iSub_beh});
    nmem(iSub_beh) = length(Modelresults.memtrials{iSub_beh});
end
nboth = nsim + nmem;
% all are 192 = 24*8 trials, excluding the pos1. So it is correct!!!
