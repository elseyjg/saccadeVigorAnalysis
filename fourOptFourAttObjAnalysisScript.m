%% 4Opt4AttObj Analysis Script
% jacob.g.elsey@gmail.com

clc
clear 
taskType = '4Opt4Att';
% Redo analyses?
redoFlag = 0;

% Which machine are you using?
computerName =  getComputerName;
if strcmp(computerName,'jacobs-macbook-pro.local')| strcmp(computerName,'jacobs-mbp') == 1;
    paths.basePath = '/Users/jacobelsey/OneDrive - Johns Hopkins/Data/madm_4Opt4Att';
    paths.dataPath = fullfile(paths.basePath,'data');
    paths.statsPath = fullfile(paths.basePath,'stats');
    paths.figurePath = fullfile(paths.basePath,'figures');
    paths.resultsPath = fullfile(paths.basePath,'results');

elseif strcmp(computerName,'mbi-vs-labpc1') == 1;
    paths.basePath   = 'C:\Users\jelsey2\Documents\Data\madm_4Opt4Att';
    paths.dataPath   = fullfile(paths.basePath,'data');
    paths.statsPath = fullfile(paths.basePath,'stats');
    paths.figurePath = fullfile(paths.basePath,'figures');
    paths.resultsPath = fullfile(paths.basePath,'results');
else
    disp('No declared path for this computer.')
    return
end

for curSubj = 44001:44040;
    if curSubj == 44001; %troubleshooting file
        continue
    end
    if curSubj == 44002; %only 80 trials, short delays
        continue
    end
    if curSubj == 44003; %only 90 trials, short delays; 44004 has 120 trials, the rest have 150
        continue
    end
    if  curSubj == 44009; %incorrect number of trials
        continue
    end
    if curSubj == 44017; %only looked at prob
        continue
    end
    if curSubj == 44022; % only looked at delay (basically) after a while, just under performance threshold
        continue
    end
    if curSubj == 44027; %no EDF files for some reason
        continue
    end
    if curSubj == 44032; % didn't pass catch trials
        continue
    end
    subjID = num2str(curSubj);
    % Find files or process new files
    if  redoFlag || ~exist(fullfile(paths.resultsPath,['Subject',subjID,'_results.mat']))
        matData_fnameStruct = dir(fullfile(paths.dataPath,['madm_4Opt4Att_',subjID,'*.mat']));
        edf_fnameStruct = dir(fullfile(paths.dataPath,['madm_4Opt4Att_',subjID,'*.edf']));
        for ii = 1:length(dir(fullfile(paths.dataPath,['madm_4Opt4Att_',subjID,'*.mat'])));
            matData_fname{ii,:} = fullfile(matData_fnameStruct(ii).folder, matData_fnameStruct(ii).name);
            edf_fname{ii,:} = fullfile(edf_fnameStruct(ii).folder, edf_fnameStruct(ii).name);
        end

        % If matData_fname is still empty, go to next subject
         if isempty(matData_fnameStruct)
            continue
        end
        derivedData = genericMADM4Opt4AttObj(edf_fname, matData_fname);
        save(fullfile(paths.resultsPath,['Subject',subjID,'_results.mat']),'derivedData');
    end
    load(fullfile(paths.resultsPath,['Subject',subjID,'_results.mat']));
    
%% Get specific analyses
     obj = specificMADM4Opt4AttObj(derivedData, subjID, paths);

%       %% General stats
%     obj = obj.genStats;
%     genStatsAll((curSubj-44000),:) = obj.gStats;
%     writetable(struct2table(genStatsAll), fullfile([paths.statsPath,'/genStats4Opt4Att.xlsx']));
%     fixNumProps((curSubj-44000),:) = [obj.gStats.fixNumProps.fixNumPropNonDomOptAttTypeMeans, obj.gStats.fixNumProps.fixNumPropNonDomOptAttTypeAmtChosenMeans,...
%         obj.gStats.fixNumProps.fixNumPropNonDomOptAttTypeProbChosenMeans, obj.gStats.fixNumProps.fixNumPropNonDomOptAttTypeLossChosenMeans, obj.gStats.fixNumProps.fixNumPropNonDomOptAttTypeDelayChosenMeans];
%     writetable(table(fixNumProps), fullfile([paths.statsPath,'/fixNumProps4Opt4Att.xlsx']));
%     fixOnAttProp((curSubj-44000),:) = struct2table([obj.gStats.fixOnAttProp]);
%     writetable(fixOnAttProp,fullfile(obj.paths.statsPath,['/fixOnAttProp.xlsx']));

%     
% %     
% %     %% General stats (normalized)
% %     obj = obj.genStatsNorm;
% %     genStatsNormAll((curSubj-44000),:) = obj.gStatsNorm;
% %     writetable(struct2table(genStatsNormAll), fullfile([paths.statsPath,'/genStatsNorm4Opt4Att.xlsx']));
% %  
% %  
%     %% Probability of choosing dom opt and nonDom opts
%     obj = obj.probChooseOpt;
%     pChooseOptAll((curSubj-44000),:) = obj.pChooseOpt;
%      writetable(struct2table(pChooseOptAll), fullfile([paths.statsPath,'/pChooseOptAll4Opt4Att.xlsx']));

%     %% Plot probability of choosing each option type AFO EV
%     obj = obj.probChooseOptTypeEVPlot;
%     chooseEVOptAll((curSubj-44000),:) = obj.pChooseOptTypeEVPlot;
%     writetable(struct2table(chooseEVOptAll), fullfile([paths.statsPath,'/chooseEVOptAll.xlsx']));

%     %% Probability of choosing nonDom opts afo fixNum to greatest attmag
%     obj = obj.probChooseOptGreatestAttMag;
%     pChooseOptGreatestAttMag((curSubj-44000),:) = obj.pChooseOptGreatestAttMag;
%     writetable(table(pChooseOptGreatestAttMag), fullfile([paths.statsPath,'/pChooseOptGreatestAttMag4Opt4Att.xlsx']));

%     %% Fixation number to each option type as a function of difference in EV
%     obj = obj.fixNumEVAnalysis;
%     fNumEVAnalysis((curSubj-44000),:) = obj.fixNumEVAnalysis;
%     writetable(table(fixNumEVAnalysis), fullfile([paths.statsPath,'/fixNumEVAnalysis.xlsx']));


    
% %%     Fixation transitions longitudinal (normalized decision time) for each transition/option/attribute type across the entire session
%     obj = obj.fixTransitionsLongitudinalByBlock;
%       transitionIndexPropOverallAll((curSubj-44000),:) = obj.fTransitionsLongitudinalByBlock.iTransitionIndexProp;
%       transitionIndexPropDom((curSubj-44000),:) = obj.fTransitionsLongitudinalByBlock.iTransitionIndexPropDom;
%       transitionIndexPropNonDom((curSubj-44000),:) = obj.fTransitionsLongitudinalByBlock.iTransitionIndexPropNonDom;
%             
%       indTrialTransitionProp_wOpt((curSubj-44000),:) = obj.fTransitionsLongitudinalByBlock.indTrialTransitionProp_wOptMean;
%       indTrialTransitionProp_wAtt((curSubj-44000),:) = obj.fTransitionsLongitudinalByBlock.indTrialTransitionProp_wAttMean;
%       indTrialTransitionProp_diag((curSubj-44000),:) = obj.fTransitionsLongitudinalByBlock.indTrialTransitionProp_diagMean;
%       indTrialTransitionProp_repeat((curSubj-44000),:) = obj.fTransitionsLongitudinalByBlock.indTrialTransitionProp_repeatMean;
%       
%       transitionIndexPropWOptAmtOpt((curSubj-44000),:) =   obj.fTransitionsLongitudinalByBlock.indTrialTransitionPropNonDom_wOptAmtOptMean;
%       transitionIndexPropWOptProbOpt((curSubj-44000),:) =  obj.fTransitionsLongitudinalByBlock.indTrialTransitionPropNonDom_wOptProbOptMean;
%       transitionIndexPropWOptLossOpt((curSubj-44000),:) =  obj.fTransitionsLongitudinalByBlock.indTrialTransitionPropNonDom_wOptLossOptMean;
%       transitionIndexPropWOptDelayOpt((curSubj-44000),:) = obj.fTransitionsLongitudinalByBlock.indTrialTransitionPropNonDom_wOptDelayOptMean;
% 
%       transitionIndexPropWAttAmtAtt((curSubj-44000),:) =   obj.fTransitionsLongitudinalByBlock.indTrialTransitionPropNonDom_wAttAmtAttMean;
%       transitionIndexPropWAttProbAtt((curSubj-44000),:) =  obj.fTransitionsLongitudinalByBlock.indTrialTransitionPropNonDom_wAttProbAttMean;
%       transitionIndexPropWAttLossAtt((curSubj-44000),:) =  obj.fTransitionsLongitudinalByBlock.indTrialTransitionPropNonDom_wAttLossAttMean;
%       transitionIndexPropWAttDelayAtt((curSubj-44000),:) = obj.fTransitionsLongitudinalByBlock.indTrialTransitionPropNonDom_wAttDelayAttMean;
% 
%       fixNumTransitionPropsChosen((curSubj-44000),:) = [obj.fTransitionsLongitudinalByBlock.propFixNumOptTypeChosenOverall, obj.fTransitionsLongitudinalByBlock.propFixNumAttTypeChosenOverall, obj.fTransitionsLongitudinalByBlock.propWOptTypeTransitionsChosenOverall, obj.fTransitionsLongitudinalByBlock.propWAttTypeTransitionsChosenOverall ,obj.fTransitionsLongitudinalByBlock.propFixNumOptAttTypeChosenOverall, obj.fTransitionsLongitudinalByBlock.propWOptWAttTypeTransitionsChosenOverall ,obj.fTransitionsLongitudinalByBlock.propFixNumTransitionOptAttTypeChosenOverall];
%       writetable(table(fixNumTransitionPropsChosen), fullfile([paths.statsPath,'/fixNumTransitionPropsChosen.xlsx']));
%         
%       allMeanOptChosenEVDiffs((curSubj-44000),:) = [obj.fTransitionsLongitudinalByBlock.allMeanOptChosenDiffs];
%       writetable(table(allMeanOptChosenEVDiffs), fullfile([paths.statsPath,'/allMeanOptChosenEVDiffs.xlsx']));
%       
%       meanAllEVsPresented((curSubj-44000),:) = [obj.fTransitionsLongitudinalByBlock.meanAllEVsPresented];
%       writetable(table(meanAllEVsPresented), fullfile([paths.statsPath,'/meanAllEVsPresented.xlsx']));     
%             
%       allMeanOptChosenEVs((curSubj-44000),:) = [obj.fTransitionsLongitudinalByBlock.allMeanOptChosenEVs];
%       writetable(table(allMeanOptChosenEVs), fullfile([paths.statsPath,'/allMeanOptChosenEVs.xlsx']));
% 
%       amtOptChosenOverallTransitionProps((curSubj-44000),:) =   obj.fTransitionsLongitudinalByBlock.amtOptChosenOverallTransitionProps;
%       probOptChosenOverallTransitionProps((curSubj-44000),:) =   obj.fTransitionsLongitudinalByBlock.probOptChosenOverallTransitionProps;
%       lossOptChosenOverallTransitionProps((curSubj-44000),:) =   obj.fTransitionsLongitudinalByBlock.lossOptChosenOverallTransitionProps;
%       delayOptChosenOverallTransitionProps((curSubj-44000),:) =   obj.fTransitionsLongitudinalByBlock.delayOptChosenOverallTransitionProps;
%       writetable(table([amtOptChosenOverallTransitionProps, probOptChosenOverallTransitionProps, lossOptChosenOverallTransitionProps, delayOptChosenOverallTransitionProps]), fullfile([paths.statsPath,'/optChosenOverallTransitionProps.xlsx']));
% 
%       EV1OptChosenOverallTransitionProps((curSubj-44000),:) =   obj.fTransitionsLongitudinalByBlock.EV1OptChosenOverallTransitionProps;
%       EV2OptChosenOverallTransitionProps((curSubj-44000),:) =   obj.fTransitionsLongitudinalByBlock.EV2OptChosenOverallTransitionProps;
%       EV3OptChosenOverallTransitionProps((curSubj-44000),:) =   obj.fTransitionsLongitudinalByBlock.EV3OptChosenOverallTransitionProps;
%       EV4OptChosenOverallTransitionProps((curSubj-44000),:) =   obj.fTransitionsLongitudinalByBlock.EV4OptChosenOverallTransitionProps;
%       writetable(table([EV1OptChosenOverallTransitionProps, EV2OptChosenOverallTransitionProps, EV3OptChosenOverallTransitionProps, EV4OptChosenOverallTransitionProps]), fullfile([paths.statsPath,'/EVOptChosenOverallTransitionProps.xlsx']));

%          %% Fixation transitions longitudinal (normalized decision time) for each transition/option/attribute type across each trial
%     obj = obj.fixTransitionsLongitudinalByTrial;
%       propFixNumAttChosen((curSubj-44000),:) =               obj.fTransitionsLongitudinalByTrial.propFixNumAttChosen;
%       propWOptTypeTransitionsChosen((curSubj-44000),:) =     obj.fTransitionsLongitudinalByTrial.propWOptTypeTransitionsChosen;
%       propWAttTypeTransitionsChosen((curSubj-44000),:) =     obj.fTransitionsLongitudinalByTrial.propWAttTypeTransitionsChosen;
%       meanDecPropTimeOfChoiceFixNumAtts((curSubj-44000),:) = obj.fTransitionsLongitudinalByTrial.meanDecPropTimeOfChoiceFixNumAtts;
%       meanDecPropTimeOfChoiceOpts((curSubj-44000),:) =       obj.fTransitionsLongitudinalByTrial.meanDecPropTimeOfChoiceOpts;
%       meanDecPropTimeOfChoiceAtts((curSubj-44000),:) =       obj.fTransitionsLongitudinalByTrial.meanDecPropTimeOfChoiceAtts;
%       propWOptWAttTypeTransitionsChosen((curSubj-44000),:) = obj.fTransitionsLongitudinalByTrial.propWOptWAttTypeTransitionsChosen;
%       propFixNumAttTypeChosenOverall ((curSubj-44000),:) =          obj.fTransitionsLongitudinalByTrial.propFixNumAttTypeChosenOverall;
%       propWOptTypeTransitionsChosenOverall((curSubj-44000),:) =     obj.fTransitionsLongitudinalByTrial.propWOptTypeTransitionsChosenOverall;
%       propWAttTypeTransitionsChosenOverall((curSubj-44000),:) =     obj.fTransitionsLongitudinalByTrial.propWAttTypeTransitionsChosenOverall;
%       propWOptWAttTypeTransitionsChosenOverall((curSubj-44000),:) = obj.fTransitionsLongitudinalByTrial.propWOptWAttTypeTransitionsChosenOverall;
%       consecTransitionMeans((curSubj-44000),:) = obj.fTransitionsLongitudinalByTrial.consecAtts;
      
%     %% Number of different options fixated longitudinal (normalized decision time)
%     obj = numOptsFixedLongitudinal(obj)    
%       meanNumOptsFixedLongitudinal{(curSubj-44000),:} =   obj.nOptsFixedLongitudinal;
%       writetable(table([meanNumOptsFixedLongitudinal]), fullfile([paths.statsPath,'/numUniqueOptsFixedLongitudinal.xlsx']));

%     %% Number of different attributes fixated
%     obj = numAttsFixed(obj);
%       nAttsFixed((curSubj-44000),:) = obj.nAttsFixed;

% 
%       
%     %% Models of choice based on individual predictors
%     obj = obj.predictModelsInd;
%         predModelIndPVals((curSubj-44000),:) = obj.predModelsInd.pVals;
%         predModelIndSlopes((curSubj-44000),:) = obj.predModelsInd.slopes;
%         predModelIndRSqs((curSubj-44000),:) = obj.predModelsInd.rSqs;
%         predModelIndCoefs((curSubj-44000),:) = obj.predModelsInd.coefs;
%         predModelIndStepwisePreds((curSubj-44000),:) = obj.predModelsInd.stepwisePreds;
% 
% %         
% %     %% Models of choice based on individual predictors (normalized)
% %     obj = obj.predictModelsIndNorm;
% %         predModelIndPValsNorm((curSubj-44000),:) = obj.predModelsIndNorm.pVals;
% %         predModelIndSlopesNorm((curSubj-44000),:) = obj.predModelsIndNorm.slopes;
% %         predModelIndRSqsNorm((curSubj-44000),:) = obj.predModelsIndNorm.rSqs;
% %         predModelIndCoefsNorm((curSubj-44000),:) = obj.predModelsIndNorm.coefs;
% %         predModelIndStepwisePredsNorm((curSubj-44000),:) = obj.predModelsIndNorm.stepwisePreds;
% % 
% % %     %% Prospect theory model fits
% % %     obj = obj.prospectTheoryModelFit;
% % 
% % %     %% Cummulative Longitudinal Means of Variables of Interest
% % %     obj = obj.cumulativeProportionsLongitudinal;

%     %% Fixation number rasters and SDFs
%     obj = obj.fixNumRasters;

%     %% Fixation number rasters and SDFs
%     obj = obj.samplingStrategy;

    %% Saccade vigor
    obj = obj.saccadeVigor;
        saccadeVigorFit((curSubj-44000),:) = {obj.saccVigor.fitresult};
%         save(fullfile(obj.paths.statsPath,'saccadeVigorFit.mat'), 'saccadeVigorFit');

%     %% Pupil diameter analysis
%     obj = obj.pupilDiameter;

end

%% Save tables for all subjects
% 
% writetable(table(transitionIndexPropOverallAll, transitionIndexPropDom, transitionIndexPropNonDom), fullfile([paths.statsPath,'/fixTransitions4Opt4AttOverall.xlsx']));
% writetable(table(indTrialTransitionProp_wOpt, indTrialTransitionProp_wAtt, indTrialTransitionProp_diag, indTrialTransitionProp_repeat), fullfile([paths.statsPath,'/indTrialTransitionProps4Opt4Att.xlsx']));
% writetable(table(transitionIndexPropWOptAmtOpt, transitionIndexPropWOptProbOpt, transitionIndexPropWOptLossOpt, transitionIndexPropWOptDelayOpt), fullfile([paths.statsPath,'/indTrialTransitionProps4Opt4Att_WOptTypes.xlsx']));
% writetable(table(transitionIndexPropWAttAmtAtt, transitionIndexPropWAttProbAtt, transitionIndexPropWAttLossAtt, transitionIndexPropWAttDelayAtt), fullfile([paths.statsPath,'/indTrialTransitionProps4Opt4Att_WAttTypes.xlsx']));

% writetable(table(propFixNumAttTypeChosenOverall, propWOptTypeTransitionsChosenOverall, propWAttTypeTransitionsChosenOverall, propWOptWAttTypeTransitionsChosenOverall, propFixNumAttChosen, propWOptTypeTransitionsChosen, propWAttTypeTransitionsChosen, propWOptWAttTypeTransitionsChosen, meanDecPropTimeOfChoiceFixNumAtts, meanDecPropTimeOfChoiceOpts, meanDecPropTimeOfChoiceAtts), fullfile([paths.statsPath,'/fixTransitions4Opt4AttByTrial.xlsx']));
% writetable(struct2table(consecTransitionMeans), fullfile([paths.statsPath,'/consecTransitionMeans.xlsx']));
% 
% writetable(table(meanNumOptsFixedLongitudinal), fullfile([paths.statsPath,'/meanNumOptsFixedLongitudinal.xlsx']));
% 
% writetable(struct2table(nAttsFixed), fullfile([paths.statsPath,'/numAttsFixed',taskType,'.xlsx']));

% 
% 
% writetable(struct2table(predModelIndPVals), 'predModelIndPVals.xlsx');
% writetable(struct2table(predModelIndSlopes), 'predModelIndSlopes.xlsx');
% writetable(struct2table(predModelIndRSqs), 'predModelIndRSqs.xlsx');
% writetable(struct2table(predModelIndCoefs), 'predModelIndCoefs.xlsx');
% writetable(struct2table(predModelIndStepwisePreds), 'predModelIndStepwisePreds.xlsx');

% writetable(struct2table(predModelIndPValsNorm), 'predModelIndPValsNorm.xlsx');
% writetable(struct2table(predModelIndSlopesNorm), 'predModelIndSlopesNorm.xlsx');
% writetable(struct2table(predModelIndRSqsNorm), 'predModelIndRSqsNorm.xlsx');
% writetable(struct2table(predModelIndCoefsNorm), 'predModelIndCoefsNorm.xlsx');
% writetable(struct2table(predModelIndStepwisePredsNorm), 'predModelIndStepwisePredsNorm.xlsx');


%% Population analyses
% obj = obj.plotProbChooseOptAll;
% obj = obj.plotPredictModelsInd;
% obj = obj.plotPredictModelsIndNorm;
% obj = obj.predictModelsPop;
% obj = obj.plotTransitionProps;
% obj = obj.plotIndTrialTransitionProps;
% obj = obj.plotIndTrialTransitionProps_wAttwAltType;
% obj = obj.plotConsecutiveTransitionsAll;
% obj = obj.plotNumAttsFixed;
% obj = obj.plotRTToOptionTypeAllSubjects;
% obj = obj.plotFixNumPropNonDomOptAttTypeHeatmapAllSubjects;
% obj = obj.plotFixNumRastersAllSubjects;
% obj = obj.plotNumUniqueOptsSampledAll;
% obj = obj.plotFixNumAttTypePropAllSubjects;
% obj = obj.plotSaccadeVigorAllSubjects;
% obj = obj.strategyBehaviorAllSubjects;




