%% 4Opt4AttObj Analysis Script
% jacob.g.elsey@gmail.com

clc
clear 
taskType = '4Opt4Att';
% Redo analyses?
redoFlag = 0;

% Which machine are you using?
computerName =  getComputerName;
if strcmp(computerName,'stephens-macbook-pro.local');
    paths.basePath = '/Users/stephenchu/StuphornMATLAB/saccadeVigor';
    paths.dataPath = fullfile(paths.basePath,'data');
    paths.statsPath = fullfile(paths.basePath,'stats');
    paths.figurePath = fullfile(paths.basePath,'figures');
    paths.resultsPath = fullfile(paths.basePath,'results');
    
elseif strcmp(computerName,'jacobs-macbook-pro.local')| strcmp(computerName,'jacobs-mbp') == 1;
    paths.basePath = '/Users/jacobelsey/OneDrive - Johns Hopkins/Matlab/saccadeVigor';
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

    %% Fixation number rasters and SDFs
    obj = obj.samplingStrategyVigor;

%     %% Saccade vigor
%     obj = obj.saccadeVigor;
%         saccadeVigorFit((curSubj-44000),:) = {obj.saccVigor.fitresult};
% %         save(fullfile(obj.paths.statsPath,'saccadeVigorFit.mat'), 'saccadeVigorFit');

end

%% Save tables for all subjects
 
% writetable(table(transitionIndexPropOverallAll, transitionIndexPropDom, transitionIndexPropNonDom), fullfile([paths.statsPath,'/fixTransitions4Opt4AttOverall.xlsx']));

%% Population analyses

% obj = obj.plotSaccadeVigorAllSubjects;
% obj = obj.strategyBehaviorAllSubjects;




