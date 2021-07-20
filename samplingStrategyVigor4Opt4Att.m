function [figureHandle] = samplingStrategy4Opt4Att(allTOD,FixDuration,FixAttributeMagnitude,whichAttribute,allResp,allStimulusOnset,subjID, dataPath, figurePath, resultsPath, statsPath,...
        iSaccAmp, iPeakVel, saccAmp, saccDur, peakVel)
%% Saccade peak velocity values for each transition type

%% First steps 
% Raster plot for each trial, greyscale for sacc vigor values
% Moving mean plot for average vigor values at each saccade over the course of each trial

%% Later steps
% Sort vigor rasters by strategy type and subsequent alignments

xLimWidth = 80; % x-axis limits for all plots
sigValue = 0.05; % criterion for significant bootstrap test
sortEVFlag = 0; % Sort trials by EV
sortEVDiffFlag = 0; % Sort trials by difference in best and second-best EV available

% Possible strategies {Within-Option to each option type, Within-Attribute to each attribute type, general strategies}
stratTypes = {'wOptAmt'; 'wOptProb'; 'wOptLoss'; 'wOptDelay'; 'wAttAmt'; 'wAttProb'; 'wAttLoss'; 'wAttDelay'; 'wOpt'; 'wAtt';};
spatialStratTypes = {'V1'; 'V2'; 'V3'; 'V4'; 'H1'; 'H2'; 'H3'; 'H4';};

%% Convert cells to mat (function at bottom of script)
if iscell(allTOD)
    allTOD = convToMat(allTOD);
end

if iscell(FixDuration)
    FixDuration = convToMat(FixDuration);
end

if iscell(FixAttributeMagnitude)
    FixAttributeMagnitude = convToMat(FixAttributeMagnitude);
end

if iscell(whichAttribute)
    whichAttribute = convToMat(whichAttribute);
end

if iscell(allResp)
    allResp = convToMat(allResp);
end

if iscell(allStimulusOnset)
    allStimulusOnset = convToMat(allStimulusOnset);
end

if iscell(iSaccAmp)
    iSaccAmp = convToMat(iSaccAmp);
end

if iscell(iPeakVel)
    iPeakVel = convToMat(iPeakVel);
end

if iscell(saccAmp)
    saccAmp = convToMat(saccAmp);
end

if iscell(saccDur)
    saccDur = convToMat(saccDur);
end

if iscell(peakVel)
    peakVel = convToMat(peakVel);
end

%% Get dominated and nonDominated trials
% ResponseRemapped = createResponseLookUpTableNew(allTOD,allResp(:,2));
ResponseRemapped = allResp(:,2);

genStats.subjNum = str2num(subjID);
subjIDChar = num2str(genStats.subjNum);

% Create allTOD magnitude rankings to determine dom and nondom
allTODRanks = allTOD(:,1:16);

% Amount/prob
allTODRanks(allTODRanks == 0.1) = 5;
allTODRanks(allTODRanks == 0.3) = 4;
allTODRanks(allTODRanks == 0.5) = 3;
allTODRanks(allTODRanks == 0.7) = 2;
allTODRanks(allTODRanks == 0.9) = 1;

% Loss
lossVals = sort(unique(allTOD(:,3)));
allTODRanks(allTODRanks == lossVals(5)) = 1;
allTODRanks(allTODRanks == lossVals(4)) = 2;
allTODRanks(allTODRanks == lossVals(3)) = 3;
allTODRanks(allTODRanks == lossVals(2)) = 4;
allTODRanks(allTODRanks == lossVals(1)) = 5;

% Delay
allTODRanks(allTODRanks == 9000) = 5;
allTODRanks(allTODRanks == 7000) = 4;
allTODRanks(allTODRanks == 5000) = 3;
allTODRanks(allTODRanks == 3000) = 2;
allTODRanks(allTODRanks == 1000) = 1;

% Get NonDominated trials
nonDomTrialsIndex = find(sum(allTODRanks(:,1:4),2) == sum(allTODRanks(:,5:8),2) & ...
                             sum(allTODRanks(:,1:4),2) == sum(allTODRanks(:,9:12),2) & ... 
                             sum(allTODRanks(:,1:4),2) == sum(allTODRanks(:,13:16),2));
nonDomTrials = allTOD(nonDomTrialsIndex,:);

% Get dominated trials
domTrialsIndex = setdiff(1:size(allTODRanks,1),nonDomTrialsIndex)';
domTrials = allTOD(domTrialsIndex,:);


%% Create variables for dominated and nondominated cases

% allTOD Ranks
allTODRanksNonDom = allTODRanks(nonDomTrialsIndex,:);
allTODRanksDom = allTODRanks(domTrialsIndex,:);

% Which attribute
nonDomWhichAttribute = whichAttribute(nonDomTrialsIndex, :);
domWhichAttribute = whichAttribute(domTrialsIndex, :);

% Fixation duration
nonDomFixDuration = FixDuration(nonDomTrialsIndex, :);
domFixDuration = FixDuration(domTrialsIndex, :);

    
% All Response times
RTs = allResp(:,1) - allStimulusOnset;
nonDomRTs = RTs(nonDomTrialsIndex);
domRTs    = RTs(domTrialsIndex);

% All Number of fixations
nFixns = sum(~isnan(FixDuration),2);
nonDomNFixns = nFixns(nonDomTrialsIndex);
domNFixns    = nFixns(domTrialsIndex);

% Saccade amplitude
saccAmpNonDom = saccAmp(nonDomTrialsIndex, :);
saccAmpDom = saccAmp(domTrialsIndex, :);

% Saccade duration
saccDurNonDom = saccDur(nonDomTrialsIndex, :);
saccDurDom = saccDur(domTrialsIndex, :);

% Peak Velocity
peakVelNonDom = peakVel(nonDomTrialsIndex, :);
peakVelDom = peakVel(domTrialsIndex, :);


% %% Load choice data
% amtOptChosen    = logical(table2array(readtable(fullfile(statsPath, ['Subject', subjID,'_amtOptChosenStats.csv']))));
% probOptChosen   = logical(table2array(readtable(fullfile(statsPath, ['Subject', subjID,'_probOptChosenStats.csv']))));
% lossOptChosen   = logical(table2array(readtable(fullfile(statsPath, ['Subject', subjID,'_lossOptChosenStats.csv']))));
% delayOptChosen  = logical(table2array(readtable(fullfile(statsPath, ['Subject', subjID,'_delayOptChosenStats.csv']))));
% allChosen       = [amtOptChosen, probOptChosen, lossOptChosen, delayOptChosen];
% optTypeEVs      = table2array(readtable(fullfile(statsPath, ['Subject', subjID,'_OptTypeNonDomEVs.xlsx'])));
% optTypeEVDiffs  = table2array(readtable(fullfile(statsPath, ['Subject', subjID,'_OptTypeNonDomEVDiffs.xlsx'])));


%% Cut early trials?
cutOff = 20; %Remove first 20 trials from analysis
nonDomTrials = nonDomTrials(cutOff:end,:);
allTODRanksNonDom = allTODRanksNonDom(cutOff:end,:);
nonDomFixDuration = nonDomFixDuration(cutOff:end,:);
nonDomWhichAttribute = nonDomWhichAttribute(cutOff:end,:);
nonDomRTs = nonDomRTs(cutOff:end,:);
nonDomNFixns = nonDomNFixns(cutOff:end,:);
saccAmpNonDom = saccAmpNonDom(cutOff:end,:);
saccDurNonDom = saccDurNonDom(cutOff:end,:);
peakVelNonDom = peakVelNonDom(cutOff:end,:);

% amtOptChosen = amtOptChosen(cutOff:end,:);
% probOptChosen = probOptChosen(cutOff:end,:);
% lossOptChosen = lossOptChosen(cutOff:end,:);
% delayOptChosen = delayOptChosen(cutOff:end,:);
% allChosen = allChosen(cutOff:end,:);
% optTypeEVs = optTypeEVs(cutOff:end,:);
% optTypeEVDiffs = optTypeEVDiffs(cutOff:end,:);



%% Average values for all RTs
genStats.meanRTs = mean(RTs);

genStats.meanNonDomRTs = mean(nonDomRTs);
genStats.meanDomRTs = mean(domRTs);

%% Average and median values for all nFixations
genStats.meanNFixns = mean(nFixns);

genStats.meanNonDomNFixns = mean(nonDomNFixns);
genStats.meanDomNFixns = mean(domNFixns);

genStats.medianNFixns = median(nFixns);

genStats.medianNonDomNFixns = median(nonDomNFixns);
genStats.medianDomNFixns = median(domNFixns);

% Separate amount and probability attributes
amts = unique([allTOD(:,1);allTOD(:,5);allTOD(:,9);allTOD(:,13);]);
amts = amts(~isnan(amts));
probs = unique([allTOD(:,2);allTOD(:,6);allTOD(:,10);allTOD(:,14);]);
probs = probs(~isnan(probs));
losses = unique([allTOD(:,3);allTOD(:,7);allTOD(:,11);allTOD(:,15);]);
losses = losses(~isnan(losses));
delays = unique([allTOD(:,4);allTOD(:,8);allTOD(:,12);allTOD(:,16);]);
delays = flipud(delays(~isnan(delays)));

%% ID Spatial locations
v1 = [651,301,651,620.333333333333,651,460.666666666667,651,780];
v2 = [857.333333333333,301,857.333333333333,620.333333333333,857.333333333333,460.666666666667,857.333333333333,780];
v3 = [1063.66666666667,301,1063.66666666667,620.333333333333,1063.66666666667,460.666666666667,1063.66666666667,780];
v4 = [1270,301,1270,620.333333333333,1270,460.666666666667,1270,780];

h1 = [1270,301,857.333333333333,301,651,301,1063.66666666667,301];
h2 = [1270,460.666666666667,857.333333333333,460.666666666667,651,460.666666666667,1063.66666666667,460.666666666667];
h3 = [1270,620.333333333333,857.333333333333,620.333333333333,651,620.333333333333,1063.66666666667,620.333333333333];
h4 = [1270,780,857.333333333333,780,651,780,1063.66666666667,780];

allPos = [v1; v2; v3; v4; h1; h2; h3; h4;];

% amtsNonDom = unique([nonDomTrials(:,1);nonDomTrials(:,5);nonDomTrials(:,9);nonDomTrials(:,13);]);
% amtsNonDom = amtsNonDom(~isnan(amtsNonDom));
% probsNonDom = unique([nonDomTrials(:,2);nonDomTrials(:,6);nonDomTrials(:,10);nonDomTrials(:,14);]);
% probsNonDom = probsNonDom(~isnan(probsNonDom));
% lossesNonDom = unique([nonDomTrials(:,3);nonDomTrials(:,7);nonDomTrials(:,11);nonDomTrials(:,15);]);
% lossesNonDom = lossesNonDom(~isnan(lossesNonDom));
% delaysNonDom = unique([nonDomTrials(:,4);nonDomTrials(:,8);nonDomTrials(:,12);nonDomTrials(:,16);]);
% delaysNonDom = delaysNonDom(~isnan(delaysNonDom));
% 
% amtsDom = unique([domTrials(:,1);domTrials(:,5);domTrials(:,9);domTrials(:,13);]);
% amtsDom = amtsDom(~isnan(amtsDom));
% probsDom = unique([domTrials(:,2);domTrials(:,6);domTrials(:,10);domTrials(:,14);]);
% probsDom = probsDom(~isnan(probsDom));
% lossesDom = unique([domTrials(:,3);domTrials(:,7);domTrials(:,11);domTrials(:,15);]);
% lossesDom = lossesDom(~isnan(lossesDom));
% delaysDom = unique([domTrials(:,4);domTrials(:,8);domTrials(:,12);domTrials(:,16);]);
% delaysDom = delaysDom(~isnan(delaysDom));

% get fix durations for A1, A2, A3, A4 
fdsA1 =getFixDurs(allTOD,whichAttribute,FixDuration, 1,amts,1);
fdsA2 =getFixDurs(allTOD,whichAttribute,FixDuration, 5,amts,5);
fdsA3 =getFixDurs(allTOD,whichAttribute,FixDuration, 9,amts,9);
fdsA4 =getFixDurs(allTOD,whichAttribute,FixDuration, 13,amts,13);

fdsByAmt = [fdsA1,fdsA2,fdsA3,fdsA4];

% get fix durations for P1, P2, P3, P4
fdsP1 =getFixDurs(allTOD,whichAttribute,FixDuration, 2,probs,2);
fdsP2 =getFixDurs(allTOD,whichAttribute,FixDuration, 6,probs,6);
fdsP3 =getFixDurs(allTOD,whichAttribute,FixDuration, 10,probs,10);
fdsP4 =getFixDurs(allTOD,whichAttribute,FixDuration, 14,probs,14);

fdsByProb = [fdsP1,fdsP2,fdsP3,fdsP4];

% get fix durations for L1, L2, L3, L4
fdsL1 =getFixDurs(allTOD,whichAttribute,FixDuration, 3,losses,3);
fdsL2 =getFixDurs(allTOD,whichAttribute,FixDuration, 7,losses,7);
fdsL3 =getFixDurs(allTOD,whichAttribute,FixDuration, 11,losses,11);
fdsL4 =getFixDurs(allTOD,whichAttribute,FixDuration, 15,losses,15);

fdsByLoss = [fdsL1,fdsL2,fdsL3,fdsL4];

% get fix durations for D1, D2, D3, D4
fdsD1 =getFixDurs(allTOD,whichAttribute,FixDuration, 4, delays,4);
fdsD2 =getFixDurs(allTOD,whichAttribute,FixDuration, 8, delays,8);
fdsD3 =getFixDurs(allTOD,whichAttribute,FixDuration, 12,delays,12);
fdsD4 =getFixDurs(allTOD,whichAttribute,FixDuration, 16,delays,16);

fdsByDelay = [fdsD1,fdsD2,fdsD3,fdsD4];


% Number and proportion of fixations by attribute type
genStats.nFixnsAmt = sum(sum(~isnan(fdsByAmt),2));
genStats.nFixnsProb = sum(sum(~isnan(fdsByProb),2));
genStats.nFixnsLoss = sum(sum(~isnan(fdsByLoss),2));
genStats.nFixnsDelay = sum(sum(~isnan(fdsByDelay),2));

genStats.propFixnsAmt = genStats.nFixnsAmt/sum(nFixns);
genStats.propFixnsProb = genStats.nFixnsProb/sum(nFixns);
genStats.propFixnsLoss = genStats.nFixnsLoss/sum(nFixns);
genStats.propFixnsDelay = genStats.nFixnsDelay/sum(nFixns);

optPositions = [1,2,3,4; 5,6,7,8; 9,10,11,12; 13,14,15,16;];


%% MAIN TRIAL LOOP - NONDOMINATED TRIALS
% find transitions by trial
nonDomWhichAttribute_Filtered = zeros(size(nonDomWhichAttribute,1), size(nonDomWhichAttribute,2));
for trialNum = 1:size(nonDomWhichAttribute,1);
    iTransitionIndexOpts(trialNum, 1) = 0;
    iTransitionIndexAtts(trialNum, 1) = 0;
    spatialPos(trialNum, 1) = 0;
    tempAtts = nonDomWhichAttribute(trialNum,:);
    tempAtts(tempAtts==0) = nan; 
    tempAtts = tempAtts(~isnan(tempAtts));
    if numel(tempAtts) <= 1;
        tempAtts = 0;
    end
    firstElementNonDom(trialNum, 1) = tempAtts(1,1);
    iAllTODRanksNonDom = allTODRanksNonDom(trialNum,:);
    iMinRank = min(allTODRanksNonDom(trialNum,[1:4]));
    iAmtOpt = find(iAllTODRanksNonDom(:,[1,5,9,13]) == iMinRank);
    iProbOpt = find(iAllTODRanksNonDom(:,[1,5,9,13]+1) == iMinRank);
    iLossOpt = find(iAllTODRanksNonDom(:,[1,5,9,13]+2) == iMinRank);
    iDelayOpt = find(iAllTODRanksNonDom(:,[1,5,9,13]+3) == iMinRank);
    
    %% Get the magnitudes of each attribute for each option type
    if iAmtOpt == 1
        amtOptMags(trialNum,:) = nonDomTrials(trialNum, 1:4);
    elseif iAmtOpt == 2
        amtOptMags(trialNum,:) = nonDomTrials(trialNum, 5:8);
    elseif iAmtOpt == 3
        amtOptMags(trialNum,:) = nonDomTrials(trialNum, 9:12);
    elseif iAmtOpt == 4
        amtOptMags(trialNum,:) = nonDomTrials(trialNum, 13:16);
    else
    end
    if iProbOpt == 1
        probOptMags(trialNum,:) = nonDomTrials(trialNum, 1:4);
    elseif iProbOpt == 2
        probOptMags(trialNum,:) = nonDomTrials(trialNum, 5:8);
    elseif iProbOpt == 3
        probOptMags(trialNum,:) = nonDomTrials(trialNum, 9:12);
    elseif iProbOpt == 4
        probOptMags(trialNum,:) = nonDomTrials(trialNum, 13:16);
    else
    end
    if iLossOpt == 1
        lossOptMags(trialNum,:) = nonDomTrials(trialNum, 1:4);
    elseif iLossOpt == 2
        lossOptMags(trialNum,:) = nonDomTrials(trialNum, 5:8);
    elseif iLossOpt == 3
        lossOptMags(trialNum,:) = nonDomTrials(trialNum, 9:12);
    elseif iLossOpt == 4
        lossOptMags(trialNum,:) = nonDomTrials(trialNum, 13:16);
    else
    end    
    if iDelayOpt == 1
        delayOptMags(trialNum,:) = nonDomTrials(trialNum, 1:4);
    elseif iDelayOpt == 2
        delayOptMags(trialNum,:) = nonDomTrials(trialNum, 5:8);
    elseif iDelayOpt == 3
        delayOptMags(trialNum,:) = nonDomTrials(trialNum, 9:12);
    elseif iDelayOpt == 4
        delayOptMags(trialNum,:) = nonDomTrials(trialNum, 13:16);
    else
    end
    
    %% Get vertical position IDs
    % 1, 2, 3, 4 = vertical IDs L to R
    % 5, 6, 7, 8 = horizontal IDs Top to Bottom
%     hold on
%     scatter(nonDomTrials(trialNum,18:2:24), nonDomTrials(trialNum,19:2:25))
%     scatter(nonDomTrials(trialNum,26:2:32), nonDomTrials(trialNum,27:2:33))
%     scatter(nonDomTrials(trialNum,34:2:40), nonDomTrials(trialNum,35:2:41))
%     scatter(nonDomTrials(trialNum,42:2:48), nonDomTrials(trialNum,43:2:49))
%     close all
    
    if all(nonDomTrials(trialNum,18:2:24) == nonDomTrials(trialNum,18))
        isVert(trialNum,:) = true;
        isHorz(trialNum,:) = false;
        iVertVals = nonDomTrials(trialNum,18:8:48);
        i1Opt = find(nonDomTrials(trialNum,18:8:49) == min(iVertVals));
            iVertVals(iVertVals == min(iVertVals)) = NaN;
        i2Opt = find(nonDomTrials(trialNum,18:8:49) == min(iVertVals));
            iVertVals(iVertVals == min(iVertVals)) = NaN;
        i3Opt = find(nonDomTrials(trialNum,18:8:49) == min(iVertVals));
            iVertVals(iVertVals == min(iVertVals)) = NaN;
        i4Opt = find(nonDomTrials(trialNum,18:8:49) == min(iVertVals));
        iOpts = [i1Opt, i2Opt, i3Opt, i4Opt];
    elseif all(nonDomTrials(trialNum,19:2:25) == nonDomTrials(trialNum,19))
        isVert(trialNum,:) = false;
        isHorz(trialNum,:) = true;
        iHorzVals = nonDomTrials(trialNum,19:8:49);
        i1Opt = find(nonDomTrials(trialNum,19:8:49) == min(iHorzVals));
            iHorzVals(iHorzVals == min(iHorzVals)) = NaN;
        i2Opt = find(nonDomTrials(trialNum,19:8:49) == min(iHorzVals));
            iHorzVals(iHorzVals == min(iHorzVals)) = NaN;
        i3Opt = find(nonDomTrials(trialNum,19:8:49) == min(iHorzVals));
            iHorzVals(iHorzVals == min(iHorzVals)) = NaN;
        i4Opt = find(nonDomTrials(trialNum,19:8:49) == min(iHorzVals));
        iOpts = [i1Opt, i2Opt, i3Opt, i4Opt];

    end
    
    % 0 pad for one attribute fixated trials
    if numel(tempAtts) <= 1;
        iTransitionIndex(trialNum, 1) = 0;
        iTransitionIndexOpts(trialNum, 1) = 0;
        iTransitionIndexAtts(trialNum, 1) = 0;
        spatialPos(trialNum, 1) = 0;
    else
    for iTrialAtt1 = 1:length(tempAtts)-1;
        iTrialAtt2 = iTrialAtt1+1;
        % 4 = repeatFix - iTriallAtt values will be equal
        if tempAtts(iTrialAtt1) == tempAtts(iTrialAtt2)
            iTransitionIndex(trialNum, iTrialAtt1) = 4;
            % 1 = withinOption - sum of iTrialAtts will be either 3 (1+2), 7 (3+4), or 11 (5+6)
            %              - amount after subtraction will be 1
        elseif     sum(tempAtts(iTrialAtt1) == 1:4) == 1 & sum(tempAtts(iTrialAtt2) == 1:4) == 1 ...
                || sum(tempAtts(iTrialAtt1) == 5:8) == 1 & sum(tempAtts(iTrialAtt2) == 5:8) == 1 ...
                || sum(tempAtts(iTrialAtt1) == 9:12) == 1 & sum(tempAtts(iTrialAtt2) == 9:12) == 1 ...
                || sum(tempAtts(iTrialAtt1) == 13:16) == 1 & sum(tempAtts(iTrialAtt2) == 13:16) == 1;
                iTransitionIndex(trialNum, iTrialAtt1) = 1;
                
                % Attribute ID for within-option transitions
                if iTrialAtt1 == 1
                    nonDomWhichAttribute_Filtered(trialNum,iTrialAtt1:iTrialAtt2) = tempAtts(iTrialAtt1:iTrialAtt2);
                else
                    nonDomWhichAttribute_Filtered(trialNum,iTrialAtt2) = tempAtts(iTrialAtt2);
                end
                % Get withinAlt types
                if sum(tempAtts(iTrialAtt2) == [1,2,3,4]) == 1 %if the wOpt transition is in opt 1
                    spatialPos(trialNum,iTrialAtt1) = find(iOpts == 1);
                    iTransitionIndexOptsPos(trialNum, iTrialAtt1) = find(allTODRanksNonDom(trialNum,[1:4]) == iMinRank); % Get the position of the best attribute
                    if iTransitionIndexOptsPos(trialNum, iTrialAtt1) == 1
                        iTransitionIndexOpts(trialNum, iTrialAtt1) = 1;
                    elseif iTransitionIndexOptsPos(trialNum, iTrialAtt1) == 2
                        iTransitionIndexOpts(trialNum, iTrialAtt1) = 2;
                    elseif iTransitionIndexOptsPos(trialNum, iTrialAtt1) == 3
                        iTransitionIndexOpts(trialNum, iTrialAtt1) = 3;
                    elseif iTransitionIndexOptsPos(trialNum, iTrialAtt1) == 4
                        iTransitionIndexOpts(trialNum, iTrialAtt1) = 4;
                    end
                elseif sum(tempAtts(iTrialAtt2) == [5,6,7,8]) == 1
                    spatialPos(trialNum,iTrialAtt1) = find(iOpts == 2);
                    iTransitionIndexOptsPos(trialNum, iTrialAtt1) = find(allTODRanksNonDom(trialNum,[5:8]) == iMinRank);
                    if iTransitionIndexOptsPos(trialNum, iTrialAtt1) == 1
                        iTransitionIndexOpts(trialNum, iTrialAtt1) = 1;
                    elseif iTransitionIndexOptsPos(trialNum, iTrialAtt1) == 2
                        iTransitionIndexOpts(trialNum, iTrialAtt1) = 2;
                    elseif iTransitionIndexOptsPos(trialNum, iTrialAtt1) == 3
                        iTransitionIndexOpts(trialNum, iTrialAtt1) = 3;
                    elseif iTransitionIndexOptsPos(trialNum, iTrialAtt1) == 4
                        iTransitionIndexOpts(trialNum, iTrialAtt1) = 4;
                    end
                elseif sum(tempAtts(iTrialAtt2) == [9,10,11,12]) == 1
                    spatialPos(trialNum,iTrialAtt1) = find(iOpts == 3);
                    iTransitionIndexOptsPos(trialNum, iTrialAtt1) = find(allTODRanksNonDom(trialNum,[9:12]) == iMinRank);
                    if iTransitionIndexOptsPos(trialNum, iTrialAtt1) == 1
                        iTransitionIndexOpts(trialNum, iTrialAtt1) = 1;
                    elseif iTransitionIndexOptsPos(trialNum, iTrialAtt1) == 2
                        iTransitionIndexOpts(trialNum, iTrialAtt1) = 2;
                    elseif iTransitionIndexOptsPos(trialNum, iTrialAtt1) == 3
                        iTransitionIndexOpts(trialNum, iTrialAtt1) = 3;
                    elseif iTransitionIndexOptsPos(trialNum, iTrialAtt1) == 4
                        iTransitionIndexOpts(trialNum, iTrialAtt1) = 4;
                    end
                elseif sum(tempAtts(iTrialAtt2) == [13,14,15,16]) == 1
                    spatialPos(trialNum,iTrialAtt1) = find(iOpts == 4);
                    iTransitionIndexOptsPos(trialNum, iTrialAtt1) = find(allTODRanksNonDom(trialNum,[13:16]) == iMinRank);
                    if iTransitionIndexOptsPos(trialNum, iTrialAtt1) == 1
                        iTransitionIndexOpts(trialNum, iTrialAtt1) = 1;
                    elseif iTransitionIndexOptsPos(trialNum, iTrialAtt1) == 2
                        iTransitionIndexOpts(trialNum, iTrialAtt1) = 2;
                    elseif iTransitionIndexOptsPos(trialNum, iTrialAtt1) == 3
                        iTransitionIndexOpts(trialNum, iTrialAtt1) = 3;
                    elseif iTransitionIndexOptsPos(trialNum, iTrialAtt1) == 4
                        iTransitionIndexOpts(trialNum, iTrialAtt1) = 4;
                    end
                end
            % 2 = withinAttribute
        elseif     sum(tempAtts(iTrialAtt1) == [1,5,9,13]) == 1 & sum(tempAtts(iTrialAtt2) == [1,5,9,13]) == 1 ...
                || sum(tempAtts(iTrialAtt1) == [1,5,9,13]+1) == 1 & sum(tempAtts(iTrialAtt2) == [1,5,9,13]+1) == 1 ...
                || sum(tempAtts(iTrialAtt1) == [1,5,9,13]+2) == 1 & sum(tempAtts(iTrialAtt2) == [1,5,9,13]+2) == 1 ...
                || sum(tempAtts(iTrialAtt1) == [1,5,9,13]+3) == 1 & sum(tempAtts(iTrialAtt2) == [1,5,9,13]+3) == 1;
                iTransitionIndex(trialNum, iTrialAtt1) = 2;
                
                % Attribute ID for within-attribute transitions
                if iTrialAtt1 == 1
                    nonDomWhichAttribute_Filtered(trialNum,iTrialAtt1:iTrialAtt2) = tempAtts(iTrialAtt1:iTrialAtt2);
                else
                    nonDomWhichAttribute_Filtered(trialNum,iTrialAtt2) = tempAtts(iTrialAtt2);
                end
                
                % Get wAtt types
                if sum(tempAtts(iTrialAtt2) == [1,5,9,13]) == 1
                    iTransitionIndexAtts(trialNum, iTrialAtt1) = 1;
                elseif sum(tempAtts(iTrialAtt2) == [1,5,9,13]+1) == 1
                    iTransitionIndexAtts(trialNum, iTrialAtt1) = 2;
                elseif sum(tempAtts(iTrialAtt2) == [1,5,9,13]+2) == 1
                    iTransitionIndexAtts(trialNum, iTrialAtt1) = 3;                        
                elseif sum(tempAtts(iTrialAtt2) == [1,5,9,13]+3) == 1
                    iTransitionIndexAtts(trialNum, iTrialAtt1) = 4;
                end
                
                % Get spatial position types
                if sum(tempAtts(iTrialAtt2) == [1,2,3,4]) == 1
                    spatialPos(trialNum,iTrialAtt1) = find(iOpts == 1);
                elseif sum(tempAtts(iTrialAtt2) == [5,6,7,8]) == 1
                    spatialPos(trialNum,iTrialAtt1) = find(iOpts == 2);
                elseif sum(tempAtts(iTrialAtt2) == [9,10,11,12]) == 1
                    spatialPos(trialNum,iTrialAtt1) = find(iOpts == 3);
                elseif sum(tempAtts(iTrialAtt2) == [13,14,15,16]) == 1
                    spatialPos(trialNum,iTrialAtt1) = find(iOpts == 4);
                end
                
            % 3 = diag - other transition types must be diagonal
        else iTransitionIndex(trialNum, iTrialAtt1) = 3;
%             % Get spatial position types
%             if sum(tempAtts(iTrialAtt2) == [1,2,3,4]) == 1
%                 spatialPos(trialNum,iTrialAtt1) = find(iOpts == 1);
%             elseif sum(tempAtts(iTrialAtt2) == [5,6,7,8]) == 1
%                 spatialPos(trialNum,iTrialAtt1) = find(iOpts == 2);
%             elseif sum(tempAtts(iTrialAtt2) == [9,10,11,12]) == 1
%                 spatialPos(trialNum,iTrialAtt1) = find(iOpts == 3);
%             elseif sum(tempAtts(iTrialAtt2) == [13,14,15,16]) == 1
%                 spatialPos(trialNum,iTrialAtt1) = find(iOpts == 4);
%             end
        end
    end
    end
    iTempTransitionIndex = iTransitionIndex(trialNum,:);
    if sum(iTempTransitionIndex) == 0;
        lastTransition(trialNum,:) = 0;
    else
        lastTransition(trialNum,:) = iTempTransitionIndex(find(iTempTransitionIndex,1,'last'));
    end
end

allOptTypeMags = [amtOptMags, probOptMags, lossOptMags, delayOptMags];

% Make sure opts and atts transition sizes are same length
% Rows
if size(iTransitionIndexOpts,1) > size(iTransitionIndexAtts,1)
    iTransitionIndexAtts(size(iTransitionIndexOpts,1),:) = 0;
elseif size(iTransitionIndexOpts,1) < size(iTransitionIndexAtts,1)
    iTransitionIndexOpts(size(iTransitionIndexAtts,1),:) = 0;
else
end

% Columns
if size(iTransitionIndexOpts,2) > size(iTransitionIndexAtts,2)
    iTransitionIndexAtts(:,size(iTransitionIndexOpts,2)) = 0;
elseif size(iTransitionIndexOpts,2) < size(iTransitionIndexAtts,2)
    iTransitionIndexOpts(:,size(iTransitionIndexAtts,2)) = 0;
else
end



% % Make sure opts and atts transition sizes are same length as spatialPos
% % Rows
% if size(iTransitionIndexOpts,1) > size(iTransitionIndexAtts,1)
%     iTransitionIndexAtts(size(iTransitionIndexOpts,1),:) = 0;
%     spatialPos(size(iTransitionIndexOpts,1),:) = 0;
% elseif size(iTransitionIndexOpts,1) < size(iTransitionIndexAtts,1)
%     iTransitionIndexOpts(size(iTransitionIndexAtts,1),:) = 0;
%     spatialPos(size(iTransitionIndexAtts,1),:) = 0;
% else
% end
% 
% % Columns
% if size(iTransitionIndexOpts,2) > size(iTransitionIndexAtts,2)
%     iTransitionIndexAtts(:,size(iTransitionIndexOpts,2)) = 0;
%     spatialPos(:,size(iTransitionIndexOpts,2)) = 0;
% elseif size(iTransitionIndexOpts,2) < size(iTransitionIndexAtts,2)
%     iTransitionIndexOpts(:,size(iTransitionIndexAtts,2)) = 0;
%     spatialPos(:,size(iTransitionIndexAtts,2)) = 0;
% else
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Quantify strategy
% 1:4 = Opt Type
% 5:8 = Att Type
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get strats
for yy = 1:size(iTransitionIndexAtts,1)
    for xx = 1:size(iTransitionIndexAtts,2)
        if iTransitionIndexOpts(yy,xx) > 0
            strat(yy,xx) = iTransitionIndexOpts(yy,xx);
        elseif iTransitionIndexAtts(yy,xx) > 0
            strat(yy,xx) = iTransitionIndexAtts(yy,xx) + 4;
        else
            strat(yy,xx) = 0;
        end
    end
end

%% Get general wOpt or wAtt strats
% Remove all zeros
for ii = 1:size(strat,1)
    x = strat(ii,:);
    x(x==0)=[];
    strat(ii,1:size(x,2)) = x;
    strat(ii,size(x,2)+1:size(strat,2)) = zeros;
end
% strat(:, size(strat,2)+1:size(strat,2)+301) = 0; % nan pad end columns for circshift


for ii = 1:size(spatialPos,1)
    x = spatialPos(ii,:);
    x(x==0)=[];
    spatialPos(ii,1:size(x,2)) = x;
    spatialPos(ii,size(x,2)+1:size(spatialPos,2)) = zeros;
end
% spatialPos(:, size(spatialPos,2)+1:size(spatialPos,2)+301) = 0; % nan pad end columns for circshift

for ii = 1:size(nonDomWhichAttribute_Filtered,1)
    x = nonDomWhichAttribute_Filtered(ii,2:end);
    x(x==0)=[];
    nonDomWhichAttribute_Filtered(ii,2:size(x,2)+1) = x;
    nonDomWhichAttribute_Filtered(ii,size(x,2)+2:size(nonDomWhichAttribute_Filtered,2)) = zeros;
end

% 
% % Condense wOpt and wAtt type lengths into single epochs
% % I.e. [1 1 1] = 1
% for ii = 1:size(strat,1)
%     stratTemp = [];
%     stratTemp = strat(ii,:);
%     stratTemp(diff(stratTemp)==0) = [];
%     strat(ii,1:size(stratTemp,2)) = stratTemp;
%     strat(ii,size(stratTemp,2)+1:size(strat,2)) = zeros;
% 
% end

% chosenEVs     = optTypeEVs(allChosen);
% 
% [~,p] = sort(chosenEVs,'descend');
% chosenEVIndex = 1:length(chosenEVs);
% chosenEVIndex(p) = chosenEVIndex;
% chosenEVIndex = chosenEVIndex';
%    
% chosenEVDiffs = optTypeEVDiffs;
% [~,p] = sort(chosenEVDiffs,'descend');
% chosenEVDiffIndex = 1:length(chosenEVDiffs);
% chosenEVDiffIndex(p) = chosenEVDiffIndex;
% chosenEVDiffIndex = chosenEVDiffIndex';

%% Sort trials if needed
if sortEVFlag == 1
    strat = strat(chosenEVIndex,:);
    spatialPos = spatialPos(chosenEVIndex,:);
elseif sortEVDiffFlag == 1
    strat = strat(chosenEVDiffIndex,:);
    spatialPos = spatialPos(chosenEVIndex,:);
else
end

% Logically index different strat types
amtOptStratLog =   strat == 1;
probOptStratLog =  strat == 2;
lossOptStratLog =  strat == 3;
delayOptStratLog = strat == 4;

amtAttStratLog =   strat == 5;
probAttStratLog =  strat == 6;
lossAttStratLog =  strat == 7;
delayAttStratLog = strat == 8;

% Gen wOpt or wAtt windows
wOptStratLog = strat > 0 & strat < 5;
wAttStratLog = strat > 4;

% wOptwAtt Strat only
stratWOptWAtt = nan(size(strat));
% Within-option strat = 1
stratWOptWAtt(strat > 0 & strat < 5) = 1;
% Within-attribute strat = 2
stratWOptWAtt(strat > 4) = 2;

%% Remove nans to the left of vigor values
for ii = 1:size(peakVelNonDom,1)
    x = peakVelNonDom(ii,:);
    x(isnan(x))=[];
    peakVelNonDom(ii,1:size(x,2)) = x;
    peakVelNonDom(ii,size(x,2)+1:size(peakVelNonDom,2)) = nan;
end

%% Fixation rasters
% Load variables
figureFolder = 'samplingStratRastersVigor';
taskType = '4Opt4Att';
close all

figure('Color',[1 1 1],'Position',[1650 -50 1200 900])


%% Peak Velocity heatmap: All NonDom Trials
subplot(2,1,1)
hold on

clims = [nanmean(reshape(peakVelNonDom,[],1)) - (nanstd(reshape(peakVelNonDom,[],1))*2), nanmean(reshape(peakVelNonDom,[],1)) + (nanstd(reshape(peakVelNonDom,[],1))*2)];
imagesc(peakVelNonDom,'AlphaData',~isnan(peakVelNonDom), clims)
colorbar
% vline(1.5, 'k')

ylabel('Trial #');
title({'All Trials',''}) 
% set(gca,'XTickLabel',[]);
% xlabel('Fixation # In Trial');
xlim([0.75 xLimWidth])

peakVelNonDomMeanByFix = nanmean(peakVelNonDom,1);
peakVelNonDomMoveMeanByFix = movmean(peakVelNonDomMeanByFix, 4);

%% Moving mean of peak velocity over time
subplot(2,1,2)
p = plot(peakVelNonDom);
colorbar
ylabel('Mean Peak Velocity');
% set(gca,'XTickLabel',[]);
xlabel('Saccade # In Trial');
xlim([0.4 xLimWidth])
ylim([0 500]);
box off


%% Save fig
figureHandle = gcf;
figure(figureHandle);
annotation('textbox','string',['Subject ', subjID, '  ', taskType] ...
    ,'position',[0.01,0.92,1,0.1],'fontsize',14 ...
    ,'fontweight','bold','linestyle','none')
% annotation('textbox','string',['CHOOSE DELAY+ OPTION \n Mean # Fixations = ', num2str(fixNumNonDomDelayOptChosenMean)] ...
%     ,'position',[0.1,0.01,1,0.1],'fontsize',14 ...
%     ,'fontweight','bold','linestyle','none')
orient landscape
set(gcf,'units','pixels' ...
    ,'units','inches'...
    ,'PaperType','usletter' ...
    ,'paperposition',[.25   .25   10.5  8])

    print(fullfile(figurePath, figureFolder,['Subject ', subjID, ' Peak Velocity Rasters ', taskType,'.pdf']),'-dpdf','-r300')

close all



end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [outMat]=convToMat(inCell)
nMaxFix = cell2mat(cellfun(@size,inCell,'uniformoutput',false));
nMaxFix = max(nMaxFix (:,2));
outMat =[];
for ii = 1:size(inCell,1)
    temp = inCell{ii,1};
    if ~isempty(temp) && size(temp,2)<nMaxFix
        temp(:,end+1:nMaxFix)=nan;
    end
    outMat=[outMat;temp];
end
end