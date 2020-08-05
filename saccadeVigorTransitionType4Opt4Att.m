function [figureHandle] = saccadeVigorTransitionType4Opt4Att(allTOD,FixDuration,FixAttributeMagnitude,whichAttribute,allResp,allStimulusOnset,subjID, dataPath, figurePath, resultsPath, statsPath, iSaccAmp, iPeakVel, saccAmp, saccDur, peakVel);

%% Fixation number rasters and bar sums for
% Attribute type
% Option type

%% Bars to sdfs
% Parse trials by type chosen
%   Update prop chosen for each alignment
%   Sort trials by diff EV
%   Sort trials by first 3 epochs wAtt vs wOpt (detect shift in strat)
%   CDF of prop chosen aligned to each raster

saccadeVigorStratsFolder = 'saccadeVigorStrats';
taskType = '4Opt4Att';
xLimWidth = 80;
sigValue = 0.05;

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

% % Get 2Opt and 3Opt trials
% allTrials2OptIndex = find(isnan(allTOD(:,18)));
% allTrials3OptIndex = setdiff(1:length(allTOD),allTrials2OptIndex)';
% allTrials2Opt = allTOD(allTrials2OptIndex,:);
% allTrials3Opt = allTOD(allTrials3OptIndex,:);

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

% Make sure opts and atts transition sizes are same length
% Rows
if size(iTransitionIndexOpts,1) > size(iTransitionIndexAtts,1)
    iTransitionIndexAtts(size(iTransitionIndexOpts,1),:) = 0;
    spatialPos(size(iTransitionIndexOpts,1),:) = 0;
elseif size(iTransitionIndexOpts,1) < size(iTransitionIndexAtts,1)
    iTransitionIndexOpts(size(iTransitionIndexAtts,1),:) = 0;
    spatialPos(size(iTransitionIndexAtts,1),:) = 0;
else
end

% Columns
if size(iTransitionIndexOpts,2) > size(iTransitionIndexAtts,2)
    iTransitionIndexAtts(:,size(iTransitionIndexOpts,2)) = 0;
    spatialPos(:,size(iTransitionIndexOpts,2)) = 0;
elseif size(iTransitionIndexOpts,2) < size(iTransitionIndexAtts,2)
    iTransitionIndexOpts(:,size(iTransitionIndexAtts,2)) = 0;
    spatialPos(:,size(iTransitionIndexAtts,2)) = 0;
else
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Quantify strategy
% 1:4 = Opt Type
% 5:8 = Att Type
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get strats
for yy = 1:size(iTransitionIndexAtts,1)
    for zz = 1:size(iTransitionIndexAtts,2)
        if iTransitionIndexOpts(yy,zz) > 0
            strat(yy,zz) = iTransitionIndexOpts(yy,zz);
        elseif iTransitionIndexAtts(yy,zz) > 0
            strat(yy,zz) = iTransitionIndexAtts(yy,zz) + 4;
        else
            strat(yy,zz) = 0;
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
stratWOptWAtt(strat > 0 & strat < 5) = 1;
stratWOptWAtt(strat > 4) = 2;

allStratLogs = {wOptStratLog; wAttStratLog; amtOptStratLog; amtAttStratLog; probOptStratLog; probAttStratLog; lossOptStratLog; lossAttStratLog; delayOptStratLog; delayAttStratLog;};

% Get peak vels for each strat
for i = 1:10
    iStratLog = allStratLogs{i,:};
    
    iPeakVelStrat = peakVel(iStratLog);
    iSaccAmpStrat = saccAmp(iStratLog);
    
    iPeakVelStrat(any(isnan(iPeakVelStrat),2), :) = [];
    iSaccAmpStrat(any(isnan(iSaccAmpStrat),2), :) = [];
    
    peakVelByStratTypeCell{i,:} = iPeakVelStrat;
    saccAmpByStratTypeCell{i,:} = iSaccAmpStrat;

    rowLengthsPeakVel(i,:) = size(iPeakVelStrat, 1);
    rowLengthsSaccAmp(i,:) = size(iSaccAmpStrat, 1);

    peakVelMeans(i,:) = nanmean(iPeakVelStrat);
    saccAmpMeans(i,:) = nanmean(iSaccAmpStrat);
end

peakVelByStratType = nan(max(rowLengthsPeakVel), i);
saccAmpByStratType = nan(max(rowLengthsSaccAmp), i);

for i = 1:10
    iPeakVelByStratTypeCell = peakVelByStratTypeCell{i,:};
    iSaccAmpByStratTypeCell = saccAmpByStratTypeCell{i,:};

    peakVelByStratType(1:size(iPeakVelByStratTypeCell,1),i) = iPeakVelByStratTypeCell;
    saccAmpByStratType(1:size(iSaccAmpByStratTypeCell,1),i) = iSaccAmpByStratTypeCell;
end
%% Box plot
figure('Color',[1 1 1],'Position',[1650 -50 1200 900])
subplot(1,2,1)
hold on
cb = load('colorblind_colormap.mat');
cb = cb.colorblind;

c = [cb(6,:); cb(12,:); 0.9290, 0.6940, 0.1250; 0.9290, 0.6940, 0.1250; 0, 0.4470, 0.7410;  0, 0.4470, 0.7410; 0.8500, 0.3250, 0.0980;  0.8500, 0.3250, 0.0980; 0.4660, 0.6740, 0.1880; 0.4660, 0.6740, 0.1880;];
C = [c; ones(1,3); c];  % this is the trick for coloring the boxes
boxplot(peakVelByStratType, 'colors', C, 'plotstyle', 'compact', ...
    'labels', {'wOpt', 'wAtt', 'wOptAmt', 'wAttAmt', 'wOptProb', 'wAttProb', 'wOptLoss', 'wAttLoss', 'wOptDelay', 'wAttDelay'});
hold on;

for ii = 1:10
    plot(NaN,1,'color', c(ii,:), 'LineWidth', 4);
end
set(gca,'box','off')

% ax = gca;
% labels = string(ax.YAxis.TickLabels); % extract
% labels(2:2:end) = nan; % remove every other one
% ax.YAxis.TickLabels = labels; % set

% title('Population Non-Dominated Choice Proportions');
ylabel('Peak Velocity');
xlabel('Transition Type');
% legend({'Win+ Option','Probability+ Option','Loss+ Option','Delay+ Option'},'location', 'northeast');

%% Save prob choose all fig
figureHandle = gcf;
figure(figureHandle);
annotation('textbox','string',[subjID,' Saccade Vigor By Strat: ', taskType] ...
    ,'position',[0.1,0.9,1,0.1],'fontsize',18 ...
    ,'fontweight','bold','linestyle','none')
orient landscape
set(gcf,'units','pixels' ...
    ,'units','inches'...
    ,'PaperType','usletter' ...
    ,'paperposition',[.25   .25   10.5  8])

    print(fullfile(figurePath, saccadeVigorStratsFolder,[subjID, ' Saccade Vigor By Strat ', taskType,'.pdf']),'-dpdf','-r300')

close all


%% Saccade Amplitude X Peak Velocity Plots
hold on
figure('Color',[1 1 1],'Position',[1650 -50 1200 900])

for i = 1:2
    
    [xData, yData] = prepareCurveData( saccAmpByStratType(:,i), peakVelByStratType(:,i) );

    % Set up fittype and options.
    ft = fittype( 'a*(1-exp(-x/c))', 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.StartPoint = [523 6.8];

    % Fit model to data.
    [fitresult, gof] = fit( xData, yData, ft, opts );
    xval = min(xData):0.1:max(xData); 
    ci = predint(fitresult,xval,0.95,'observation','off');


    % figure( 'Name', 'untitled fit 1' );
    h = plot( fitresult, xData, yData);
    set(h,'color',c(i,:));

    hold on        
%     plot(xval,ci,'k--')
    l = legend( h, 'peakVel vs. iSaccAmp', 'fit', 'Location', 'NorthEast', 'Interpreter', 'none' );
    set(l,'visible','off')

    % Label axes
    xlabel( 'Saccade Amplitude', 'Interpreter', 'none' );
    ylabel( 'Peak Velocity', 'Interpreter', 'none' );
    grid off
end
set(gca,'box','off')

% saccVigor.fitresult = fitresult;
% saccVigor.gof = gof;

%% Save fig
figureFolder = 'saccadeVigor';
figureHandle = gcf;
figure(figureHandle);
annotation('textbox','string',['Subject ', subjID, ' Saccade Vigor Main Sequence By Strat',taskType] ...
    ,'position',[0.1,0.92,1,0.1],'fontsize',14 ...
    ,'fontweight','bold','linestyle','none')
% annotation('textbox','string',['CHOOSE DELAY+ OPTION \n Mean # Fixations = ', num2str(fixNumNonDomDelayOptChosenMean)] ...
%     ,'position',[0.1,0.01,1,0.1],'fontsize',14 ...
%     ,'fontweight','bold','linestyle','none')
orient landscape
set(gcf,'units','pixels' ...
    ,'units','inches'...
    ,'PaperType','usletter' ...
    ,'paperposition',[.25   .25   10.5  8])

    print(fullfile(figurePath, figureFolder,['Subject ', subjID, ' Saccade Vigor By Strat ', taskType,'.pdf']),'-dpdf','-r300')

close all




end