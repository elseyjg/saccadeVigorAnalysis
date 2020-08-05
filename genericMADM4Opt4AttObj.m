classdef (ConstructOnLoad = false)genericMADM4Opt4AttObj
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess = private)
%         edf2MatObj
%         matData % not sure we need this property
        derivedData
    end
    
    methods
        function obj = genericMADM4Opt4AttObj(edf_fname, matData_fname)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            for ii = 1:length(edf_fname)
                % extract edf data
                curEdf = Edf2Mat(edf_fname{ii});
                [iSaccBegin_, iSaccEnd_, iBlinkStart, iBlinkEnd,...
                    pos_adjusted, iFixStart, iFixEnd, FixPosInPixles,...
                    iSaccAmp,iSaccDur, iPeakVel]...
                    = SaccadeDetectionEDF2mat(curEdf);
                
                % add saccade data to object
                obj.derivedData.eyeMovementData.iSaccBegin_{ii,1}     = iSaccBegin_;
                obj.derivedData.eyeMovementData.iSaccEnd_{ii,1}       = iSaccEnd_;
                obj.derivedData.eyeMovementData.iBlinkStart{ii,1}     = iBlinkStart;
                obj.derivedData.eyeMovementData.iBlinkEnd{ii,1}       = iBlinkEnd;
                obj.derivedData.eyeMovementData.pos_adjusted{ii,1}    = pos_adjusted;
                obj.derivedData.eyeMovementData.iFixStart{ii,1}       = iFixStart;
                obj.derivedData.eyeMovementData.iFixEnd{ii,1}         = iFixEnd;
                obj.derivedData.eyeMovementData.FixPosInPixles{ii,1}  = FixPosInPixles;
                obj.derivedData.eyeMovementData.iSaccAmp{ii,1}        = iSaccAmp;
                obj.derivedData.eyeMovementData.iSaccDur{ii,1}        = iSaccDur;
                obj.derivedData.eyeMovementData.iPeakVel{ii,1}        = iPeakVel;
                edf2MatObj{ii,1}      = curEdf;
            end
            disp('Getting trial parameters and events from .mat files')
            for ii = 1:length(matData_fname)
%                 obj.matData{ii,1} = load(matData_fname{ii})
                curRD = load(matData_fname{ii},'RECORD_DATA');
                curRD = curRD.RECORD_DATA;
                [dataPath,datafile,EXT] = fileparts(matData_fname{ii});
                % get trial parameters
                [TOD trialResults] = getTrialData4Opt4Att(dataPath, [datafile,EXT]);
                %% Remove duplicate rows due to calibration bug
                if size(TOD,1) ~= size(trialResults,1)
                    if size(trialResults,1) == 50
                        [~,~,jj] = unique(TOD,'rows', 'stable');
                        TOD = TOD([true;diff(jj)~=0],:);   
                    elseif size(trialResults,1) == 80
                        [~,~,jj] = unique(TOD,'rows', 'stable');
                        TOD = TOD([true;diff(jj)~=0],:);   
                    end
                end
                
                if size(TOD,1) ~= size(trialResults,1)
                    error('TOD and TrialStart are not the same size')
                end
                
                % add trial parameters to object
                obj.derivedData.trialParams.TOD{ii,1}          = TOD;
                obj.derivedData.trialParams.trialResults{ii,1} = trialResults;
                 
                % get behavioral events
                [trialStart stimulusOnset trialEnd buttonPress calibrationTrials ] ...
                    = getEdfEvents4Opt4Att(dataPath, [datafile,EXT],edf2MatObj{ii,1});
                
                % add behavioral events to object
                obj.derivedData.behavioralEvents.trialStart{ii,1}         = trialStart;
                obj.derivedData.behavioralEvents.stimulusOnset{ii,1}      = stimulusOnset;
                obj.derivedData.behavioralEvents.trialEnd{ii,1}           = trialEnd;
                obj.derivedData.behavioralEvents.buttonPress{ii,1}        = [buttonPress buttonPress(:,1)-stimulusOnset];
                obj.derivedData.behavioralEvents.calibrationTrials{ii,1}  = calibrationTrials;
                obj.derivedData.behavioralEvents.RewardFeedbackMatData{ii,1} = ...
                    cell2mat(curRD(strmatch('Feedback Screen Onset', curRD(:,3)),4));
                obj.derivedData.behavioralEvents.RewardFeedbackTime{ii,1} = ...
                    cell2mat(curRD(strmatch('Feedback Screen Onset', curRD(:,3)),2));
                obj.derivedData.behavioralEvents.SearchArrayOn{ii,1} = ...
                    cell2mat(curRD(strmatch('Display Screen Onset', curRD(:,3)),2));
                
                obj.derivedData.eyeMovementData.maskFixationStop{ii,1} = ...
                    cell2mat(curRD(strmatch('Mask Fixation Stop', curRD(:,3)),2));
                obj.derivedData.eyeMovementData.maskFixationStart{ii,1} = ...
                    cell2mat(curRD(strmatch('Mask Fixation Start', curRD(:,3)),2));
%                 
                % Get more saccade data
                [fixDuration, fixAttributeMagnitude, ...
                    whichAttribute,fixStartTime,fixEndTime,pupilSize,...
                    saccAmp, saccDur, peakVel] ...
                    = getFixationsOfInterest4Opt4Att ...
                    (TOD,trialStart,buttonPress,...
                    obj.derivedData.eyeMovementData.iFixStart{ii,1}, ...
                    obj.derivedData.eyeMovementData.iFixEnd{ii,1}, ...
                    obj.derivedData.eyeMovementData.FixPosInPixles{ii,1}, ...
                    obj.derivedData.eyeMovementData.iSaccBegin_{ii,1}, ...
                    obj.derivedData.eyeMovementData.iSaccEnd_{ii,1}, ...
                    obj.derivedData.eyeMovementData.iSaccAmp{ii,1}, ...
                    obj.derivedData.eyeMovementData.iSaccDur{ii,1}, ...
                    obj.derivedData.eyeMovementData.iPeakVel{ii,1}, ...
                    edf2MatObj{ii,1});
                
                % add saccade data to object
                obj.derivedData.eyeMovementData.fixDuration{ii,1}           = fixDuration;
                obj.derivedData.eyeMovementData.fixAttributeMagnitude{ii,1} = fixAttributeMagnitude;
                obj.derivedData.eyeMovementData.whichAttribute{ii,1}        = whichAttribute;
                obj.derivedData.eyeMovementData.fixStartTime{ii,1}          = fixStartTime;
                obj.derivedData.eyeMovementData.fixEndTime{ii,1}            = fixEndTime; % JE 11/18/19
                obj.derivedData.eyeMovementData.pupilSize{ii,1}             = pupilSize;
                obj.derivedData.eyeMovementData.saccAmp{ii,1}               = saccAmp;
                obj.derivedData.eyeMovementData.saccDur{ii,1}               = saccDur;
                obj.derivedData.eyeMovementData.peakVel{ii,1}               = peakVel;                
                
                %% Get pupil diameter data
                edfObj = edf2MatObj{ii,1};
                
                % initialize output vars
                % Pupil diameter, timestamps, Left eye PD, Right eye PD
                diameter      = struct('t_ms',[],'L',[],'R',[]);
                % Pupil diameter units
                diameterUnit  = 'pixels';
                % intervals of interest, a table with column labels
                zeroTime_ms = 0; % edfObj.Samples.time(1);   % first timestamp of recording
                
                diameter.t_ms = edfObj.Samples.time-edfObj.Samples.time(1);
                % find which eye was being used
                whichEye = edfObj.Events.Start.eye{1};
                if isequal(whichEye,'RIGHT')
                    diameter.R = edfObj.Samples.pupilSize;
                    diameter.L = [];
                else
                    diameter.L = edfObj.Samples.pupilSize;
                    diameter.R = [];
                end
                segmentStart = []; segmentEnd    = []; 
                segmentName  = {}; segmentSource = {};
                ct = 1;
                for trl = 1:size(fixStartTime,1)
                    idx    = find(~isnan(fixStartTime(trl,:)));
                    trlFix = fixStartTime(trl, idx);
                    for curFix = 1:length(idx)
                        segmentStart(ct,1)  = fixStartTime(trl,idx(curFix))-edfObj.Samples.time(1);
                        segmentEnd(ct,1)    = segmentStart(ct,1) + fixDuration(trl,idx(curFix));
                        segmentSource{ct,1} = 'EyeLink1000';
                        switch whichAttribute(trl,idx(curFix))
                            case 1 % A1
                                curAtt = ['A ', num2str(fixAttributeMagnitude(trl,idx(curFix)))];
                            case 2 % P1
                                curAtt = ['P ', num2str(fixAttributeMagnitude(trl,idx(curFix)))];
                            case 3 % L1
                                curAtt = ['L ', num2str(fixAttributeMagnitude(trl,idx(curFix)))];
                            case 4 % D1
                                curAtt = ['D ', num2str(fixAttributeMagnitude(trl,idx(curFix)))];
                            case 5 % A2
                                curAtt = ['A ', num2str(fixAttributeMagnitude(trl,idx(curFix)))];
                            case 6 % P2
                                curAtt = ['P ', num2str(fixAttributeMagnitude(trl,idx(curFix)))];
                            case 7 % L2
                                curAtt = ['L ', num2str(fixAttributeMagnitude(trl,idx(curFix)))];
                            case 8 % D2
                                curAtt = ['D ', num2str(fixAttributeMagnitude(trl,idx(curFix)))]; 
                            case 9 % A3
                                curAtt = ['A ', num2str(fixAttributeMagnitude(trl,idx(curFix)))];
                            case 10 % P3
                                curAtt = ['P ', num2str(fixAttributeMagnitude(trl,idx(curFix)))];
                            case 11 % L3
                                curAtt = ['L ', num2str(fixAttributeMagnitude(trl,idx(curFix)))];
                            case 12 % D3
                                curAtt = ['D ', num2str(fixAttributeMagnitude(trl,idx(curFix)))];
                            case 13 % A4
                                curAtt = ['A ', num2str(fixAttributeMagnitude(trl,idx(curFix)))];
                            case 14 % P4
                                curAtt = ['P ', num2str(fixAttributeMagnitude(trl,idx(curFix)))];
                            case 15 % L4
                                curAtt = ['L ', num2str(fixAttributeMagnitude(trl,idx(curFix)))];
                            case 16 % D4
                                curAtt = ['D ', num2str(fixAttributeMagnitude(trl,idx(curFix)))];                        
                        end
                        segmentName{ct,1}   = curAtt;        
                        ct = ct + 1;
                    end
                end
                colNames = {'segmentStart','segmentEnd', 'segmentName','segmentSource'};
                segmentsTable = table(segmentStart,segmentEnd,segmentName,segmentSource ...
                    ,'VariableNames', colNames);
                % save output vars to file 
                obj.derivedData.eyeMovementData.pupilData.diameter{ii,1} = diameter;
                obj.derivedData.eyeMovementData.pupilData.zeroTime_ms{ii,1} = zeroTime_ms;
                obj.derivedData.eyeMovementData.pupilData.diameterUnit{ii,1} = diameterUnit;
                obj.derivedData.eyeMovementData.pupilData.segmentsTable{ii,1} = segmentsTable;
                
%                 save(fullfile(outputDataPath,[NAME,'_PupilD.mat']) ...
%                     , 'diameter', 'zeroTime_ms', 'diameterUnit', 'segmentsTable')
            end
            obj.derivedData.trialParams.infos_MAD_EPdata = getMADInfos4Opt4Att(obj.derivedData.trialParams.TOD...
                , cell2mat(obj.derivedData.behavioralEvents.buttonPress)...
                , obj.derivedData.eyeMovementData.whichAttribute ...
                , cell2mat(obj.derivedData.behavioralEvents.RewardFeedbackMatData));
        end
        

    end
end
