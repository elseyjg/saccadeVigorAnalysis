classdef (ConstructOnLoad = false)specificMADM4Opt4AttObj
    % specificMADM4Opt4AttObj
    
    properties
        madmData
        subjID
        paths
        gStats
        gStatsNorm
        fDomChosen
        pChooseOpt
        pChooseOptTypeEVPlot
        pChooseOptGreatestAttMag
        fNumEVAnalysis
        fTransitionsLongitudinal
        fTransitionsLongitudinalByBlock
        fTransitionsLongitudinalByTrial
        nOptsFixedLongitudinal
        nAttsFixed
        pABLongitudinal
        compEffect
        analysisResults
        predModelsInd
        predModelsIndNorm
        plotPredModelsInd
        plotPredModelsIndNorm
        ptModelFit
        fNumRasters
        samplingStrat
        plotPChooseOptAll
        plotFNumAttTypePropAllSubjects
        predModelsPop
        plotTransProps
        plotIndTrialTransProps
        plotIndTrialTransProps_wAttwAltType
        plotConsecTransitionsAll
        plotNAttsFixed
        plotRTToOptTypeAllSubjects
        plotFNumPropNonDomOptAttTypeHeatmapAllSubjects
        plotNUniqueOptsSampledAll
        saccVigor
        pInfoSampled
        plotSaccVigorAllSubjects
        pupilDiam
        samplingStratVigor
    end
    
    methods
        % Constructor
        function obj = specificMADM4Opt4AttObj(genericMADM4Opt4AttObj, subjID, paths)            
            if isa(genericMADM4Opt4AttObj,'genericMADM4Opt4AttObj')
                obj.madmData = genericMADM4Opt4AttObj.derivedData;
                obj.subjID = subjID;
                obj.paths = paths;
            else
                error('input must be of class: genMADMObj')
            end
        end
        % Call all experiment specific functions
        function obj = batchAnalyses(obj)
            obj = genStats4Opt4Att(obj.madmData.trialParams.TOD ,obj.madmData.eyeMovementData.fixDuration,...
                                    obj.madmData.eyeMovementData.fixAttributeMagnitude, obj.madmData.eyeMovementData.whichAttribute,...
                                    obj.madmData.behavioralEvents.buttonPress, obj.madmData.behavioralEvents.stimulusOnset,...
                                    obj.subjID, obj.paths.dataPath, obj.paths.figurePath, obj.paths.resultsPath);
                                

        end
        
        %Experiment specific functions  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GENERAL SUMMARY STATS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function obj = genStats(obj)
    obj.gStats = genStats4Opt4Att(obj.madmData.trialParams.TOD ,obj.madmData.eyeMovementData.fixDuration,...
        obj.madmData.eyeMovementData.fixAttributeMagnitude, obj.madmData.eyeMovementData.whichAttribute,...
        obj.madmData.behavioralEvents.buttonPress, obj.madmData.behavioralEvents.stimulusOnset,...
        obj.subjID, obj.paths.dataPath, obj.paths.figurePath, obj.paths.resultsPath, obj.paths.statsPath);
end

function obj = fixNumRasters(obj)
    obj.fNumRasters = fixNumRasters4Opt4Att(obj.madmData.trialParams.TOD ,obj.madmData.eyeMovementData.fixDuration,...
        obj.madmData.eyeMovementData.fixAttributeMagnitude, obj.madmData.eyeMovementData.whichAttribute,...
        obj.madmData.behavioralEvents.buttonPress, obj.madmData.behavioralEvents.stimulusOnset,...
        obj.subjID, obj.paths.dataPath, obj.paths.figurePath, obj.paths.resultsPath, obj.paths.statsPath);
end

function obj = saccadeVigor(obj)
    obj.saccVigor = saccadeVigor4Opt4Att(obj.madmData.trialParams.TOD ,obj.madmData.eyeMovementData.fixDuration,...
        obj.madmData.eyeMovementData.fixAttributeMagnitude, obj.madmData.eyeMovementData.whichAttribute,...
        obj.madmData.behavioralEvents.buttonPress, obj.madmData.behavioralEvents.stimulusOnset,...
        obj.subjID, obj.paths.dataPath, obj.paths.figurePath, obj.paths.resultsPath, obj.paths.statsPath,...
        obj.madmData.eyeMovementData.iSaccAmp, obj.madmData.eyeMovementData.iPeakVel,...
        obj.madmData.eyeMovementData.saccAmp, obj.madmData.eyeMovementData.saccDur, obj.madmData.eyeMovementData.peakVel);
end

function obj = samplingStrategyVigor(obj)
    obj.samplingStratVigor = samplingStrategyVigor4Opt4Att(obj.madmData.trialParams.TOD ,obj.madmData.eyeMovementData.fixDuration,...
        obj.madmData.eyeMovementData.fixAttributeMagnitude, obj.madmData.eyeMovementData.whichAttribute,...
        obj.madmData.behavioralEvents.buttonPress, obj.madmData.behavioralEvents.stimulusOnset,...
        obj.subjID, obj.paths.dataPath, obj.paths.figurePath, obj.paths.resultsPath, obj.paths.statsPath,...
        obj.madmData.eyeMovementData.iSaccAmp, obj.madmData.eyeMovementData.iPeakVel,...
        obj.madmData.eyeMovementData.saccAmp, obj.madmData.eyeMovementData.saccDur, obj.madmData.eyeMovementData.peakVel);
end

function obj = plotSaccadeVigorAllSubjects(obj)
    obj.plotSaccVigorAllSubjects = plotSaccadeVigorAllSubjects4Opt4Att(obj.subjID, obj.paths.figurePath, obj.paths.statsPath);
end

function obj = strategyBehaviorAllSubjects(obj)
    obj.stratBehaviorAllSubjects = strategyBehaviorAllSubjects4Opt4Att(obj.subjID, obj.paths.figurePath, obj.paths.statsPath);
end



    end
end