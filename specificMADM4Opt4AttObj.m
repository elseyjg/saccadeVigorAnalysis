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
% function obj = genStats4Opt4Att(obj, TOD, fixDuration,...
%                 fixAttributeMagnitude, whichAttribute,...
%                 buttonPress, stimulusOnset,...
%                 subjID, dataPath, figurePath, resultsPath)

function obj = genStats(obj)
    obj.gStats = genStats4Opt4Att(obj.madmData.trialParams.TOD ,obj.madmData.eyeMovementData.fixDuration,...
        obj.madmData.eyeMovementData.fixAttributeMagnitude, obj.madmData.eyeMovementData.whichAttribute,...
        obj.madmData.behavioralEvents.buttonPress, obj.madmData.behavioralEvents.stimulusOnset,...
        obj.subjID, obj.paths.dataPath, obj.paths.figurePath, obj.paths.resultsPath, obj.paths.statsPath);
end

function obj = genStatsNorm(obj)
    obj.gStatsNorm = genStatsNorm4Opt4Att(obj.madmData.trialParams.TOD ,obj.madmData.eyeMovementData.fixDuration,...
        obj.madmData.eyeMovementData.fixAttributeMagnitude, obj.madmData.eyeMovementData.whichAttribute,...
        obj.madmData.behavioralEvents.buttonPress, obj.madmData.behavioralEvents.stimulusOnset,...
        obj.subjID, obj.paths.dataPath, obj.paths.figurePath, obj.paths.resultsPath, obj.paths.statsPath);
end
        
function obj = freqDomChosen(obj)
    obj.fDomChosen = freqDomChosen4Opt4Att(obj.madmData.trialParams.TOD ,obj.madmData.eyeMovementData.fixDuration,...
        obj.madmData.eyeMovementData.fixAttributeMagnitude, obj.madmData.eyeMovementData.whichAttribute,...
        obj.madmData.behavioralEvents.buttonPress, obj.madmData.behavioralEvents.stimulusOnset,...
        obj.subjID, obj.paths.dataPath, obj.paths.figurePath, obj.paths.resultsPath);
end

function obj = probChooseOpt(obj)
    obj.pChooseOpt = probChooseOpt4Opt4Att(obj.madmData.trialParams.infos_MAD_EPdata, obj.madmData.trialParams.TOD ,obj.madmData.eyeMovementData.fixDuration,...
        obj.madmData.eyeMovementData.fixAttributeMagnitude, obj.madmData.eyeMovementData.whichAttribute,...
        obj.madmData.behavioralEvents.buttonPress, obj.madmData.behavioralEvents.stimulusOnset,...
        obj.subjID, obj.paths.dataPath, obj.paths.figurePath, obj.paths.resultsPath, obj.paths.statsPath);
end

function obj = probChooseOptTypeEVPlot(obj)
    obj.pChooseOptTypeEVPlot = probChooseOptTypeEVPlot4Opt4Att(obj.madmData.trialParams.infos_MAD_EPdata, obj.madmData.trialParams.TOD ,obj.madmData.eyeMovementData.fixDuration,...
        obj.madmData.eyeMovementData.fixAttributeMagnitude, obj.madmData.eyeMovementData.whichAttribute,...
        obj.madmData.behavioralEvents.buttonPress, obj.madmData.behavioralEvents.stimulusOnset,...
        obj.subjID, obj.paths.dataPath, obj.paths.figurePath, obj.paths.resultsPath, obj.paths.statsPath);
end

function obj = probChooseOptGreatestAttMag(obj)
    obj.pChooseOptGreatestAttMag = probChooseOptGreatestAttMag4Opt4Att(obj.madmData.trialParams.infos_MAD_EPdata, obj.madmData.trialParams.TOD ,obj.madmData.eyeMovementData.fixDuration,...
        obj.madmData.eyeMovementData.fixAttributeMagnitude, obj.madmData.eyeMovementData.whichAttribute,...
        obj.madmData.behavioralEvents.buttonPress, obj.madmData.behavioralEvents.stimulusOnset,...
        obj.subjID, obj.paths.dataPath, obj.paths.figurePath, obj.paths.resultsPath, obj.paths.statsPath);
end

function obj = fixNumEVAnalysis(obj)
    obj.fNumEVAnalysis = fixNumEVAnalysis4Opt4Att(obj.madmData.trialParams.infos_MAD_EPdata, obj.madmData.trialParams.TOD ,obj.madmData.eyeMovementData.fixDuration,...
        obj.madmData.eyeMovementData.fixAttributeMagnitude, obj.madmData.eyeMovementData.whichAttribute,...
        obj.madmData.behavioralEvents.buttonPress, obj.madmData.behavioralEvents.stimulusOnset,...
        obj.subjID, obj.paths.dataPath, obj.paths.figurePath, obj.paths.resultsPath, obj.paths.statsPath);
end

function obj = fixTransitionsLongitudinal(obj)
    obj.fTransitionsLongitudinal = fixTransitionsLongitudinal4Opt4Att(obj.madmData.trialParams.TOD,...
        obj.madmData.eyeMovementData.whichAttribute,obj.subjID)
end
function obj = fixTransitionsLongitudinalByBlock(obj)
    obj.fTransitionsLongitudinalByBlock = fixTransitionsLongitudinal4Opt4AttByBlock(obj.madmData.trialParams.TOD,obj.madmData.trialParams.infos_MAD_EPdata,...
        obj.madmData.eyeMovementData.whichAttribute, obj.madmData.eyeMovementData.fixDuration, obj.subjID, obj.paths.figurePath, obj.paths.statsPath)
end
function obj = fixTransitionsLongitudinalByTrial(obj)
    obj.fTransitionsLongitudinalByTrial = fixTransitionsLongitudinal4Opt4AttByTrial(obj.madmData.trialParams.TOD,...
        obj.madmData.eyeMovementData.whichAttribute,obj.subjID, obj.paths.figurePath, obj.paths.statsPath)
end
function obj = numOptsFixedLongitudinal(obj)
    obj.nOptsFixedLongitudinal = numOptsFixedLongitudinal4Opt4AttByBlock(obj.madmData.trialParams.TOD,...
        obj.madmData.eyeMovementData.whichAttribute,obj.subjID, obj.paths.dataPath, obj.paths.figurePath)
end   
function obj = numAttsFixed(obj)
    obj.nAttsFixed = numAttsFixed4Opt4Att(obj.madmData.trialParams.TOD,...
        obj.madmData.eyeMovementData.whichAttribute,obj.subjID, obj.paths.dataPath, obj.paths.figurePath, obj.paths.statsPath)
end   

function obj = preferenceABLongitudinal(obj)
    obj.pABLongitudinal = preferenceABLongitudinalBlock(obj.madmData.trialParams.infos_MAD_EPdata ,obj.subjID)
end

function obj = compromiseEffect(obj)
    obj.compEffect = compromiseEffectBlock(obj.madmData.trialParams.infos_MAD_EPdata ,obj.subjID)
end

function obj = predictModelsInd(obj)
    obj.predModelsInd = predictModelsInd4Opt4Att(obj.madmData.trialParams.TOD ,obj.madmData.eyeMovementData.fixDuration,...
        obj.madmData.eyeMovementData.fixAttributeMagnitude, obj.madmData.eyeMovementData.whichAttribute,...
        obj.madmData.behavioralEvents.buttonPress, obj.madmData.behavioralEvents.stimulusOnset,...
        obj.subjID, obj.paths.dataPath, obj.paths.figurePath, obj.paths.resultsPath, obj.paths.statsPath);
end

function obj = predictModelsIndNorm(obj)
    obj.predModelsIndNorm = predictModelsInd4Opt4AttNorm(obj.madmData.trialParams.TOD ,obj.madmData.eyeMovementData.fixDuration,...
        obj.madmData.eyeMovementData.fixAttributeMagnitude, obj.madmData.eyeMovementData.whichAttribute,...
        obj.madmData.behavioralEvents.buttonPress, obj.madmData.behavioralEvents.stimulusOnset,...
        obj.subjID, obj.paths.dataPath, obj.paths.figurePath, obj.paths.resultsPath, obj.paths.statsPath);
end

function obj = plotPredictModelsInd(obj)
    obj.plotPredModelsInd = plotPredictModelsInd4Opt4Att(obj.madmData.trialParams.TOD ,obj.madmData.eyeMovementData.fixDuration,...
        obj.madmData.eyeMovementData.fixAttributeMagnitude, obj.madmData.eyeMovementData.whichAttribute,...
        obj.madmData.behavioralEvents.buttonPress, obj.madmData.behavioralEvents.stimulusOnset,...
        obj.subjID, obj.paths.dataPath, obj.paths.figurePath, obj.paths.resultsPath);
end

function obj = plotPredictModelsIndNorm(obj)
    obj.plotPredModelsIndNorm = plotPredictModelsInd4Opt4AttNorm(obj.madmData.trialParams.TOD ,obj.madmData.eyeMovementData.fixDuration,...
        obj.madmData.eyeMovementData.fixAttributeMagnitude, obj.madmData.eyeMovementData.whichAttribute,...
        obj.madmData.behavioralEvents.buttonPress, obj.madmData.behavioralEvents.stimulusOnset,...
        obj.subjID, obj.paths.dataPath, obj.paths.figurePath, obj.paths.resultsPath);
end

function obj = prospectTheoryModelFit(obj)
    obj.ptModelFit = prospectTheoryModelFit4Opt4Att(obj.madmData.trialParams.infos_MAD_EPdata,...
        obj.subjID, obj.paths.dataPath, obj.paths.figurePath, obj.paths.resultsPath, 1);
end

function obj = cumulativeProportionsLongitudinal(obj)
    obj.cumProportionsLongitudinal = cumulativeProportionsLongitudinal4Opt4Att(obj.madmData.trialParams.TOD ,obj.madmData.eyeMovementData.fixDuration,...
        obj.madmData.eyeMovementData.fixAttributeMagnitude, obj.madmData.eyeMovementData.whichAttribute,...
        obj.madmData.behavioralEvents.buttonPress, obj.madmData.behavioralEvents.stimulusOnset,...
        obj.subjID, obj.paths.dataPath, obj.paths.figurePath, obj.paths.resultsPath, obj.paths.statsPath);
end

function obj = plotProbChooseOptAll(obj)
    obj.plotPChooseOptAll = plotProbChooseOptAll4Opt4Att(obj.subjID, obj.paths.dataPath, obj.paths.figurePath, obj.paths.resultsPath, obj.paths.statsPath);
end

function obj = predictModelsPop(obj)
    obj.predModelsPop = predictModelsPop4Opt4Att(obj.subjID, obj.paths.dataPath, obj.paths.figurePath, obj.paths.resultsPath, obj.paths.statsPath);
end

function obj = plotTransitionProps(obj)
    obj.plotTransProps = plotTransitionProps4Opt4Att(obj.subjID, obj.paths.statsPath, obj.paths.figurePath);
end

function obj = plotIndTrialTransitionProps(obj)
    obj.plotIndTrialTransProps = plotIndTrialTransitionProps4Opt4Att(obj.subjID, obj.paths.statsPath, obj.paths.figurePath);
end

function obj = plotIndTrialTransitionProps_wAttwAltType(obj)
    obj.plotIndTrialTransProps_wAttwAltType = plotIndTrialTransitionProps_wAttwAltType4Opt4Att(obj.subjID, obj.paths.figurePath, obj.paths.statsPath);
end

function obj = plotConsecutiveTransitionsAll(obj)
    obj.plotConsecTransitionsAll = plotConsecutiveTransitionsAll4Opt4Att(obj.subjID, obj.paths.figurePath, obj.paths.statsPath);
end

function obj = plotNumAttsFixed(obj)
    obj.plotNAttsFixed = plotNumAttsFixed4Opt4Att(obj.subjID, obj.paths.figurePath, obj.paths.statsPath);
end

function obj = plotRTToOptionTypeAllSubjects(obj)
    obj.plotRTToOptTypeAllSubjects = plotRTToOptionTypeAllSubjects4Opt4Att(obj.subjID, obj.paths.dataPath, obj.paths.figurePath, obj.paths.resultsPath, obj.paths.statsPath);
end

function obj = plotFixNumPropNonDomOptAttTypeHeatmapAllSubjects(obj)
    obj.plotFNumPropNonDomOptAttTypeHeatmapAllSubjects = plotFixNumPropNonDomOptAttTypeHeatmapAllSubjects4Opt4Att(obj.subjID, obj.paths.figurePath, obj.paths.statsPath);
end

function obj = plotFixNumAttTypePropAllSubjects(obj)
    obj.plotFNumAttTypePropAllSubjects = plotFixNumAttTypePropAllSubjects4Opt4Att(obj.subjID, obj.paths.dataPath, obj.paths.figurePath, obj.paths.resultsPath, obj.paths.statsPath);
end

function obj = fixNumRasters(obj)
    obj.fNumRasters = fixNumRasters4Opt4Att(obj.madmData.trialParams.TOD ,obj.madmData.eyeMovementData.fixDuration,...
        obj.madmData.eyeMovementData.fixAttributeMagnitude, obj.madmData.eyeMovementData.whichAttribute,...
        obj.madmData.behavioralEvents.buttonPress, obj.madmData.behavioralEvents.stimulusOnset,...
        obj.subjID, obj.paths.dataPath, obj.paths.figurePath, obj.paths.resultsPath, obj.paths.statsPath);
end

function obj = samplingStrategy(obj)
    obj.samplingStrat = samplingStrategy4Opt4Att(obj.madmData.trialParams.TOD ,obj.madmData.eyeMovementData.fixDuration,...
        obj.madmData.eyeMovementData.fixAttributeMagnitude, obj.madmData.eyeMovementData.whichAttribute,...
        obj.madmData.behavioralEvents.buttonPress, obj.madmData.behavioralEvents.stimulusOnset,...
        obj.subjID, obj.paths.dataPath, obj.paths.figurePath, obj.paths.resultsPath, obj.paths.statsPath);
end

function obj = plotNumUniqueOptsSampledAll(obj)
    obj.plotNUniqueOptsSampledAll = plotNumUniqueOptsSampledAll4Opt4Att(obj.subjID, obj.paths.figurePath, obj.paths.statsPath);
end

function obj = saccadeVigor(obj)
    obj.saccVigor = saccadeVigor4Opt4Att(obj.madmData.trialParams.TOD ,obj.madmData.eyeMovementData.fixDuration,...
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

function obj = pupilDiameter(obj)
    obj.pupilDiam = pupilDiameter4Opt4Att(obj.subjID, obj.paths.figurePath, obj.paths.statsPath, obj.paths.dataPath, obj.madmData.eyeMovementData.pupilData, 1);
end

    end
end