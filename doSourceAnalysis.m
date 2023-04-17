function sourceTimeCourse = doSourceAnalysis(s,sourceModel,outlierInfo,trialTime,varargin)
%function sourceAnalysisOutput = doSourceAnalysis(sourceModel,outlierInfo,trialTime,varargin)
%
% 1st:  Perform Source Localization and extract the time course of each Source.
% Last: Average all sources according to label (per atlas).
%
% Optional: Create 'general
% 
%
% Mandatory Inputs:
%
%           s               -   the sib object
%           outlierInfo  	-   information about outliers from egiPreprocess
%           sourceModel     -   the result from prepareSourceModel.m
%           triaTime        -   data to be extracted, relative to what 's'
%                               is aligned to (i.e. [0 1000])
%
%
%   AS - Aug 2020
%

    %check that the inputs are valid
        %is s a sib object
            if ~ischar(s.file)
                error("'s' is likely not a sib object.")
            end
    
        %does sourceModel contain all necessary fields
            essentialSourceModelFields = {'alignedElectrodesWithOutFid'; 'leadfield'; 'atlas'; 'sourceIdx2AtlasLabel'};
            for essentialSourceModelFieldsCntr = 1:length(essentialSourceModelFields)
                if ~isfield(sourceModel,essentialSourceModelFields{essentialSourceModelFieldsCntr})
                    error(['The sourceModel is not complete. The field ' essentialSourceModelFields{essentialSourceModelFieldsCntr} ' is missing!'])
                end
            end
    
        %is trialTime defined correctly
            if ~isnumeric(trialTime) || ~isequal(numel(trialTime),2)
                error('trialTime is not specified correctly. Provide [startTime stopTime].')
            end

            
    defaultSourceCombinationType = 'average';
    sourceCombinationTypeOptions = {'average'; 'best'; 'both'};
    
    %defaultSourcePowerType = 'none';
    %sourcePowerTypeOptions = {'none';'averageModel'}; 

    
    %parse all of the user's input parameters
    p = inputParser;
    p.StructExpand = true;
    p.KeepUnmatched = true;

     	addParameter(p,'byCondition', false,@islogical);    %source analysis separate by condition of for all conditions together
        addParameter(p,'stimulusTime', []);
        addParameter(p,'conditions', []);                   %identifyer for conditions if byCondition==true. If empty, then ft conditionIdentifyer will be used
       	addParameter(p,'nrSources', 1);                     %if by condition is true, then we can also keep multiple sources
        addParameter(p,'sourceCombinationType',defaultSourceCombinationType, @(x) any(validatestring(lower(x),sourceCombinationTypeOptions))); %average over sources or use the source with the highest power at a particular frequency
        addParameter(p,'sourceGain',false,@islogical);      %for steadyState paradigms
        addParameter(p,'targetFrequency',[]);               %if sourceCominationType is winner, then we can specify either a (vector of) frequency(/ies) 
                                                            %with one for all conditions (i.e.: [5]), one per condition (i.e.: [5 10 20 33]) 
                                                            %or a cell (array) with two values per cell to specify a range of
                                                            %frequencies over which to average {i.e: {[5 10]; [30 33]; [41 60], ...})
        %addParameter(p,'sourcePowerType',defaultSourcePowerType, @(x) any(validatestring(lower(x),sourcePowerTypeOptions)));%if targetFrequency is specified
                                                            %then we can save the response at that frequency either from the average or per trial, for each source
                                                            %requesting 'average' will allow localization of the peak response within an ROI

        
        
	parse(p,varargin{:});
    
    sourceTimeCourse.inputArgs = p.Results;
    
  	%does target frequency make sense in regard to sourceCombinationType chosen
        if strcmpi(p.Results.sourceCombinationType, 'best') || strcmpi(p.Results.sourceCombinationType, 'both')% || ~strcmpi(p.Results.sourcePowerType,'none')
            if isempty(p.Results.targetFrequency)
                error('If you choose to select a particular source or to locate the Peak within an ROI, then you need to specify a frequency as the criterion!')
            elseif  ~isequal(length(p.Results.targetFrequency),numel(unique(p.Results.conditions))) && ~isequal(length(p.Results.targetFrequency),1)
                error('Either specify one target frequency(-range) for all conditions or specify one for each condition!')
            end
        end

    
    %prepare the data for source analysis
        ftData = fieldtrip(s,'lfp',1:257,'time',trialTime);
        time = trialTime(1):diff(ftData.time{1}(1:2))*1000:trialTime(2);
        
        %replace fields of ftData with actual information
        ftData.label = sourceModel.alignedElectrodesWithOutFid.label';
        ftData.cfg.channel = ftData.label;
        
                     
        %find trials where we have the remaining electrodes
        useElectrodes = zeros([size(ftData.trial{1},1),size(ftData.trial,2)]);
        for trialCntr = 1:size(ftData.trial,2)
             useElectrodes((sum(isnan(ftData.trial{trialCntr}),2)==0),trialCntr) = 1;
        end

        
        if p.Results.byCondition
            if ~isempty(p.Results.conditions)
                conditions = p.Results.conditions;
            else
                conditions = ftData.trialinfo;
            end
        else
            conditions = ones(1,numel(ftData.trialinfo));
        end
        
        uCond = unique(conditions);
        
        if length(p.Results.targetFrequency) == 1
            targetFrequency = repmat(p.Results.targetFrequency,numel(uCond));
        else
            targetFrequency = p.Results.targetFrequency;
        end
        
        if ~iscell(targetFrequency)
            if size(targetFrequency,1)==1
                targetFrequency = targetFrequency';
            end
          	targetFrequency = mat2cell(targetFrequency,ones(1,size(targetFrequency,1)));
        end
        
        
        %create the struct will the source timecourse data beforehand, so
        %that we can make sure that the final result has trials atrranged
        %the same way, regardless of how many conditions we split the up in
        
        sourceTimeCourse.conditions = conditions;
        
        
        atlasNames = fieldnames(sourceModel.atlas);
        for atlasCntr = 1:length(atlasNames)
            if strcmpi(p.Results.sourceCombinationType, 'best') || strcmpi(p.Results.sourceCombinationType, 'both')
                sourceTimeCourse.(atlasNames{atlasCntr}).bestSource.data = nan([size(ftData.trial{1},2),numel(conditions),length(sourceModel.sourceIdx2AtlasLabel.(atlasNames{atlasCntr})),p.Results.nrSources]);
                %if only one source is requested, then remove that unnecessary 4th dimension of 1    
                sourceTimeCourse.(atlasNames{atlasCntr}).bestSource.data = squeeze(sourceTimeCourse.(atlasNames{atlasCntr}).bestSource.data);
            end
            if strcmpi(p.Results.sourceCombinationType, 'average') || strcmpi(p.Results.sourceCombinationType, 'both')
                sourceTimeCourse.(atlasNames{atlasCntr}).averageSource.data = nan([size(ftData.trial{1},2),numel(conditions),length(sourceModel.sourceIdx2AtlasLabel.(atlasNames{atlasCntr}))]);
          	end
        end
        if p.Results.sourceGain
            sourceTimeCourse.sourceGain = nan(length(ftData.trial),sum(sourceModel.leadfield.inside));
        end
        
        
        generalSource = cell(1,numel(uCond));
        
        for conditionCntr = 1:numel(uCond)
            
            trialLogIdx = conditions==uCond(conditionCntr);
            %trialIdx = find(trialLogIdx);
            ftCondData = ftData;
            targetFields = {'trialinfo';'sampleinfo';'trial';'time'};
            for targetFieldCntr = 1:length(targetFields)
                if strcmpi(targetFields{targetFieldCntr},'sampleinfo')
                    ftCondData.(targetFields{targetFieldCntr})(~trialLogIdx,:) = [];
                else
                    ftCondData.(targetFields{targetFieldCntr})(~trialLogIdx) = [];
                end
            end
        

            %create an average signal to create a spatial filter
            tempData = nan(size(ftCondData.trial{1},2),length(ftCondData.trial),size(ftCondData.trial{1},1));
            for trialCntr = 1:length(ftCondData.trial)
                for electrodeCntr = 1:size(ftCondData.trial{1},1)
                    tempData(:,trialCntr,electrodeCntr) = ftCondData.trial{trialCntr}(electrodeCntr,:);
                end
            end

            %create ftData object just for the average across all trials
            aveFtData = ftCondData;
            aveFtData.trialinfo = 1;
            aveFtData.trial = {(squeeze(nanmean(tempData,2)))'};
            aveFtData.time = aveFtData.time(1);
            aveFtData.sampleinfo = [1 numel(aveFtData.time{1})];
            aveFtData.label = aveFtData.label;

            clearvars tempData

            electrodeSelection= sourceModel.alignedElectrodesWithOutFid.label(logical(sum(useElectrodes,2)>0));

            %create timelocked data
            cfg                	= [];
            cfg.covariance     	= 'yes';
            cfg.channel         = electrodeSelection;
            timelock           	= ft_timelockanalysis(cfg, aveFtData);


        %source analysis on the average data to generate a filter that all
        %conditions then will share               
        
         	cfg                             = [];
            cfg.method                      = 'lcmv';
            cfg.channel                     = electrodeSelection';
            cfg.elec                        = sourceModel.alignedElectrodesWithOutFid;
            cfg.lcmv.normalize              = 'yes'; %lf being computed during this step, so need to specify this
            cfg.headmodel                   = sourceModel.headmodel; % volume conduction model (headmodel)
            cfg.grid                        = sourceModel.leadfield;
            cfg.keepleadfields              = 'yes';
            cfg.lcmv.keepfilter             = 'yes';
            cfg.lcmv.fixedori               = 'yes'; % project on axis of most variance using SVD
            cfg.lcmv.lambda                 = 5;
            cfg.lcmv.projectnoise           = 'yes';
            generalSource{conditionCntr}  	= ft_sourceanalysis(cfg, timelock);
        end
        
        if	strcmpi(p.Results.sourceCombinationType, 'best') || strcmpi(p.Results.sourceCombinationType, 'both')% || ~strcmpi(p.Results.sourcePowerType,'none')
            
            if isempty(p.Results.stimulusTime)
               	stimulusTime = trialTime;
            else
                stimulusTime = p.Results.stimulusTime;
            end
            signalLength =numel(stimulusTime(1):diff(ftData.time{1}(1:2))*1000:stimulusTime(2));
            
            fftSignalSource = cell(1,numel(uCond));
            fftFrequencies = s.iLFP.sf*(0:(signalLength/2))/signalLength;
            
            for atlasCntr = 1:length(atlasNames)
                for labelCntr = 1:length(sourceModel.sourceIdx2AtlasLabel.(atlasNames{atlasCntr}))
                    signalStrength = nan(numel(sourceModel.sourceIdx2AtlasLabel.(atlasNames{atlasCntr}){labelCntr}),numel(uCond));
                    if ~isempty(sourceModel.sourceIdx2AtlasLabel.(atlasNames{atlasCntr}){labelCntr})
                        
                        if numel(sourceModel.sourceIdx2AtlasLabel.(atlasNames{atlasCntr}){labelCntr})>1
                            for sourceCntr = 1:numel(sourceModel.sourceIdx2AtlasLabel.(atlasNames{atlasCntr}){labelCntr})
                                for conditionCntr = 1:numel(uCond)
                                    
%['Condition: ' num2str(conditionCntr) '/' num2str(numel(uCond)) ' - ' 'Source: ' num2str(sourceCntr) '/' num2str(numel(sourceModel.sourceIdx2AtlasLabel.(atlasNames{atlasCntr}){labelCntr})) ' - ' 'Label: ' num2str(labelCntr) '/' num2str(length(sourceModel.sourceIdx2AtlasLabel.(atlasNames{atlasCntr}))) ' - ' 'Atlas: ' num2str(atlasCntr) '/' num2str(length(atlasNames))]

                                    fftSignalTemp = fft(generalSource{conditionCntr}.avg.mom{sourceModel.sourceIdx2AtlasLabel.(atlasNames{atlasCntr}){labelCntr}(sourceCntr)}(time>stimulusTime(1) & time<=stimulusTime(2)));

                                    fftSignalTemp = (fftSignalTemp/signalLength);
                                    fftSignalTemp = fftSignalTemp(1:signalLength/2+1);
                                    fftSignalTemp(2:end-1) = 2*fftSignalTemp(2:end-1);

                                    fftSignalSource{conditionCntr}(:,sourceCntr) = abs(fftSignalTemp);

                                    if numel(targetFrequency{conditionCntr})>1
                                        targetFreq = knnsearch(fftFrequencies',targetFrequency{conditionCntr}(1)):knnsearch(fftFrequencies',targetFrequency{conditionCntr}(2));
                                    elseif ischar(targetFrequency{conditionCntr})
                                        targetFreq = NaN;
                                    else
                                        targetFreq = knnsearch(fftFrequencies',targetFrequency{conditionCntr}(1));
                                    end

                                    if isnan(targetFreq)
                                        signalStrength(sourceCntr,conditionCntr) = NaN;
                                    else
                                        signalStrength(sourceCntr,conditionCntr) = nanmean(fftSignalSource{conditionCntr}(targetFreq,sourceCntr),1);
                                    end

                                end

                            end

                            %find the best source for each condition
                            [powerVal, bestSourceIdx] = sort(signalStrength,'descend');

                            %replace nan conditions with the source that does best
                            %across all remaining conditions
                            nanCondition = find(sum(isnan(powerVal),1) == size(powerVal,1),size(powerVal,2));
                            if ~isempty(nanCondition)
                                if isequal(numel(nanCondition),size(powerVal,2))
                                    error('Looks like all target frequencies you spcified are characters. Change that to numbers and try again.')
                                end
                                %in this case we use the other conditions to choose a 'best' source

                                goodSources = bestSourceIdx;
                                goodSources(:,nanCondition) = [];
                                sourceScoreMatrix = reshape(repmat(1:size(goodSources,1),[1 size(goodSources,2)]),size(goodSources));

                                sourceOrder = nan(size(sourceScoreMatrix));
                                for conditionCntr = 1:size(goodSources,2)
                                    sourceOrder(goodSources(:,conditionCntr),conditionCntr) = sourceScoreMatrix(:,conditionCntr);
                                end

                                [~,replacementSourceOrder] = sort(sum(sourceOrder,2));

                                bestSourceIdx(:,nanCondition) = replacementSourceOrder;
                            end

                            sourceTimeCourse.(atlasNames{atlasCntr}).sourceIdxOrder{labelCntr} = bestSourceIdx;
                        else
                            sourceTimeCourse.(atlasNames{atlasCntr}).sourceIdxOrder{labelCntr} = ones(1,numel(uCond));
                        end
                    
                    else
                        sourceTimeCourse.(atlasNames{atlasCntr}).sourceIdxOrder{labelCntr} = NaN;
                    end
                        
                end
            end
        end     
                
             
        
        sourceCalcIdx =find(sourceModel.leadfield.inside);

        for conditionCntr = 1:numel(uCond) 

            trialLogIdx = conditions==uCond(conditionCntr);
            trialIdx = find(trialLogIdx);
            ftCondData = ftData;
            targetFields = {'trialinfo';'sampleinfo';'trial';'time'};
            for targetFieldCntr = 1:length(targetFields)
                if strcmpi(targetFields{targetFieldCntr},'sampleinfo')
                    ftCondData.(targetFields{targetFieldCntr})(~trialLogIdx,:) = [];
                else
                    ftCondData.(targetFields{targetFieldCntr})(~trialLogIdx) = [];
                end
            end


            singleTrialData = aveFtData;
            for trialCntr = 1:length(ftCondData.trial)
['Trial: ' num2str(trialCntr) '/' num2str(length(ftCondData.trial)) ' - ' 'Condition: ' num2str(conditionCntr) '/' num2str(numel(uCond))]
                clearvars sourcePerTrialTemp
                singleTrialData.trial{1} = ftCondData.trial{trialCntr};
                singleTrialData.sampleinfo = [1 numel(singleTrialData.time{1})];

                %which electrodes to exclude?
                trialElectrodeSelection= sourceModel.alignedElectrodesWithOutFid.label(logical(useElectrodes(:,trialIdx(trialCntr))));
                %some of the electrodes have been rejected because of bad impedance.
                %we need to remove those for good because they have not been included in the average model
                useElectrodesForSelection = useElectrodes(~sum(useElectrodes,2)==0,:);
                
                
                %and adjust the leadfield and filter
                trialFilter = generalSource{conditionCntr}.avg.filter;
                trialLeadfield = sourceModel.leadfield.leadfield;
                insideSource = find(sourceModel.leadfield.inside);
                for sourceCntr = 1:numel(insideSource)
                    trialFilter{insideSource(sourceCntr)} = trialFilter{insideSource(sourceCntr)}(logical(useElectrodesForSelection(:,trialIdx(trialCntr))));
                    trialLeadfield{insideSource(sourceCntr)} = trialLeadfield{insideSource(sourceCntr)}(logical(useElectrodesForSelection(:,trialIdx(trialCntr))));
                end

                tic
                cfg                         = [];
                cfg.covariance              = 'yes';
                cfg.channel                 = trialElectrodeSelection';
                singleTrialData.cfg.channel = trialElectrodeSelection';
                timelockByTrial             = ft_timelockanalysis(cfg, singleTrialData);


                % create spatial filter using the lcmv beamformer
                cfg                     = [];
                cfg.method              = 'lcmv';
                cfg.channel             = trialElectrodeSelection';
                cfg.elec                = sourceModel.alignedElectrodesWithOutFid;
                cfg.lcmvnormalize       = 'yes'; %lf being computed during this step, so need to specify this
                cfg.headmodel           = sourceModel.headmodel; % volume conduction model (headmodel)
                %cfg.grid                = sourceModel.leadfield;
                %cfg.grid.filter         = generalSource{conditionCntr}.avg.filter;
                cfg.sourcemodel.pos     = sourceModel.leadfield.pos;
                cfg.sourcemodel.inside  = sourceModel.leadfield.inside;
                cfg.sourcemodel.label   = trialElectrodeSelection;
                cfg.sourcemodel.dim     = sourceModel.leadfield.dim;
                cfg.sourcemodel.filter  = trialFilter;
                cfg.sourcemodel.leadfield = trialLeadfield;
                cfg.lcmv.keepfilter     = 'no';
                cfg.lcmv.fixedori       = 'no'; % project on axis of most variance using SVD
                sourcePerTrialTemp   	= ft_sourceanalysis(cfg, timelockByTrial);



            	%average all the sources for each label
                for atlasCntr = 1:length(atlasNames)
                    for labelCntr = 1:length(sourceModel.sourceIdx2AtlasLabel.(atlasNames{atlasCntr}))
                        if ~isempty(sourceModel.sourceIdx2AtlasLabel.(atlasNames{atlasCntr}){labelCntr})

                            if strcmpi(p.Results.sourceCombinationType, 'best') || strcmpi(p.Results.sourceCombinationType, 'both')
                                
                                tempSourceCourse = nan(numel(sourcePerTrialTemp.time),p.Results.nrSources);
                                for sourceCntr = 1:p.Results.nrSources
                                    tempSourceCourse(:,sourceCntr) = sourcePerTrialTemp.avg.mom{sourceModel.sourceIdx2AtlasLabel.(atlasNames{atlasCntr}){labelCntr}(sourceTimeCourse.(atlasNames{atlasCntr}).sourceIdxOrder{labelCntr}(sourceCntr,conditionCntr))};
                                end
                                sourceTimeCourse.(atlasNames{atlasCntr}).bestSource.data(:,trialIdx(trialCntr),labelCntr,:) = tempSourceCourse;
                            end
                            if strcmpi(p.Results.sourceCombinationType, 'average') || strcmpi(p.Results.sourceCombinationType, 'both')

                                tempSourceCourse = nan(numel(sourcePerTrialTemp.time),length(sourceModel.sourceIdx2AtlasLabel.(atlasNames{atlasCntr}){labelCntr}));
                                for sourceCntr = 1:length(sourceModel.sourceIdx2AtlasLabel.(atlasNames{atlasCntr}){labelCntr})
                                    tempSourceCourse(:,sourceCntr) = sourcePerTrialTemp.avg.mom{sourceModel.sourceIdx2AtlasLabel.(atlasNames{atlasCntr}){labelCntr}(sourceCntr)};
                                end
                                sourceTimeCourse.(atlasNames{atlasCntr}).averageSource.data(:,trialIdx(trialCntr),labelCntr) = nanmean(tempSourceCourse,2);
                                
                            end

                        else
                            if strcmpi(p.Results.sourceCombinationType, 'best') || strcmpi(p.Results.sourceCombinationType, 'both')
                                
                                if p.Results.nrSources>1
                                    sourceTimeCourse.(atlasNames{atlasCntr}).bestSource.data(:,trialIdx(trialCntr),labelCntr,:) = nan(numel(sourcePerTrialTemp.time),p.Results.nrSources);
                                else
                                    sourceTimeCourse.(atlasNames{atlasCntr}).bestSource.data(:,trialIdx(trialCntr),labelCntr) = nan(numel(sourcePerTrialTemp.time),1);
                                end
                            end
                            if strcmpi(p.Results.sourceCombinationType, 'average') || strcmpi(p.Results.sourceCombinationType, 'both')
                                sourceTimeCourse.(atlasNames{atlasCntr}).averageSource.data(:,trialIdx(trialCntr),labelCntr) = nan(numel(sourcePerTrialTemp.time),1);
                            end
                                
                        end
                    end
                    sourceTimeCourse.(atlasNames{atlasCntr}).labels = sourceModel.atlas.(atlasNames{atlasCntr}).tissuelabel;
                end


                if p.Results.sourceGain
                    for sourceCntr = 1:numel(sourceCalcIdx)
                        sourceTimeCourse.sourceGain(trialIdx(trialCntr),sourceCntr) = generalSource{conditionCntr}.avg.mom{sourceCalcIdx(sourceCntr)}'\sourcePerTrialTemp.avg.mom{sourceCalcIdx(sourceCntr)}';
                    end
                end


            end
            sourceTimeCourse.avgSourceLoc{conditionCntr} = generalSource{conditionCntr}; 
        end
                    
                

end

       