function [s, outlierInfo] = egiPreprocess(s,varargin)
%function [s, avoidSegments] = egiPreprocess(s,varargin)
%
% - Filter, detrend and re-reference the raw signal for eeg analysis
% - Identify and replace bad electrodes (via impedances) 
%   and trials (via artifact detection) (No proper method for interpolation
%   identified (yet) to handle the job satisfactory, so repacing with NaNs
% - perform any other preprocessing step that fieldtrip's 'preproc.m' allows
%
% AS 2019
%
%INFO - INTERPOLATION/REPLACEMENT - INFO
%one of the most imporant questions is about what to do with bad channes / bad segments.
%
%The option to interpolate will be kept, but is (at the moment) not recommended.
%None of the inpterpolation methods that fieldtrip offers seem to do a good
%job. Instead, huge artifacts can become 'infectious' and might spread to other
%electrodes, rather than imporiving anything. In addition, for source
%localization that information should not be necessary, since 257
%electrodes offer redundancy and interpolation will improve source localization.
%The main reason to interpolate trials would be for the conveniance to be able to just
%analyze everything, without worrying about which trials/electrodes to exclude. 
%SO, unless somebody finds a methods that works very well to replace
%electrodes/trials with data instead of NaNs, we should not do it.


replacemenOptions = {'reject','interpolate'};
defaultReplacementMethod = 'reject';    %interpolation works very poorly, it's best to replace bad electrodes/trials with nans


%INFO FOR - INTERPOLATION/REPLACEMENT - INFO

bpfilttypeOptions = {'but';'firws';'fir';'firls'};
defaultBpfilttype = 'but';
bpfiltwintypeOptions = {'hamming';'hann';'blackman';'kaiser'};
defaultBpfiltwintype = 'hamming';
bpfilterOptions = {'no';'yes'};
defaultBpfilter = 'yes';
dftfilterOptions = {'no';'dft';'interpolation'};
defaultDftfilter = 'dft';
detrendOptions = {'yes'; 'no'};
defaultDetrend = 'yes';
demeanOptions = {'yes'; 'no'};
defaultDemean = 'yes';
inpterpolationMethodOptions = {'weighted'; 'average'; 'spline'; 'slap'; 'nan'};
defaultInterpolationMethod = 'spline';
rerefOptions = {'yes'; 'no'};
defaultReref = 'yes';
refmethodOptions = {'avg';'median';'bipolar'};
defaultRefmethod = 'avg';

  p = inputParser;
    p.StructExpand = true;
    p.KeepUnmatched = true;
       	addParameter(p,'screenElectrodes',false,@islogical); 	%use impedances to identify bad electrodes and keep that information to not use them for later steps in the analysis (currently the same as replace electrodes)
        addParameter(p,'replaceElectrodes',true,@islogical);    %use impedances to identify bad electrodes and replace them with NaNs or by interpolation (Default is NaN)
      	addParameter(p,'preplacemenMethod',defaultReplacementMethod, @(x) any(validatestring(lower(x),replacemenOptions)));
        addParameter(p,'screenTrials',false,@islogical);        
        addParameter(p,'replaceTrials',false,@islogical);
        addParameter(p,'maxPercBadElectrodes',30,@isnumeric); 	%the percentage of electrodes that can have artifacts for a specific trial, so that we interpolate it with good data and not just other rubbish
     	addParameter(p,'outlierThreshold',4);
        addParameter(p,'trialTime',[],@isnumeric);
      	addParameter(p,'alignEvent',s.info.egi.BTRLNS);
        addParameter(p,'impedanceThreshold',99,@isnumeric);
        addParameter(p,'zcheckFiles',-1,@isnumeric);                         %-1  default: use first and last zcheck file, assuming htere are only 2 valid ones)
                                                                             %0   use all available files
                                                                             %indeces (i.e. 1/ [1 5]/ [2 3 4]);
        addParameter(p,'badElectrodeInterpolationMethod',defaultInterpolationMethod, @(x) any(validatestring(lower(x),inpterpolationMethodOptions)));
        addParameter(p,'badElectrodeInterpolationLambda',[]);
      	addParameter(p,'badTrialInterpolationMethod',defaultInterpolationMethod, @(x) any(validatestring(lower(x),inpterpolationMethodOptions)));
        addParameter(p,'badTrialInterpolationLambda',[]);
        %these are specific cfg inputs that fieldtrip allows for the
        %preproc function. We anticipate users to request (some of) these 
        %and therefore want to keep these options separate from whatever else 
        %they might want to 'preproc.m' to do and therefore use the same
        %names as fieldtrip does;
        %Those additional user requests will be passed to preproc.m via
        %parameters caught by p.Unmatched
        %For additional input arguments, see 'preproc.m' in the fieldtrip toolbox
        addParameter(p,'bpfilter',defaultBpfilter, @(x) any(validatestring(lower(x),bpfilterOptions)));
        addParameter(p,'bpfreq',[0.3 100],@isnumeric);
       	addParameter(p,'bpfilttype',defaultBpfilttype, @(x) any(validatestring(lower(x),bpfilttypeOptions)));
        addParameter(p,'bpfiltwintype',defaultBpfiltwintype, @(x) any(validatestring(lower(x),bpfiltwintypeOptions)));
      	addParameter(p,'bpfiltord',2,@isnumeric);
        addParameter(p,'dftfilter',defaultDftfilter, @(x) any(validatestring(lower(x),dftfilterOptions)));
        addParameter(p,'dftfreq',60,@isnumeric);        
      	addParameter(p,'demean',defaultDemean, @(x) any(validatestring(lower(x),demeanOptions)));
        addParameter(p,'detrend',defaultDetrend, @(x) any(validatestring(lower(x),detrendOptions)));
        addParameter(p,'reref',defaultReref, @(x) any(validatestring(lower(x),rerefOptions)));
        addParameter(p,'refchannel','all');    %default is all, alterrnatively user can specify (a) specific electrode(s)
        addParameter(p,'refmethod',defaultRefmethod, @(x) any(validatestring(lower(x),refmethodOptions)));  
	parse(p,varargin{:});
    
        tic;
        
        %fieldtrip should be on the path
        ft_defaults
        
        
        %if we want to screen for artifacts in trials, then we need to know
        %the time of interest and the user needs to specify a different align Event
        if p.Results.replaceTrials && isempty(p.Results.trialTime)
        	error("Specify 'trialTime' to be able to screen trial relevant periods for artifacts.");
        end
        
        if isequal(p.Results.alignEvent,s.info.egi.BTRLNS) && p.Results.trialTime(1)<0
           warning("You chose 'BTRLNS' as the align event and specified at least one value in 'trialTime' to be smaller than 0. Not sure if that will work...");
        end
        
        if (p.Results.screenTrials || p.Results.replaceTrials) && strcmpi(p.Results.bpfilter,'no')
           	warning("Artifact detection will likely not work on signal that is not high-pass filtered, since trials will be more different from each other than any artifact within trials would be from a trial's mean.")
        end
    
     	%if trials got removed, then some fields of the egi data will not be accessible anymore
        %So we read the sibf file again
        if isempty(s.info.egi.mfffile.data)
            sInfo = sib(s.file);
            mffFilePath = fullfile(fileparts(s.file),sInfo.info.egi.mfffile.data);
            electrodePath = fullfile(fileparts(s.file),sInfo.info.egi.COORDINATESFILE.data);
          	zBeforePath = fullfile(fileparts(s.file),sInfo.info.egi.ZBEFOREFILE.data);
            zAfterPath = fullfile(fileparts(s.file),sInfo.info.egi.ZAFTERFILE.data);
            [zAllZPathTemp, ~, ~, ~]= retrieve(sInfo.info,'allwithdata',s.info.egi.zAllFiles);
            zAllZPath = cell(1,length(zAllZPathTemp{1}));
            for zFileCntr = 1:length(zAllZPathTemp{1})
                zAllZPath{zFileCntr} = fullfile(fileparts(s.file),zAllZPathTemp{1}{zFileCntr});
            end
            viFlaggedElectrodes = fullfile(fileparts(s.file),sInfo.info.egi.VIFLAGGEDELECTRODESFILE.data);
            
            clearvars sInfo
        else
          	mffFilePath = fullfile(fileparts(s.file),s.info.egi.mfffile.data);
            electrodePath = fullfile(fileparts(s.file),s.info.egi.COORDINATESFILE.data);
          	zBeforePath = fullfile(fileparts(s.file),s.info.egi.ZBEFOREFILE.data);
            zAfterPath = fullfile(fileparts(s.file),s.info.egi.ZAFTERFILE.data);
            [zAllZPathTemp, ~, ~, ~]= retrieve(s.info,'allwithdata',s.info.egi.zAllFiles);
            zAllZPath = cell(1,length(zAllZPathTemp{1}));
            for zFileCntr = 1:length(zAllZPathTemp{1})
                zAllZPath{zFileCntr} = fullfile(fileparts(s.file),zAllZPathTemp{1}{zFileCntr});
            end
            viFlaggedElectrodes = fullfile(fileparts(s.file),s.info.egi.VIFLAGGEDELECTRODESFILE.data);
        end

        %read electrode information, so that we know labels, etc.
        try
            electrodeLocations = ft_read_sens(electrodePath);
        catch
            fieldtripPath = fileparts(which('ft_defaults'));
            electrodePath = [fieldtripPath '\template\electrode\GSN-HydroCel-257.sfp'];
            electrodeLocations = ft_read_sens(electrodePath);
          	warning('No GPS file found. You should be hesitant about interpolating electrodes and trials if you do not have accurate electrode positions. BUT DEFINITELY DO NOT USE THIS DATASET FOR SOURCE LOCALIZATION!')
        end

        %remove the fiducials and their locations.
        %At this point information about fiducials will only confuse fieldtrip
        fiducialIdx = contains(electrodeLocations.label,{'FidNz';'FidT9';'FidT10'});
        targetFields = {'chanpos';'chantype';'chanunit';'elecpos';'label'};
        for targetFieldCntr = 1:length(targetFields)
           	electrodeLocations.(targetFields{targetFieldCntr})(fiducialIdx,:) = [];
        end
        
        %load the raw signal (if it is not loaded yet)
        if ~s.iLFP.isLoaded
             s=readlfp(s,'e',true); % Read all channels into the iLFP object
        end

        
        %the continuous signal is now in the format that we need for
        %filtering the signal (if requested)
        if strcmpi(p.Results.bpfilter,'yes') || ~strcmpi(p.Results.dftfilter,'no') || strcmpi(p.Results.demean,'yes') || strcmpi(p.Results.detrend,'yes')  %if either one of these got requested, then that's what we'll do first with the continuous signal
           
         	%get the continuous data
            [trialTimes, continuousSignal] = raw2Continuous(s);

            %perform preprocessing steps
            time = (1:size(continuousSignal,2))/s.iLFP.sf;    
            
            %step into folder where preproc.m lives
            currentPath = pwd;
            try
                fieldtripPath = fileparts(which('ft_defaults'));
                if contains(fieldtripPath,'/')
                    %we are on the cluster
                    fieldtripPrivatePath = [fieldtripPath '/private/'];
                else
                    %we are not on the cluster
                    fieldtripPrivatePath = [fieldtripPath '\private\'];
                end
                cd(fieldtripPrivatePath)
               
                
        %bandpass filtering:
                %check the following fields and put those into the cfg object
                clearvars cfg
                %first check if any bandpass filtering is requested
                checkFields = {'bpfilter'; 'bpfreq'; 'bpfilttype'; 'bpfiltwintype'; 'bpfiltord'};

                for checkFieldCntr = 1:length(checkFields)
                    %provide the user with information about what's happening
                    if strcmpi(checkFields{checkFieldCntr},'bpfilter') && strcmpi(p.Results.bpfilter,'yes')
                        sibwarn(['Bandpass filter between: ' num2str(p.Results.bpfreq(1)) ' and ' num2str(p.Results.bpfreq(2)) ' Hz']);
                        sibwarn(['Filter type is: ' p.Results.bpfilttype]);
                        sibwarn(['Filter window type is: ' p.Results.bpfiltwintype]);
                        sibwarn(['Filter order is: ' num2str(p.Results.bpfiltord)]);
                    elseif strcmpi(checkFields{checkFieldCntr},'bpfilter') && strcmpi(p.Results.bpfilter,'no')
                        sibwarn('No bandpass filtering requested');
                    end
                    try
                        cfg.(checkFields{checkFieldCntr}) = p.Results.(checkFields{checkFieldCntr});
                    catch

                    end
                end
            
                %then do the bandpass filtering, if requested
                if exist('cfg','var')
                 	[continuousSignal, ~, ~, ~] = preproc(continuousSignal, electrodeLocations.label, time, cfg);
                end
            
        %detrending/demeaning (1st pass):
                %fieldtrip recommends to demean the signal before and after line noise filtering
                clearvars cfg
                checkFields = {'demean'; 'detrend'};
              	for checkFieldCntr = 1:length(checkFields)
                  	if strcmpi(checkFields{checkFieldCntr},'demean')
                        sibwarn(['Demean Signal: ' num2str(p.Results.demean)]);
                    elseif strcmpi(checkFields{checkFieldCntr},'detrend')
                        sibwarn(['Detrend Signal: ' num2str(p.Results.detrend)]);
                    end
                    try
                        cfg.(checkFields{checkFieldCntr}) = p.Results.(checkFields{checkFieldCntr});
                    catch

                    end
                end
               
                %do the detrending and demeaning if requested
              	if exist('cfg','var')
                 	[continuousSignal, ~, ~, ~] = preproc(continuousSignal, electrodeLocations.label, time, cfg);
                end
                
                
        %Line Noise Filtering:
                %next thing is to check if people want to do line noise
                %filtering and what kind they want
                clearvars cfg
                checkFields = {'dftfilter'; 'dftfreq'};
              	for checkFieldCntr = 1:length(checkFields)
                    if strcmpi(checkFields{checkFieldCntr},'dftfilter') && ~strcmpi(p.Results.dftfilter,'no')
                        sibwarn(['Line noise filter at: ' num2str(p.Results.dftfreq) ' Hz']);
                        if strcmpi(p.Results.dftfilter,'dft')
                          	sibwarn('Line noise filtering will NOT be done via fieldtrip, but instead by designing our own filter and using dft.');
                        elseif strcmpi(p.Results.dftfilter,'neighbour')
                            sibwarn('Line noise filtering will be done via fieldtrip, by using Neighbour Interpolation.');
                        end
                    end
                    
                    try
                        cfg.(checkFields{checkFieldCntr}) = p.Results.(checkFields{checkFieldCntr});
                    catch

                    end
                end
                
               	%then do the line noise filtering, if requested
              	if exist('cfg','var')
                    
                    %the performance of the fieldtrip 'discrete fourier transform ('dft')' line noise filter is
                    %not the best, so if a dft filter is requested (default), then we are using our own filter
                    
                    if strcmpi(cfg.dftfilter,'dft')
                        %However, first we need to go back to where we were and
                        %remove fildtrip from the path for a moment

                        %go back to old path
                        cd(currentPath)
                        %remove fieldtrip from path
                        warning('off','all');
                        fieldtripRemovePath = fileparts(which('ft_defaults'));
                        rmpath(genpath(fieldtripRemovePath));
                        warning('on','all');

                        %create filter
                        for filterFreqCntr = 1:numel(cfg.dftfreq)   %it is possible that people want multiple frequencies to be filtered out

                            lineNoiseFilter = designfilt('bandstopiir','FilterOrder',p.Results.bpfiltord, ...
                           'HalfPowerFrequency1',cfg.dftfreq(filterFreqCntr)-1,'HalfPowerFrequency2',cfg.dftfreq(filterFreqCntr)+1, ...
                           'DesignMethod','butter','SampleRate',s.iLFP.sf);

                            for electrodeCntr = 1:size(continuousSignal,1)
                                tempSignal = filtfilt(lineNoiseFilter,double(continuousSignal(electrodeCntr,:)));
                                continuousSignal(electrodeCntr,:) = tempSignal;
                            end
                        end
                        %if we used matlab's build-in filtfilt function to do the
                        %line noise filtering, then we will have to add fieldtrip
                        %to the path again
                        
                    	addpath(fieldtripPath)
                        ft_defaults
                
                        %step into fieldtrip folder
                        try
                            cd(fieldtripPrivatePath)
                        catch
                            cd(currentFolder)
                        end

                %use fieltrip instead
                    elseif strcmpi(cfg.dftfilter,'interpolation')
                        %this means we will use fieldtrip also for the line
                        %noise filtering, but with the neighbour
                        %interpolation method
                        cfg.dftfilter = 'yes';
                        cfg.dftreplace    = 'neighbour';
                        
                        [continuousSignal, ~, ~, ~] = preproc(continuousSignal, electrodeLocations.label, time, cfg);
                    end
                end
                
        %Demeaning/detrending (2nd pass)
                %next step is demeaning and detrending, if requested
                clearvars cfg
                checkFields = {'demean'; 'detrend'};
              	for checkFieldCntr = 1:length(checkFields)
                  	if strcmpi(checkFields{checkFieldCntr},'demean')
                        sibwarn(['Demean Signal: ' num2str(p.Results.demean)]);
                    elseif strcmpi(checkFields{checkFieldCntr},'detrend')
                        sibwarn(['Detrend Signal: ' num2str(p.Results.detrend)]);
                    end
                    try
                        cfg.(checkFields{checkFieldCntr}) = p.Results.(checkFields{checkFieldCntr});
                    catch

                    end
                end

                %do the detrending and demeaning
                if exist('cfg','var')
                    [continuousSignal, ~, ~, ~] = preproc(continuousSignal, electrodeLocations.label, time, cfg);
                end
                    
            catch
                cd(currentPath)
                warning('It was not possible to step into the folder where preproc.m is located')
            end
          	%go back to previous folder
          	cd(currentPath)
            
            %overwrite raw data with the results
            s = setRaw(s,continuousSignal,trialTimes,'continuous',true);
            clearvars continuousSignal
        end
        
        %did the user request electrode impedances and/or to replace
        %electrodes beyond the tolarence (default: 99 kOhm)
        if p.Results.screenElectrodes || p.Results.replaceElectrodes
             
            if p.Results.zcheckFiles>=0
                
                if p.Results.zcheckFiles == 0
                    zfileIdx = 1:length(zAllZPath);
                else
                    zfileIdx = p.Results.zcheckFiles;
                end
                tempImpedanceFile = mff_importinfon(zAllZPath{zfileIdx(zFileCntr)},1);
                impedanceFileImpedances = nan(numel(tempImpedanceFile.impedance)-1,numel(zfileIdx));
                for zFileCntr = 1:numel(zfileIdx)
                    tempImpedanceFile = mff_importinfon(zAllZPath{zfileIdx(zFileCntr)},1);
                    impedanceFileImpedances(:,zFileCntr) = tempImpedanceFile.impedance(1:end-1);
                end
                    
            	badElectrodeIdx = find(sum(impedanceFileImpedances>p.Results.impedanceThreshold,2)>0);
                
            elseif isequal(p.Results.zcheckFiles,-1)
        
                %check if we can find impedance measurements. If we find none, we have a problem.
                %If we find at least one, then we can at least keep working on this data set, but should inform the user that there is an issue 
                try
                    preImpedance =  mff_importinfon(zBeforePath,1);
                    preImpedance = preImpedance.impedance(1:end-1); %the last one is none of the recording channels
                catch
                    preImpedance = NaN;
                end

                try
                    postImpedance = mff_importinfon(zAfterPath,1);
                    postImpedance = postImpedance.impedance(1:end-1); %the last one is none of the recording channels
                catch
                    postImpedance = NaN;
                end

                if isnan(sum(preImpedance)) && isnan(sum(postImpedance))
                    error('No impedance Measures found!')
                elseif isnan(sum(preImpedance))
                    preImpedance = postImpedance;
                    warning('No preImpedance measure found. PreImpedance will be approximated by using postImpedance instead.')
                elseif isnan(sum(postImpedance))
                    postImpedance = preImpedance;
                    warning('No postImpedance measure found. PostImpedance will be approximated by using preImpedance instead.')
                end

                badElectrodeIdx = find(preImpedance>p.Results.impedanceThreshold | postImpedance>p.Results.impedanceThreshold);
            end
        
            %check if there is a file that indicates electrodes that should
            %be excluded based on previous visual inspection of the gps
            if ~isempty(s.info.egi.VIFLAGGEDELECTRODESFILE.data)
                viFlaggedTable = readtable(viFlaggedElectrodes);
                try
                    viFLagged = viFlaggedTable.Var1;
               	catch
                    %first entry could be misidentified as column header if
                    %it was entered as string rather than number
                    if ~isfield(viFlaggedTable,'Var1')
                        viFLagged = [viFlaggedTable.(viFlaggedTable.Properties.VariableNames{1}); str2double(viFlaggedTable.Properties.VariableDescriptions{1})];
                    end
                end 
               	badElectrodeIdx = unique([badElectrodeIdx; viFLagged]);
            end
                

            %also look for electrodes that might be bad based on their
            %voltage values
            [data,~] =get(s.iLFP,'e',s.iLFP.electrodes,'time',p.Results.trialTime);
            
            variation = zeros(1,size(data,3));
          	for electrodeCntr = 1:size(data,3)
              	tempData = data(:,:,electrodeCntr);
              	lowestVals = min(tempData);
               	highestVals = max(tempData);
                 
              	minMax = diff([lowestVals; highestVals]);
               	variation(electrodeCntr) = iqr(minMax);
            end
            additionalBadElectrodeIdx = find(variation> iqr(variation)*10);
            
            try
                badElectrodeIdx = cat(2,badElectrodeIdx,additionalBadElectrodeIdx);
            catch
                badElectrodeIdx = cat(2,badElectrodeIdx',additionalBadElectrodeIdx);
            end
            
            badElectrodeIdx = unique(badElectrodeIdx);
            clearvars data
            

            badElectrodeLabels = num2str(badElectrodeIdx'); %the labels should be strings
            badElectrodeLabels = [repmat('E',[size(badElectrodeLabels,1),1]) badElectrodeLabels]; %the leading character should be 'E' for all electrode entries
            badElectrodeLabels = cellstr(badElectrodeLabels); %cell array, not character array, to allow for different lengths of labels
            badElectrodeLabels = strrep(badElectrodeLabels,' ',''); %remove spaces
       
            % now we know which electrodes to replace
            sibwarn(['Electrodes: ' num2str(badElectrodeIdx) '  - cross the impedance threshold of: ' num2str(p.Results.impedanceThreshold) ' kOhm - or have particularly poor signal'])
            
            %warn the user in case the the reference is among those electrodes
           	if any(badElectrodeIdx==257)
                warning('The Reference (CZ) has bad impedances... . This might have serious consequences for the data quality of all your electrodes since they were referenced to CZ during data acquisiton')
            end
            
            %read continuous data (in case that we did not request it yet)
          	[trialTimes, ~] = raw2Continuous(s);
            %remove CZ from the bad
            
            if p.Results.replaceElectrodes
                
                %what should we do with bad electrodes? (interpolate or
                %replace with nans
                
                switch p.Results.preplacemenMethod
                    case 'interpolate'
                        
                        %replace those electrodes with the average of its neighbors
                        sibwarn(['Replacing bad electrodes via interpolation. The interpolation method is: ' p.Results.badElectrodeInterpolationMethod])
                        if any(badElectrodeIdx==257)
                            warning("CZ will be removed from the set of bad electrodes, so that we don't interpolate it before re-referencing")
                            badElectrodeIdx(badElectrodeIdx==257) = [];
                        end
                        %first we need information about the electrode labels
                        keepFiducials = 0;
                        electrodeLocations = readElectrodeLocations(electrodePath,keepFiducials);
                        %for that we first have to find those neighbors
                        neighborInformation = findElectrodeNeighbors(mffFilePath,badElectrodeIdx);


                        s=align(s,s.info.egi.BTRLNS);
                        data = fieldtrip(s,'lfp',1:257);
                        data.elec = electrodeLocations;

                        %now we need will interpolate bad channels
                        cfg = [];
                        cfg.method    	= p.Results.badElectrodeInterpolationMethod;
                        cfg.lambda      = p.Results.badElectrodeInterpolationLambda;
                        cfg.badchannel	= badElectrodeLabels;
                        cfg.neighbours 	= neighborInformation;  %fieldtrio uses english spelling for neighbor
                        cfg.trials     	= 'all';% or a selection given as a 1xN vector (default = 'all')

                        [data] = ft_channelrepair(cfg, data);   %repair channels
                        s = setRaw(s,data.trial,trialTimes);        %exchange 'raw data' with 'repaired' data
                        clearvars data
                        
                    case 'reject' %replace bad data with nans
                        sibwarn('Replacing bad electrodes with NaNs')
                        
                        ft_warning off     	%fieldtrip will complain because of the nan's. But there is only fake preprocessing done to get a fieldtrip like data structure.
                                          	%and we therefore don't care about fieldtrip's warnings.
                        data = fieldtrip(s,'lfp',1:257);
                        ft_warning on       %turn warning on again

                        if ~isempty(find(badElectrodeIdx==257,1)) %do not remove the reference electrode
                            badElectrodeIdx(find(badElectrodeIdx==257,1)) = [];
                        end
                        
                        for badElectrodesCntr = 1:numel(badElectrodeIdx)
                            for trialCntr = 1:length(data.trial)
                                data.trial{trialCntr}(badElectrodeIdx(badElectrodesCntr),:) = NaN;
                            end
                        end
                        s = setRaw(s,data.trial,trialTimes);        %exchange 'raw data' with 'repaired' data
                end
            end

        end
        
        keepFiducials = 0;
        %screen trials for artificts and replace them via interpolation from neighboring electrodes
        if p.Results.screenTrials || p.Results.replaceTrials
            [s, ~, trialOutliers, trialsKept]= removeOutliers(s,p,trialTimes,'electrodePath',electrodePath,'mffFilePath',mffFilePath,'keepFiducials',keepFiducials,'outlierThreshold',p.Results.outlierThreshold);
        else
            trialsKept = ones(1,s.nrTrials);
        end

      	%this last step allows users to run one last round of preprocessing
        %with any parameter that preproc.m allows. Those parameters are
        %passed via p.Unmatched.
        %Perform those steps before re-referencing
        if ~isempty(fieldnames(p.Unmatched))
            %created confic struct for preproc.m
            sibwarn('You requested some additional preprocessing steps from preproc.m. Those are done now.');
            cfg = [];
            
            parmFields = fieldnames(p.Unmatched);
            for parmFieldCntr = 1:length(parmFields)
                cfg.(parmFields{parmFieldCntr}) = p.Unmatched.(parmFields{parmFieldCntr});
            end
            
            %load raw data as continuous data
            [trialTimes, continuousSignal] = raw2Continuous(s);
            
            %create vector with time information, since it might be needed
            time = (1:size(continuousSignal,2))/s.iLFP.sf;
          	%read electrode locations, since they might be needed
            electrodeLocations = readElectrodeLocations(electrodePath,0);
            
        	%step into folder where preproc.m lives
          	currentPath = pwd;
            try
              	fieldtripPath = fileparts(which('ft_defaults'));
                if contains(fieldtripPath,'/')
                    %we are on the cluster
                    fieldtripPrivatePath = [fieldtripPath '/private/'];
                else
                    %we are not on the cluster
                    fieldtripPrivatePath = [fieldtripPath '\private\'];
                end
                cd(fieldtripPrivatePath)
               
                [continuousSignal, ~, ~, ~] = preproc(continuousSignal, electrodeLocations.label, time, cfg);

                %go back to previous folder
                cd(currentPath)
            catch
                cd(currentPath)
                warning('It was not possible to step into the folder where preproc.m is located')
            end
        
            %overwrite raw data with the results
            s = setRaw(s,continuousSignal,trialTimes,'continuous',true);
            clearvars continuousSignal
        end
        
        
        %re-reference data to all electrodes
        if strcmpi(p.Results.reref,'yes')
         	sibwarn('Re-referencing EEG signal.');
            
            %load raw data as continuous data
            [trialTimes, continuousSignal] = raw2Continuous(s);
 
            if ~strcmpi( p.Results.refmethod,'avg')
                warning('Currently, the only re-referencing method implemented is "avg". Will use the average of all electrodes instead.')
            end
            
            %re-reference all electrodes to the grand average
            averageSignal = nanmean(continuousSignal,1);
            continuousSignal = continuousSignal-averageSignal;
            
            %assign those results to the sib Object
            s = setRaw(s,continuousSignal,trialTimes,'continuous',true);
            
            clearvars continuousSignal
        end
        

        %combine information about bad electrodes and bad trials into
        %'useSegments'
        
        try
            trialOutliers(:,badElectrodeIdx) = true;
            outlierInfo.outlierTrials = ~trialOutliers;
            outlierInfo.keptTrials = trialsKept;
        catch
            warning('No Information about unusable trials or electrodes was requested. Be careful using these data for your analysis!')
            outlierInfo = [];
        end
        
        
        thisTookXMinutes = toc;
        thisTookXMinutes = thisTookXMinutes/60;
        sibwarn(['Done preprocessing. This took ' num2str(thisTookXMinutes) ' minutes'])
end




function [trialTimes, continuousSignal] = raw2Continuous(s)

    %Stitch trials together
    s=align(s,s.info.egi.BTRLNS);   %align to begin trial events
    tempRaw = s.iLFP.raw;           %get the raw signal
    startIdx = 1;
    %keep start and stop indeces of trials, so that we can go back to raw
    %data the way sib expects once we performed all operations we want on
    %the continuous signal
    for trialCntr = 1:size(tempRaw,2)
        trialTimes.startTimes(trialCntr) = startIdx;
        trialTimes.stopTimes(trialCntr) = startIdx-1 +size(tempRaw{trialCntr},1);
        continuousSignal(:,trialTimes.startTimes(trialCntr):trialTimes.stopTimes(trialCntr)) = tempRaw{trialCntr}';
        startIdx = trialTimes.stopTimes(trialCntr)+1;
    end

end

function s = setRaw(s,inputSignal,trialTimes,varargin)

    p = inputParser;
       	addParameter(p,'continuous',false,@(x) islogical(x) && isequal(length(logical(x)),1));    %use impedances to identify bad electrodes and keep that information to not use them for later steps in the analysis
	parse(p,varargin{:});

    s=align(s,s.info.egi.BTRLNS);   %align to the begining of the trial again
    tempRaw= cell(1,size(trialTimes.startTimes,2));
    if p.Results.continuous
        for trialCntr = 1:size(trialTimes.startTimes,2)
            tempRaw{trialCntr} = inputSignal(:,trialTimes.startTimes(trialCntr):trialTimes.stopTimes(trialCntr))';
        end
    else
        for trialCntr = 1:size(trialTimes.startTimes,2)
            tempRaw{trialCntr} = inputSignal{trialCntr}';
        end
    end
    s.iLFP = set(s.iLFP,tempRaw,'e',1:257,'sf',s.iLFP.sf,'rawAlignEvent','BTRLNS','rawAlignSource','egi','rawAlignOffset',0,'replace',true);

end


function neighborInformation = findElectrodeNeighbors(mffFilePath,electrodeIdx)
%load information about electrodeneighbors and put that information into a
%fieldtrip compatible format

    if nargin<2
        electrodeIdx = 1:257;
    end
    
    neighborInformationTemp =mff_importsensorlayout(mffFilePath);
    neighborInformation(length(neighborInformationTemp.neighbors)).label = [];
    
    %but we want that information to resemble what ft expects
    for electrodeCntr = 1:length(neighborInformationTemp.neighbors)
        neighborInformation(electrodeCntr).label = ['E' num2str(neighborInformationTemp.neighbors(electrodeCntr).channelnumber)];
        tempNeighblabel = num2str(neighborInformationTemp.neighbors(electrodeCntr).neighbors');
        tempNeighblabel = [repmat('E',[size(tempNeighblabel,1),1]) tempNeighblabel]; %#ok
        tempNeighblabel = cellstr(tempNeighblabel);
        neighborInformation(electrodeCntr).neighblabel = strrep(tempNeighblabel,' ','');
    end

    neighborInformation = neighborInformation(electrodeIdx);
    
 	%replace 'E257' with 'Cz' for the label of electrode 257 and all neighbors
        for electrodeCntr = 1:length(neighborInformation)
            if strcmpi(neighborInformation(electrodeCntr).label,'E257')
                neighborInformation(electrodeCntr).label = 'Cz';
            end
            if any(strcmpi(neighborInformation(electrodeCntr).neighblabel,'E257'))
                neighborInformation(electrodeCntr).neighblabel{find(strcmpi(neighborInformation(electrodeCntr).neighblabel,'E257'),1)} = 'Cz';
            end
        end

end


function electrodeLocations = readElectrodeLocations(electrodePath,keepFiducials)
	%read electrodeLocations
  	try
        electrodeLocations = ft_read_sens(electrodePath);
    catch
     	electrodePath = 'L:\github\fieldtrip\template\electrode\GSN-HydroCel-257.sfp';
      	electrodeLocations = ft_read_sens(electrodePath);
        warning('No GPS file found. You should be hesitant about interpolating electrodes and trials if you do not have accurate electrode positions. BUT DEFINITELY DO NOT USE THIS DATASET FOR SOURCE LOCALIZATION!')
  	end
    

    if nargin <2 || ~keepFiducials
    	fiducialIdx = contains(electrodeLocations.label,{'FidNz';'FidT9';'FidT10'});
        targetFields = {'chanpos';'chantype';'chanunit';'elecpos';'label'};
        for targetFieldCntr = 1:length(targetFields)
           	electrodeLocations.(targetFields{targetFieldCntr})(fiducialIdx,:) = [];
          	if strcmpi(targetFields{targetFieldCntr},'tra') && isfield(electrodeLocations,'tra')
              	electrodeLocations.(targetFields{targetFieldCntr})(:,fiducialIdx) = [];
          	end
        end
    end
    
end 
    
function outliers = detectOutliers(inputSignal,threshold)


    outliers = nan(size(inputSignal,2),size(inputSignal,3));
    for electrodeCntr = 1:size(inputSignal,3)
        minMaxTrial =zeros(1,size(inputSignal,2));
        
        %find outliers
        workSet = squeeze(inputSignal(:,:,electrodeCntr));

        for trialCntr = 1:size(workSet,2)
            minMaxTrial(trialCntr) = diff([min(workSet(:,trialCntr)) max(workSet(:,trialCntr))]);
        end
    	lb = median(minMaxTrial)-threshold*iqr(minMaxTrial);
      	ub = median(minMaxTrial)+threshold*iqr(minMaxTrial);
        badTrialLB = minMaxTrial<lb;
        badTrialUB = minMaxTrial>ub;

        outliers(:,electrodeCntr) = (badTrialLB + badTrialUB) >0;
    end

end

function [cleanS, trialTimes, trialOutliers, trialsKept] = removeOutliers(s,mainInputParser,trialTimes,varargin)


	p = inputParser;
    p.StructExpand = false;
    p.KeepUnmatched = true;
       	addParameter(p,'electrodePath',[]);
        addParameter(p,'mffFilePath',[]); 
        addParameter(p,'keepFiducials',[]);
        addParameter(p,'outlierThreshold',10); %10 is the default for the second pass
 	parse(p,varargin{:});

%we might want to repeat this action at the end of the preprocessing, to
%with a very high threshold, to make sure that no bad trials escaped 
    if ~mainInputParser.Results.detrend
        warning('You chose not to bandpass filter the signal. However, it is unlikely that this method of artifact detection will work without high-pass filtering beforehand.')
    end

    %the easiest way to accomplish this is by loading data via 'get'
    s = align(s,mainInputParser.Results.alignEvent);
    [data,~] =get(s.iLFP,'e',s.iLFP.electrodes,'time',mainInputParser.Results.trialTime);

    %identify outliers for electrodes
    outliers = detectOutliers(data,p.Results.outlierThreshold);

    % check if there is a trial that every electrode is thrown off by,
    %if that is the case, then outlier detection in general might be thrown off
    badTrials = ((sum(outliers,2)/size(outliers,1))*100)>mainInputParser.Results.maxPercBadElectrodes;
    allTrialIdx = 1:size(data,2);
    trialsKept = allTrialIdx;
    
    if sum(badTrials)>0
        trialsKept(badTrials) = [];
        s=removetrials(s,badTrials);
        trialTimes.startTimes(badTrials) = [];
        trialTimes.stopTimes(badTrials) = [];
    end
    
    %now repeat outlier detection
    [data,~] =get(s.iLFP,'e',s.iLFP.electrodes,'time',mainInputParser.Results.trialTime);
    outliers = detectOutliers(data,p.Results.outlierThreshold);
    
    
    %if Cz is the reference, then every trial will be marked as an outlier. 
    %Fix that by saying that none of those trials are outliers
    if isequal(sum(outliers(:,257)),size(outliers,1))
        outliers(:,257) = 0;
    end

    
    %no we know the outliers and can replace bad trials
    if mainInputParser.Results.replaceTrials

        %If there are still trials during which a lot of electrodes suffer from artifacts, then we'll remove those trials.
        %The threshold determining our tolerance of what we are willing to deal with is defined via maxPercBadElectrodes
        badTrials = ((sum(outliers,2)/size(outliers,1))*100)>mainInputParser.Results.maxPercBadElectrodes;
        if sum(badTrials)>0
            trialsKept(badTrials) = [];

            sibwarn(['There were ' num2str(sum(badTrials)) ' trials during which more than ' num2str(mainInputParser.Results.maxPercBadElectrodes) ' % of Electrodes showed artifacts'])
            s=removetrials(s,badTrials);
            outliers(badTrials,:) = [];
            sibwarn('Those trials were removed from your data set')

            %also remove those trials from 'trialTImes' for when we
            %want to bring them over to the raw data again
            trialTimes.startTimes(badTrials) = [];
            trialTimes.stopTimes(badTrials) = [];
        end

        data = fieldtrip(s,'lfp',1:257);
        
        %replace by interpolation or with nans
      	switch mainInputParser.Results.preplacemenMethod
            
            case 'reject' %replace bad data with nans
                sibwarn('Replacing bad trials with NaNs')

                for trialCntr = 1:length(data.trial)
                    badElectrodes = find(outliers(trialCntr,:));
                    if ~isempty(badElectrodes)
                        for badElectrodesCntr = 1:numel(badElectrodes)
                            data.trial{trialCntr}(badElectrodes(badElectrodesCntr),:) = NaN;
                        end
                    end
                end
                
                %assign results to output
                trialOutliers = outliers;
                cleanS = setRaw(s,data.trial,trialTimes);        %exchange 'raw data' with 'repaired' data
                
                
            case 'interpolate'
        
                %get information on electrode locations and all neighbors
                electrodeLocations = readElectrodeLocations(p.Results.electrodePath,p.Results.keepFiducials);
                neighborInformation = findElectrodeNeighbors(p.Results.mffFilePath,1:257);
                data.elec = electrodeLocations;

                %try doing this parallel, but since it is likely that there
                %is not enough memory to run this on 8 cores, we'll
                %iteratively adjust the number until it works
                nrCores = feature('numcores'); %how many cores do we have
                if nrCores > 5      %more than 6 is often too much to handle anyway
                    nrCores = 5;
                end

                %determine which trials must be replaced
                outlierTrials = find(sum(outliers,2)>0);

                sibwarn(['Starting to interpolate data for ' num2str(numel(outlierTrials)) ' Trials. This will be done via parallel processing to save time, but may make information displyed confusing.'])
                succesful = 0;
                while succesful<1
                    try
                        for trialCntr = size(outliers,1):-1:1
                            poorlyInterpolatedTrial(trialCntr).trial = [];
                            poorlyInterpolatedTrial(trialCntr).warning = [];
                        end

                        tempData = cell(1,size(outliers,1));
                        tempWarning = cell(1,size(outliers,1));
                        lastWarning = cell(1,size(outliers,1));

                        parfor (trialCntr = 1:numel(outlierTrials),nrCores) %do as many electrodes as possibles, at the same time

                            warning(['Working on trial ' num2str(trialCntr) ' out of ' num2str(numel(outlierTrials)) ' trials with artifacts.'])
                            %which, if any, electrodes show artifacts during this trial?
                            replaceTrialsForElectrodes = find(outliers(outlierTrials(trialCntr),:)); %#ok

                            if ~isempty(replaceTrialsForElectrodes)
                                tempNeighborInformation = neighborInformation(replaceTrialsForElectrodes); %#ok
                                tempBadElectrodeLabels =  {neighborInformation(replaceTrialsForElectrodes).label}';

                                replaceTrial = false(1,size(outliers,1));
                                replaceTrial(outlierTrials(trialCntr)) = true;

                                %there is a chance that neighboring electrodes are also bad.
                                %That is one potential reason why the interpolation might be bad.
                                %In that case fieldtrip will display an error message.
                                %If that happens, then we'll safe the error message and flag that trial for removal.
                                %Otherwise, we'll jsut take it as good enough and accept the trial.


                                %now we need will interpolate bad channels
                                cfg = [];
                                cfg.method    	= mainInputParser.Results.badElectrodeInterpolationMethod; %#ok
                                cfg.lambda      = mainInputParser.Results.badElectrodeInterpolationLambda;
                                cfg.badchannel	= tempBadElectrodeLabels;
                                cfg.neighbours 	= tempNeighborInformation;  %fieldtrio uses english spelling for neighbor
                                cfg.trials     	= replaceTrial;% or a selection given as a 1xN vector (default = 'all')

                                [tempData{trialCntr}] = ft_channelrepair(cfg, data);   %repair channels

                                %try to keep track of potentially worrysome error
                                %messages, so that we know that we should exclude
                                %particular trials
                                try
                                    tempWarning{trialCntr} = ft_warning('last'); %#ok
                                    lastWarning{trialCntr} = tempWarning{trialCntr}.message; %#ok
                                catch

                                end

                                if contains(lastWarning{trialCntr},'No good channels found') || contains(lastWarning{trialCntr},'bad spherical fit')
                                    poorlyInterpolatedTrial(trialCntr).trial = true;
                                    poorlyInterpolatedTrial(trialCntr).warning = lastWarning{trialCntr};
                                else
                                    poorlyInterpolatedTrial(trialCntr).trial = false;
                                end

                            else
                                poorlyInterpolatedTrial(trialCntr).trial = 0;
                            end
                        end
                        succesful = 1;
                    catch
                        %a failed pool will not release all memory
                        %properly, so the best thing to do is to close it
                        badPool = gcp();
                        delete(badPool)
                        %reduce number of cores used
                        nrCores = nrCores - 1; %if the parfor loop failed (supposedly because we ran out of memory), then reduce the number cores used by 1
                        %...try again
                    end
                end


                %check if there are any trials that we should remove from
                %the data set alltogether

                removeTheseTrials = find([poorlyInterpolatedTrial.trial]);
                if ~isempty(removeTheseTrials)
                    tempData(outlierTrials(removeTheseTrials)) = [];
                    s=removetrials(s,outlierTrials(removeTheseTrials));
                end

                %replace the trials that we chose to interpolate
                replaceTheseTrials = find(~cellfun(@isempty,tempData));
                for replaceTheseTrialsCntr = 1:numel(replaceTheseTrials)
                    data.trial{outlierTrials(replaceTheseTrials(replaceTheseTrialsCntr))} = tempData{replaceTheseTrials(replaceTheseTrialsCntr)}.trial{1};
                end
                
                 %assign results to output
                trialOutliers = false(size(outliers));          %if interpolation would work, then there would be no outliers left
                cleanS = setRaw(s,data.trial,trialTimes);     	%exchange 'raw data' with 'repaired' data
        end
        
        
    end
end




