function [sourceModel, sourceModelFilePath]= prepareSourceModel(s,varargin)
%function [sourceModel, sourceModelFilePath] = prepareSourceModel(s,varargin)
%
%
% - 12 steps to prepare a struct that has all necessary information to
%   perform a sourceAnalysis and most steps are kep as well to allow the
%   user to (re-)start from any of those 12 steps
%
% - If sourceModel exists and is complete, this function will load that existing model
%  
% - If sourceModel is not complete, or if the user requests to redo a
%   particular step, then that step and all consecutive steps will be repeated
%   To assure that no downstream field is ever kept even though
%
%
%   AS - Nov 2019
%

    %general parameters that specifiy the type of sourceModel and will
    %specify the (save)name of the resulting file
    mriOptions = {'individualScan'; 'ftTemplateScan'};
    defaultMri = 'individualScan';
    
    electrodeNetOptions = {'gpsNet'; 'templateNet'};
    defaultNet = 'gpsNet';
    
   	headmodelTypeOptions = {'threeTissues'; 'fiveTissues'};
    defaultHeadmodelType = 'threeTissues';
    
    %paramters for creating the headmodel/sourceGrid
    defaultSourceGridType = 'headmodel';
    sourceGridTypeOptions = {'headmodel';'ftTemplateGrid';'powerGrid'};
    
    defaultWarpSourceGrid2MNI = 'yes';
    warpSourceGrid2MNIOptions = {'yes';'no'};
    
  	defaultWarpSourceGridNonLinear = 'no';
    warpSourceGridNonLinearOptions = {'yes';'no'};
    
    %parameters for creating a mniVersion of the participants brain
 	defaultSpmVersionForMNIWarp = 'spm12';
    spmVersionForMNIWarpOptions = {'spm12';'spm8'};
    %if we choose spm12 (default) then we also get to choose nonlinear (default: yes)
  	defaultWarpMRINonLinear = 'yes';
    warpMRINonLinearOptions = {'yes';'no'};
    
    
    %parse all of the user's input parameters
    p = inputParser;
    p.StructExpand = true;
    p.KeepUnmatched = false;
        %the most basic
        addParameter(p,'ask', true, @islogical); %if (it looks like) the model is not complete or if the user requests to redo steps, confirmation is needed so that no information is lost.
       	addParameter(p,'saveSteps', true,@islogical);
        addParameter(p,'justPrepare', false,@islogical); %just complete up to selecting fiducials so that we can continue on the cluster afterwards
        addParameter(p,'outlierInfo', []);  %the second (optional) output from egiPreprocess showing which electrodes/trials not to use (bad electrodes should not be part of the model)
        addParameter(p,'mriTemplate', defaultMri, @(x) any(validatestring(lower(x),mriOptions)));                         %use the individual mri or the template
        addParameter(p,'individiualMRIPath',[]);                                                                          %specify where that participant's mri can be found, since we do not have this information in the sibfile (yet?)
        addParameter(p,'electrodeNet', defaultNet, @(x) any(validatestring(lower(x),electrodeNetOptions)));               %use the individual electrode locations or the template   
        addParameter(p,'headmodelType',defaultHeadmodelType, @(x) any(validatestring(lower(x),headmodelTypeOptions)));    %threeTissues [default]: three tissue types modeled (mesh is triangles and headmodel is concentric spheres)
                                                                                                                          %fiveTissues           : five tissue types modeled (mesh is hexahedrals and headmodel is based on simBio)
        %parameters that create the sourcegrid (using the headmodel takes more time, but is advised)
        addParameter(p,'sourceGridType',defaultSourceGridType, @(x) any(validatestring(lower(x),sourceGridTypeOptions))); %use the headmodel [default] (which I recommend), fieldtrip's template or the power et al (2011) network map
        addParameter(p,'sourceGridResolution',3, @isnumeric);                                                             %resolution in mm. If a template is chosen then this number should be the resolution of the template,
                                                                                                                          %otherwise (headmodel) it should be the resolution desired.
        %parameters in case a template sourcegrid is chosen
      	addParameter(p,'warpSourcegrid2MNI',defaultWarpSourceGrid2MNI, @(x) any(validatestring(lower(x),warpSourceGrid2MNIOptions)));
        addParameter(p,'warpSourcegridNonLinear',defaultWarpSourceGridNonLinear, @(x) any(validatestring(lower(x),warpSourceGridNonLinearOptions)));
        
        %parameters to create/inspect the electrode alignment to the headmodel
        addParameter(p,'inspectElectrodeAlignment',false,@islogical);                                                     %after the model has been created we can check if we are happy with the electrode placement
       	addParameter(p,'projectElectrodes',false,@islogical);                                                             %project electrodes to nearest vertex on the mesh
       	addParameter(p,'interactiveElectrodeAlignment',false,@islogical);                                                 %best to use together with inspectElectrodeAlignment, but only locally, sice this will not work on the cluster
        
        %to create a region labels we are going to create an mni version of
        %the participants brain. The following parameters relate to that process
        addParameter(p,'spmVersion',defaultSpmVersionForMNIWarp, @(x) any(validatestring(lower(x),spmVersionForMNIWarpOptions)));
        addParameter(p,'warpMNINonLinear',defaultWarpMRINonLinear, @(x) any(validatestring(lower(x),warpMRINonLinearOptions)));
        
        %which brainAtlas should we use to crease labels for the sources/dipoles
        addParameter(p,'brainAtlasPath',[]);                                                                              %provide a path to the atlasTemplate, otherwise we will use the aal atlas
     	addParameter(p,'brainAtlasName',[]);                                                                              %provide a name of the atlas if you provide a path, otherwise it will be named 'unknownAtlas'
        
        %are there any measures that we want redone
        addParameter(p,'redoMri', false, @islogical);                       %step1
        addParameter(p,'redoAlignment', false, @islogical);                 %step2
        addParameter(p,'redoSegmenting', false, @islogical);                %step3
        addParameter(p,'redoMesh', false, @islogical);                      %step4
        addParameter(p,'redoHeadmodel', false, @islogical);                 %step5
        addParameter(p,'redoSourceGrid', false, @islogical);                %step6
      	addParameter(p,'redoElectrodeAlignment', false, @islogical);        %step7
        addParameter(p,'redoLeadfield', false, @islogical);                 %step8
        addParameter(p,'redoMniMRI', false, @islogical);                    %step9
        addParameter(p,'redoMniSourceCoordinates', false, @islogical);      %step10
        addParameter(p,'redoAtlas', false, @islogical);                     %step11
      	addParameter(p,'redoSourceIdx2AtlasLabel', false, @islogical);      %step12
        
        
 	parse(p,varargin{:});

    
    %first check if we do already have a sourceModelStruct (of this type) for this dataSet
        sourceFileFolder = fileparts(s.file);
        participantNr = s.subject;
        
        %what is the content of our target/source folder
      	sourceFolderContent = dir(sourceFileFolder);
        sourceFolderContent = {sourceFolderContent.name};
        
        %check if a file with the source model information exists
            %we might me on the cluster
            currentFunctionName = dbstack;
            currentFunctionName = currentFunctionName.name;
            
            if contains(which(currentFunctionName),'/')	%we are on the cluster
                onCluster = true;
            else
             	onCluster = false;                      %we are local
            end
            
            %change the path depending on whether we are on the cluster or are if we are local
            %(we could do the same for the case that we expect and individual structural scan, but did not encounter that, 
            %but since we should not really create a source model based on a template mri we will not do that.)
            sourceModelFilePathIntended = fullfile(sourceFileFolder,[participantNr '_sourceModel_' p.Results.mriTemplate '_' p.Results.headmodelType '_' p.Results.sourceGridType '_' p.Results.electrodeNet '.mat']);
           	sourceModelFilePathAlternative = fullfile(sourceFileFolder,[participantNr '_sourceModel_' p.Results.mriTemplate '_' p.Results.headmodelType '_' p.Results.sourceGridType '_' electrodeNetOptions{2} '.mat']);
             
            %there actually is no difference between slashes/backslashes
            %between local and cluster when it comes to where files are saved
%             if onCluster
%                 sourceModelFilePathIntended = strrep(sourceModelFilePathIntended,'/','\');  %we need to replace slashes with backslashes if we are on the cluster
%                 sourceModelFilePathAlternative = strrep(sourceModelFilePathAlternative,'/','\');  %we need to replace slashes with backslashes if we are on the cluster
%             end
            %finally check if a file with source model information exisst (wherever we are)
            [~, subjectSourceModelFilename, ~] =fileparts(sourceModelFilePathIntended);
            
         	sourceModelExists = sum(contains(sourceFolderContent,subjectSourceModelFilename));
            
            %it is possible the model was intended to use an gpsNet, but
            %failed doing so because no solved gps file exists
            if ~sourceModelExists
                %try looking for a model with the netTemplate, instead
                [~, subjectSourceModelFilename, ~] =fileparts(sourceModelFilePathAlternative);
                sourceModelExists = sum(contains(sourceFolderContent,subjectSourceModelFilename));
                if sourceModelExists %if it exists now, then we'll change the filePath
                   	m=input("No sourcemodel based on the subject's solved gps file exists, however, there is a model based on the standard electrode locations. Do you want to use the existing (Y) or create a new model (N)? Y/N [N]:",'s');

                    if strcmpi(m,'n')
                        disp('Creating new source model')
                        sourceModelFilePath = sourceModelFilePathIntended;  %stay with the intended path
                    elseif  strcmpi(m,'y')
                        disp('Using existing model based on the template electrode locations')
                        sourceModelFilePath = sourceModelFilePathAlternative;
                    end
                    
                else
                    sourceModelFilePath = sourceModelFilePathIntended;
                end
            else
                sourceModelFilePath = sourceModelFilePathIntended;  %stay with the intended path
            end
                
            
           	sourceModelFields =     {'mri'    ; 'aligned'      ; 'segmented'     ; 'mesh'    ; 'headmodel'    ; 'sourceGrid'    ; 'alignedElectrodes';      ...
                                     'leadfield';      'mniMRI';     'mniSourceCoordinates';    'atlas';     'sourceIdx2AtlasLabel'};
          	redoSourceModelFields = {'redoMri'; 'redoAlignment'; 'redoSegmenting'; 'redoMesh'; 'redoHeadmodel'; 'redoSourceGrid'; 'redoElectrodeAlignment'; ...
                                     'redoLeadfield'; 'redoMniMRI'; 'redoMniSourceCoordinates'; 'redoAtlas'; 'redoSourceIdx2AtlasLabel'};
            %if a source model exists, then we want to load it
            if sourceModelExists >= 1
                sourceModel = load(sourceModelFilePath);
                if isfield(sourceModel,'sourceModel')
                    sourceModel = sourceModel.sourceModel;
                end
                
                %in addition, if it does exist, the we want to check if it is complete
                doSourceModelField = zeros(1,length(sourceModelFields));
              	for sourceModelFieldCntr = 1:length(sourceModelFields)
                  	doSourceModelField(sourceModelFieldCntr) = ~isfield(sourceModel,sourceModelFields{sourceModelFieldCntr}) | p.Results.(redoSourceModelFields{sourceModelFieldCntr});
                end
                %if we request any step to be redone, then all the
                %following steps should be redone as well
                if sum(doSourceModelField)>0
                    doSourceModelField(find(doSourceModelField,1):end) = 1;
                  	%remove all analysis fields from the one on that we decided
                    %to (re-)do, so that if something fails, we don't end up
                    %with previous parts of the source model steps
                    doStartIdx = find(doSourceModelField,1);
                    for removeCntr = doStartIdx:length(sourceModelFields)
                        if isfield(sourceModel,sourceModelFields{removeCntr})
                            sourceModel = rmfield(sourceModel,sourceModelFields{removeCntr});
                            %if we are working on the electrode net then we
                            %also have to make sure that we are removing
                            %the version without the fiducials, so that we
                            %do not end up with an older version in that field
                            if strcmpi(sourceModelFields{removeCntr},'alignedElectrodes')
                            	sourceModel = rmfield(sourceModel,'alignedElectrodesWithOutFid');
                            end
                        end
                    end
                end
            else
                doSourceModelField(1:length(sourceModelFields)) = 1;
            end
          

    	%make sure that the user indeed wants to create a new source model, so that we do not unintentionally overwrite an existing one
        if ~isequal(sum(doSourceModelField),0)
            
            %create a safety net to prevent users from accidentially overwriting 
            %existing source models, but at the same time make it possible 
            %to run this where supervision and inputs cannot be provided
            if ~p.Results.ask
                m = 'y';
            else
                targetField = sourceModelFields{find(doSourceModelField,1)};
                
                if isequal(find(doSourceModelField,1),1)
                    m=input("Either no sourcemodel could be found or you requested to recalculate an entirely new one. \n If you proceed any existing sourcemodel with those parameters will be overwritten. \n Do you want to proceed? (Y) or abort (N)? Y/N [N]:",'s');
                else
                 	m=input(['You requested to redo the model starting with the step: ' targetField '\n That step and all following steps will be overwritten. \n Do you want to proceed? (Y) or abort (N)? Y/N [N]:'],'s');
                end
            end
            
            
%% Start creating/re-doing the source model
            if strcmpi(m,'y')
           
                %if trials got removed, then some fields of the egi data will not be
                %accessible anymore
                if isempty(s.info.egi.mfffile.data)
                    s = sib(s.file);
                end


    %% Step 1:  Load the MRI
                if doSourceModelField(1)
                    
                    sibwarn('Loading MRI')

                    %preemptively set path to template, in case we need one
                    fieldtripPath = fileparts(which('ft_defaults'));

                    mriScanPath = [fieldtripPath '/external/spm8/templates/T1.nii'];
                    if onCluster
                        mriScanPath = strrep(mriScanPath,'\','/');
                    end

                    %check if the user wants the subjects scan and if it exists
                    switch lower(p.Results.mriTemplate)

                        case 'individualscan'

                            try
                                %mriPath = s.info.egi.mri; %this might be a field on day, but that info has to be gathered by hand at this point
                                mriPath = p.Results.individiualMRIPath;
                                sourceModel.mri = ft_read_mri(mriPath);
                            catch
                                warning("Failed loading the participant's structural scan. Using the fieldtrip template of a structural scan instead.")
                                warning("Be aware that you will have to change the input 'mriTemplate' to 'ftTemplateScan' or it will not be recognized that a model already exists when you are trying to work with the resulting source model!")
                                sourceModel.mri = ft_read_mri(mriScanPath);

                                %in this case we will have to change the name of the source model to reflect 
                                %the fact that we are actually not using an individual scan
                                mriTemplateUsed = mriOptions{2};
                                delete(sourceModelFilePath);    %first delete the previous model

                                %then rename the current one to show that we are not using an individual scan
                                sourceModelFilePath = fullfile(sourceFileFolder,[participantNr '_sourceModel_' mriTemplateUsed '_' p.Results.headmodelType '_' p.Results.sourceGridType '_' p.Results.electrodeNet '.mat']);
                            end

                        case 'fttemplatescan'   %use a participant's scan

                            sourceModel.mri = ft_read_mri(mriScanPath);
                    end

                    %overwrite 'sourceModel' with the new one, to update the mri
                    sourceModel.mri = ft_convert_units(sourceModel.mri,'mm'); %we should still be in mm but better safe than sorry
                    if p.Results.saveSteps
                        save(sourceModelFilePath,'sourceModel','-v7.3')
                    end
                end     



    %% Step 2: Align the MRI and make the user select ficucials
                if doSourceModelField(2)

                  	sibwarn('Aligning MRI')
                    
                    %The anatomical scan might be oriented weird, if so flip the orientation 
                    %(once is enough and twice doesn't do anything anyway)
                    fig = figure;
                    questionChoice = 0;
                    while questionChoice==0
                        midDim = round(size(sourceModel.mri.anatomy)/2);

                        subplot(2,2,1)
                            imagesc(imrotate(squeeze(sourceModel.mri.anatomy(:,midDim(2),:)),90));


                        subplot(2,2,2)
                            imagesc(imrotate(squeeze(sourceModel.mri.anatomy(midDim(1),:,:)),90));

                        subplot(2,2,3)
                            imagesc(imrotate(squeeze(sourceModel.mri.anatomy(:,:,midDim(3))),90));

                        colormap(gray)

                        m=input('Is the structural scan oriented correctly?, Y/N [N]:','s');

                        if strcmpi(m,'n')
                            disp('Flipping the anatomical image')
                            cfg = [];
                            cfg.method     =  'flip';
                            sourceModel.mri = ft_volumereslice(cfg,sourceModel.mri);
                        elseif  strcmpi(m,'y')
                            questionChoice = 1;
                            close(fig)
                        end
                    end

                    %select fiducals
                    cfg = [];
                    cfg.interactive='yes';
                    cfg.coordsys = 'ctf';
                    sourceModel.aligned = ft_volumerealign(cfg,sourceModel.mri);

                    %we could reslice the mri, and at the same time change the dimensions.
                    %However, I assume that the original mri does have isotropic voxels, so there should be no need for this.
%                     cfg            = [];
%                     cfg.resolution = 1;
%                     cfg.dim        = [256 256 256];
%                     sourceModel.aligned            = ft_volumereslice(cfg, sourceModel.aligned);
%                     

                    %overwrite 'sourceModel' with the new one, to update the aligned mri
                    if p.Results.saveSteps
                        save(sourceModelFilePath,'sourceModel','-v7.3')
                    end

                end
                
                if p.Results.justPrepare
                    warning('Succesfully completed placement of fiducials. Finish on cluster')
                    return
                end



    %% Step 3: Segment the MRI
                if doSourceModelField(3)
                    
                  	sibwarn('Segmenting MRI')                    
                    
                    switch lower(p.Results.headmodelType)
                        case 'threetissues'
                            cfg           = [];
                            cfg.output = {'gray','skull','scalp'};
                            cfg.brainsmooth             = 2;
                            cfg.skullsmooth             = 2;
                            cfg.scalpsmooth             = 2;
                            sourceModel.segmented  = ft_volumesegment(cfg, sourceModel.aligned);

                        case 'fivetissues'
                            cfg           = [];
                            cfg.brainsmooth             = 2;
                            cfg.skullsmooth             = 2;
                            cfg.scalpsmooth             = 2;
                            %need to segment into tissue types for FEM model
                            cfg.output = {'gray','white','csf','skull','scalp'};
                            sourceModel.segmented  = ft_volumesegment(cfg, sourceModel.aligned);
                    end

                    %overwrite 'sourceModel' with the new one, to update the segmented mri
                    if p.Results.saveSteps
                        save(sourceModelFilePath,'sourceModel','-v7.3')
                    end           
                end



    %% Step 4: Create the MESH
                if doSourceModelField(4)

                  	sibwarn('Creating Tissue Mesh')                    
                    
                    %create a mesh fo the number tissue types we requested via headmodelType
                    switch lower(p.Results.headmodelType)
                        case 'threetissues'
                            cfg=[];
                            cfg.tissue= {'gray','skull','scalp'};
                            cfg.numvertices = [3000 2000 1000];
                            sourceModel.mesh=ft_prepare_mesh(cfg,sourceModel.segmented);

                        case 'fivetissues'
                            cfg = [];
                            cfg.shift = 0.3; %**range 0-0.49, shift parameter controlling boundaries between tissue types; default = 0.3
                            cfg.method = 'hexahedral';
                            sourceModel.mesh = ft_prepare_mesh(cfg,sourceModel.segmented);
                    end

                    %overwrite 'sourceModel' with the new one, to update the segmented mri
                    if p.Results.saveSteps
                        save(sourceModelFilePath,'sourceModel','-v7.3')
                    end
                end



    %% Step 5: Create HEADMODEL    
    
                %if we are local and want to do the three tissue model, then we
                %can proceed either way.
                %If we want to use the five tissue model, then we have to be on
                %the cluster, since we can only use simBio on the cluster
                if doSourceModelField(5) && strcmpi(p.Results.headmodelType,'threeTissues')
                    
                    sibwarn('Creating Headmodel')

                    %create headmodel
                    cfg = [];
                    cfg.method='concentricspheres';
                    cfg.conductivity = [0.33 0.01 0.43];
                    sourceModel.headmodel = ft_prepare_headmodel(cfg, sourceModel.mesh);
                    sourceModel.headmodel = ft_convert_units(sourceModel.headmodel,'mm'); %we should still be in mm but better safe than sorry

                    %overwrite 'sourceModel' with the new one, to update the headmodel
                    if p.Results.saveSteps
                        save(sourceModelFilePath,'sourceModel','-v7.3')
                    end 

                elseif doSourceModelField(5) && strcmpi(p.Results.headmodelType,'fiveTissues') && onCluster

                    cfg = [];
                    cfg.method ='simbio'; %FEM method
                    cfg.conductivity = [0.33 0.14 1.79 0.01 0.43];   % order follows mesh.tissuelabel
                    sourceModel.headmodel = ft_prepare_headmodel(cfg, sourceModel.mesh);
                    sourceModel.headmodel = ft_convert_units(sourceModel.headmodel,'mm'); %we should still be in mm but better safe than sorry

                    %overwrite 'sourceModel' with the new one, to update the headmodel
                    if p.Results.saveSteps
                        save(sourceModelFilePath,'sourceModel','-v7.3')
                    end

                elseif doSourceModelField(5) && strcmpi(p.Results.headmodelType,'fiveTissues') && ~onCluster
                    warning("Everything so far seems to have gone well. You should proceed on the Cluster by calling thise function there if you want to use a five tissue type model for this participant, though.")
                    doSourceModelField(6:end) = 0;
                end



    %% Step 6: Create SourceGrid                
                if doSourceModelField(6)

                  	sibwarn('Creating Source Grid')                    
                    
                    switch lower(p.Results.sourceGridType)

                        case 'fttemplategrid'
                            %we are going to use the 4mm fieldtrip template for as the source grid

                            %where is fieldtrip
                            fieldtripPath = fileparts(which('ft_defaults'));
                            sourceGridTemplatePath = [fieldtripPath '\template\sourcemodel\standard_sourcemodel3d4mm.mat'];

                            if onCluster
                                sourceGridTemplatePath = strrep(sourceGridTemplatePath,'\','/');
                            end

                            templateGrid = load(sourceGridTemplatePath);


                        case 'powergrid'

                            %we currently do not know what merrit the power map
                            %has, so once we know if it is worth is and where
                            %we are going to put the template, this piece can
                            %easily be adjusted, but until then this is just a
                            %placeholder
                            %and will actually just use the fieldtrip template
                            %insted

                            fieldtripPath = fileparts(which('ft_defaults'));
                            sourceGridTemplatePath = [fieldtripPath '\template\sourcemodel\standard_sourcemodel3d4mm.mat'];

                            if onCluster
                                sourceGridTemplatePath = strrep(sourceGridTemplatePath,'\','/');
                            end

                            templateGrid = load(sourceGridTemplatePath);

                            warning('Currently the use of the power-grid is not implemented... Using fieldtrip template, instead!')
                            
                        otherwise
                            %nothign to be done, because in this case we
                            %are using the headmodel instead

                    end
                    
                    cfg = [];
                    
                    if strcmpi(p.Results.sourceGridType,'fttemplategrid') || strcmpi(p.Results.sourceGridType,'powergrid')

                        templateGrid = templateGrid.sourcemodel;
                        cfg.warpmni = p.Results.warpSourcegrid2MNI;
                        cfg.nonlinear = p.Results.warpSourcegridNonLinear;
                        templateGrid = ft_convert_units(templateGrid,'mm'); %we should still be in mm but better safe than sorry
                        cfg = [];
                        cfg.template = templateGrid; %Note that even though template_grid is in mm, this step yields a grid in cm (that is the same as if we had converted units of template_grid to cm before)
                        cfg.mri = sourceModel.aligned;
                        
                        
                    elseif strcmpi(p.Results.sourceGridType,'headmodel')
                        
                        cfg.headmodel = sourceModel.headmodel; 
                       
                    else
                       %there should be no way that we end up here 
                        
                    end
                    
                    %warpmni and nonlinear should not matter if we choose headmodel 
                    
                    cfg.resolution = p.Results.sourceGridResolution;
                    cfg.inwardshift = 1;
                    cfg.unit = 'mm';

                    sourceModel.sourceGrid = ft_prepare_sourcemodel(cfg);

                    %overwrite 'sourceModel' with the new one, to update the sourcegrid
                    if p.Results.saveSteps
                        save(sourceModelFilePath,'sourceModel','-v7.3')
                    end

                end



    %% Step 7: Align the electrode net to the aligned MRI
                if doSourceModelField(7)
                    
                  	sibwarn('Aligning Electrode Net to MRI')                     

                    %read electrode locations
                        %do we have a solved gps or does the user intend to use a template for electrode locations
                        if strcmpi(p.Results.electrodeNet,'templateNet') || isempty(s.info.egi.COORDINATESFILE.data)

                            fieldtripPath = fileparts(which('ft_defaults'));
                            electrodeTemplatePath = [fieldtripPath '\template\sourcemodel\standard_sourcemodel3d4mm.mat'];

                            if onCluster
                                electrodeTemplatePath = strrep(electrodeTemplatePath,'\','/');
                            end

                            if isempty(s.info.egi.COORDINATESFILE.data) && strcmpi(p.Results.electrodeNet,'gpsNet')
                                warning("You specified to use the individual's individual electrode locations, but none were found. Will use standard electrode locations instead and will rename source model file to reflect that!")

                                sourceModelFilePathAlternative = fullfile(sourceFileFolder,[participantNr '_sourceModel_' p.Results.mriTemplate '_' p.Results.headmodelType '_' p.Results.sourceGridType '_' electrodeNetOptions{2} '.mat']);

                                if onCluster
                                    sourceModelFilePathAlternative = strrep(sourceModelFilePathAlternative,'/','\');  %we need to replace slashes with backslashes if we are on the cluster
                                end

                                delete(sourceModelFilePath);%remove wrong file(name) and replace sourceModel(name) with updated one

                                sourceModelFilePath = sourceModelFilePathAlternative; %replace sourceModel name with new name
                                if p.Results.saveSteps
                                    save(sourceModelFilePath,'sourceModel','-v7.3')
                                end
                            end

                        else
                            electrodeTemplatePath = fullfile(fileparts(s.file), s.info.egi.COORDINATESFILE.data);
                        end

                        electrodeLocations = ft_read_sens(electrodeTemplatePath);
                        electrodeLocations = ft_convert_units(electrodeLocations,'mm'); % should be the same unit as MRI

%                       Not removing 'bad electrodes' anymore. Instead the
%                       user can exclude them in the source analysis
%                     %remove all electrodes that we can identify as 'bad electrodes'
%                     if ~isempty(p.Results.outlierInfo)
%                         badElectrodeIdx = find(sum(p.Results.outlierInfo.outlierTrials==0,1)==size(p.Results.outlierInfo.outlierTrials,1));
%                         %add 'E' to the elxtrode index
%                         for electrodeCntr = 1:numel(badElectrodeIdx)
%                             badElectrodeLabel = ['E' num2str(badElectrodeIdx(electrodeCntr))];
%                             badElectrodeIdxTemp = strcmpi(electrodeLocations.label,badElectrodeLabel);
% 
%                             %remove that electrode from all fields in 'electrodeLocations'
%                             targetFields = {'chanpos';'chantype';'chanunit';'elecpos';'label'};
%                             for targetFieldCntr = 1:length(targetFields)
%                                 electrodeLocations.(targetFields{targetFieldCntr})(badElectrodeIdxTemp,:) = [];
%                             end
%                         end
%                     end


                    %change label names to what fieldtrip expects
                    electrodeLocations.label{strcmp(electrodeLocations.label, 'FidNz')} = 'nasion';
                    electrodeLocations.label{strcmp(electrodeLocations.label, 'FidT9')} = 'left';
                    electrodeLocations.label{strcmp(electrodeLocations.label, 'FidT10')} = 'right';
                    % these are expressed in the coordinate system of the electrode position capture device
    %                 Nas = electrodeLocations.chanpos(strcmp(electrodeLocations.label, 'nasion'),:);
    %                 Lpa = electrodeLocations.chanpos(strcmp(electrodeLocations.label, 'left'),:);
    %                 Rpa = electrodeLocations.chanpos(strcmp(electrodeLocations.label, 'right'),:);

                %co-register eeg electrode locations to structural scan    
                    vox2head = sourceModel.aligned.transform; % transformation matrix of individual MRI

                    % transform voxel indices to MRI head coordinates
                    vox_Nas = sourceModel.aligned.cfg.fiducial.nas;  % fiducials saved in mri structure
                    vox_Lpa = sourceModel.aligned.cfg.fiducial.lpa;
                    vox_Rpa = sourceModel.aligned.cfg.fiducial.rpa;
                    head_Nas = ft_warp_apply(vox2head, vox_Nas, 'homogenous'); % nasion
                    head_Lpa = ft_warp_apply(vox2head, vox_Lpa, 'homogenous'); % Left preauricular
                    head_Rpa = ft_warp_apply(vox2head, vox_Rpa, 'homogenous'); % Right preauricular

                    %put the results into a struct that fieldtrip expects
                    mriFiducials.chanpos = [
                      head_Nas
                      head_Lpa
                      head_Rpa
                      ];
                    mriFiducials.label = {'nasion', 'left', 'right'};
                    mriFiducials.unit  = 'mm';
                    mriFiducials.elecpos = mriFiducials.chanpos;    %make sure that the field elecpos exists, but it is the same as chanpos

                    % coregister the electrodes to the MRI using fiducials
                    cfg = [];
                    cfg.method   = 'fiducial';
                    cfg.target   = mriFiducials;
                    cfg.elec     = electrodeLocations;
                    cfg.fiducial = {'nasion', 'left', 'right'};

                    sourceModel.alignedElectrodes = ft_electroderealign(cfg);
                    
                    
                    %electrode should be projected to the surface of the headmodel.
                    %There are more ways, but this seems to do a decent job
                    
                    if p.Results.projectElectrodes
                           
                        cfg = [];
                        if strcmpi(p.Results.headmodelType,'threetissues')
                            cfg.headshape = sourceModel.mesh(3);
                        elseif strcmpi(p.Results.headmodelType,'fivetissues')
                            cfg.headshape = sourceModel.headmodel;    %maybe using the headmodel will work here?
                        end

                        cfg.method    = 'project';
                        cfg.elec      = sourceModel.alignedElectrodes;
                        cfg.fiducial  = {'nasion', 'left', 'right'};
                        tempReAligned = ft_electroderealign(cfg);
                    
                        if p.Results.inspectElectrodeAlignment
                            
                            switch p.Results.headmodelType
                                
                                case ' threetissues'

                                    mesh = sourceModel.mesh;
                                    
                                case 'fivetissues'
                                    
                                    mesh = sourceModel.headmodel;
                            end
                            
                            %plot the 
                            plotElectrodeAlignment(mesh,sourceModel.alignedElectrodes,p.Results.headmodelType,'preReAlignment')

                          	plotElectrodeAlignment(mesh,tempReAligned,p.Results.headmodelType,'postReAlignment')


                            m=input('Did Re-Alignment Improve the Electrode Locations?, Y/N [N]:','s');

                            if strcmpi(m,'n')
                                disp("Discarding Re-Aligned Electrode Locations. - Using the 'old' Electrode Locations, instead.")

                            elseif  strcmpi(m,'y')

                                %update electrode locations with re-aligned coordinates
                                sourceModel.alignedElectrodes = tempReAligned;
                            end
                            
                        else
                            %update electrode locations with re-aligned coordinates
                            sourceModel.alignedElectrodes = tempReAligned;
                            
                        end
                    end
                    

                    %electrode locations without fiducials
                    sourceModel.alignedElectrodesWithOutFid = sourceModel.alignedElectrodes;

                    fiducialIdx = contains(electrodeLocations.label,{'nasion';'left';'right'});
                    targetFields = {'chanpos';'chantype';'chanunit';'elecpos';'label';'tra'};
                    for targetFieldCntr = 1:length(targetFields)
                        if strcmpi(targetFields{targetFieldCntr},'tra') && isfield(sourceModel.alignedElectrodesWithOutFid,'tra')
                            sourceModel.alignedElectrodesWithOutFid.(targetFields{targetFieldCntr})(:,fiducialIdx) = [];
                            sourceModel.alignedElectrodesWithOutFid.(targetFields{targetFieldCntr})(fiducialIdx,:) = [];
                        elseif strcmpi(targetFields{targetFieldCntr},'tra') && ~isfield(sourceModel.alignedElectrodesWithOutFid,'tra')
                            %do nothing
                        else
                            sourceModel.alignedElectrodesWithOutFid.(targetFields{targetFieldCntr})(fiducialIdx,:) = [];
                        end
                    end

                    %overwrite 'sourceModel' with the new one, to update the
                    %aligned Electrodes as well as the electrode net without
                    %fiducials (necessary for creating the leadmodel)
                    if p.Results.saveSteps
                        save(sourceModelFilePath,'sourceModel','-v7.3')
                    end
                end
                
                
                                
    %% Step 8: Calculate the leadfield
                if doSourceModelField(8)
                    
                  	sibwarn('Creating Leadfield')                     

                    %make sure that the headmodel and the elecotrode net are ready to create the leadfield
                    [tempHeadmodel, tempElectrodes] = ft_prepare_vol_sens(sourceModel.headmodel, sourceModel.alignedElectrodesWithOutFid);
                    sourceModel.headmodel = tempHeadmodel;
                    sourceModel.alignedElectrodesWithOutFid  =tempElectrodes;
                    
                    clearvars tempHeadmodel tempElectrodes

                    %Create the leadfield. Most importantly, this step will
                    %also calculate a transfer matrix.
                    cfg               = [];
                    cfg.sourcemodel    = sourceModel.sourceGrid;            % source points
                    cfg.headmodel = sourceModel.headmodel;                  % volume conduction model
                    cfg.elec    = sourceModel.alignedElectrodesWithOutFid;  % electrode information
                    cfg.reducerank    = 3;
                    sourceModel.leadfield = ft_prepare_leadfield(cfg);
                    
                    
                	%overwrite 'sourceModel' with the new one, to update the
                    %aligned Electrodes as well as the headmodel and to add
                    %the leadfield (from now on saving each step will
                    %unfortunatly take much longer because of the size of
                    %the leadfield(/transfer matrix)
                    if p.Results.saveSteps
                        save(sourceModelFilePath,'sourceModel','-v7.3')
                    end
                    
                end
                
                
                
    %% Step 9: Create a MNI-space version of the MRI
                if doSourceModelField(9)

                  	sibwarn('Creating MNI Version of Individual Brain')                     
                    
                    %we need to know where the modeled sources for the participant
                    %are located in respect to mni space. Creating a MNI version 
                    %of the participants brain will give us transformation
                    %matrices (with linear and, if requested, nonliner
                    %components).
                    
                    if strcmpi(p.Results.spmVersion,'spm12') 
                        cfg = [];
                        cfg.spmversion = 'spm12';
                        cfg.nonlinear = p.Results.warpMNINonLinear;
                        cfg.spmmethod = 'new';
                        sourceModel.mniMRI = ft_volumenormalise(cfg, sourceModel.aligned);
            
                    elseif strcmpi(p.Results.spmVersion,'spm8')

                        sourceModel.mniMRI = ft_volumenormalise([], sourceModel.aligned);
                       
                    else
                        error('There should be no way the we end up here')
                    end
                    
                    
                 	%overwrite 'sourceModel' with the new one, to add the
                    %MNI space version of the participant's brain
                    if p.Results.saveSteps
                        save(sourceModelFilePath,'sourceModel','-v7.3')
                    end
                    
                end                
                
                
                
    %% Step 10: Calculate MNI-space coordinates of the sources (leadfield) 
                if doSourceModelField(10)
                    
                  	sibwarn('Creating MNI Version of Leadfield')                     
                    
                    %with trnasformation information from the mni-space
                    %version of the participant's brain its now possible to
                    %morph the source locations of the individual('s)
                    %leadfield to mni space
                	mniSourceCoordinates = ft_warp_apply(sourceModel.mniMRI.params,ft_warp_apply(sourceModel.mniMRI.initial, sourceModel.leadfield.pos), 'individual2sn');
                        
                	%keep only whats inside the brain
                    sourceModel.mniSourceCoordinates = mniSourceCoordinates(sourceModel.leadfield.inside,:);
                        
                    
                    %overwrite 'sourceModel' with the new one, to add the
                    %MNI space version of the source coordinates
                    if p.Results.saveSteps
                        save(sourceModelFilePath,'sourceModel','-v7.3')
                    end
                    
                end
                
                
                
    %% Step 11: Load a brain atlas (mni-based)
                if doSourceModelField(11)
                    
                  	sibwarn('Adding Brain Atlases to SourceModel')                    
                    
                    %The atlas provides coordinates for labels that we can match the source coordinates to.
                    % we are going to use one atlas (that the user can define)
                    % as the 'whole brain' atlas. And a second one
                    % specifically for a 'high resolution' of visual areas
                    
                   %find the atlas
                    fieldtripPath = fileparts(which('ft_defaults'));

                    if onCluster
                        brainAtlasDir = [fieldtripPath '/template/atlas/aal/ROI_MNI_V4.nii'];
                        visualAtlasDir =  [fieldtripPath '/template/atlas/vtpm/vtpm.mat'];

                    else
                        brainAtlasDir = [fieldtripPath '\template\atlas\aal\ROI_MNI_V4.nii'];
                        visualAtlasDir =  [fieldtripPath '\template\atlas\vtpm\vtpm.mat'];
                    end

                    %load the atlases
                    if ~isempty(p.Results.brainAtlasPath)
                        try
                            brainAtlasTemplate = ft_read_atlas(p.Results.brainAtlasPath);
                            if ~isempty(p.Results.brainAtlasName)
                                sourceModel.atlas.(p.Results.brainAtlasName) = brainAtlasTemplate;
                            else
                                warning("If you choose you own atlas, then you should specifiy the name as well. Atlas will be called 'unknownAtlas', instead")
                                sourceModel.atlas.unknownAtlas = brainAtlasTemplate;
                            end
                        catch
                            warning("The path to the brain atlas you provided is wrong/does not exist. Using the 'aal' atlas instead")
                            brainAtlasTemplate = ft_read_atlas(brainAtlasDir);
                            sourceModel.atlas.aal = brainAtlasTemplate;

                        end
                    else
                        brainAtlasTemplate = ft_read_atlas(brainAtlasDir);
                        sourceModel.atlas.aal = brainAtlasTemplate;
                    end

                    %also load the visual atlas
                    visualAtlasTemplate = load(visualAtlasDir);
                    sourceModel.atlas.vtpm = visualAtlasTemplate.vtpm;


                  	%overwrite 'sourceModel' with the new one, 
                    %to add the brain/visual atlas
                    if p.Results.saveSteps
                        save(sourceModelFilePath,'sourceModel','-v7.3')
                    end
                    
                end
                
                
                
    %% Step 12: Label sources according to the chosen atlas (mni-based)
                if doSourceModelField(12)
                    
                  	sibwarn("Creating Map from Individual's Sources to MNI-Space")                     
                    
                    %step into fieldtrips private folder so that we can use atlas_lookup
                    currentPath = pwd;

                    atlasNames = fieldnames(sourceModel.atlas);
                    
                    for atlasCntr = 1:length(atlasNames)
                        sourceAtlas.(atlasNames{atlasCntr}).label = cell(1,size(sourceModel.mniSourceCoordinates,1));
                    end
                    
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

                        for posIdx = 1:size(sourceModel.mniSourceCoordinates,1)
                            for atlasCntr = 1:length(atlasNames)
                                sourceAtlas.(atlasNames{atlasCntr}).label{posIdx} = atlas_lookup(sourceModel.atlas.(atlasNames{atlasCntr}), sourceModel.mniSourceCoordinates(posIdx,:),'inputcoord','mni');
                            end
                        end   
                        cd(currentPath)
                    catch
                        cd(currentPath)
                    end

                    
                    for atlasCntr = 1:length(atlasNames)

                        % Some coordinates might be assigned to multiple labels
                        %Assign each Coordinate to the label that it has been
                        containsLabels = ~cellfun(@isempty,sourceAtlas.(atlasNames{atlasCntr}).label);
                        labelIdx = find(containsLabels);
                        probabilisticLabels = cell(1,numel(labelIdx));

                        for labelIdxCntr = 1:numel(labelIdx)
                            if length(unique(sourceAtlas.(atlasNames{atlasCntr}).label{labelIdx(labelIdxCntr)}))>1
                               [whichLabel, howMany] = unique(sourceAtlas.(atlasNames{atlasCntr}).label{labelIdx(labelIdxCntr)});
                               [~, where] = max(howMany);
                               probabilisticLabels{labelIdxCntr} = whichLabel{where};
                            else
                                probabilisticLabels{labelIdxCntr} = sourceAtlas.(atlasNames{atlasCntr}).label{labelIdx(labelIdxCntr)}{1};
                            end
                        end

                        isInside = find(sourceModel.leadfield.inside);
                        insideSourceWithLabelIdx = isInside(labelIdx);

                        %the final label 
                      	for atlasLabelCntr = 1:length(sourceModel.atlas.(atlasNames{atlasCntr}).tissuelabel)
                            sourceModel.sourceIdx2AtlasLabel.(atlasNames{atlasCntr}){atlasLabelCntr} = insideSourceWithLabelIdx(cellfun(@contains,repmat(sourceModel.atlas.(atlasNames{atlasCntr}).tissuelabel(atlasLabelCntr),[1, length(probabilisticLabels)]),probabilisticLabels));
                        end
                    end
                    
                 	%overwrite 'sourceModel' with the new one, 
                    %to add mapping from source coordinates to tissuelabels
                    if p.Results.saveSteps
                        save(sourceModelFilePath,'sourceModel','-v7.3')
                    end 

                    
                end                



            else
                warning('You chose to abort recalculating/generating a new source model') %done with creating/updating the source model
            end
            
        elseif isequal(sum(doSourceModelField),0)
            warning('Done Loading Source Model')
        end
        
        if sum(doSourceModelField)>0 && strcmpi(m,'y')
            warning('Done Generating Source Model')
        end
        

end






%% helper functions

function plotElectrodeAlignment(mesh,electrodeLocations,headmodelType,titleName)


        switch lower(headmodelType)
            
            case 'threetissues'

                figure;
                    %plot brain mesh
                    subplot(2,2,1)
                        ft_plot_sens(electrodeLocations, 'elec', true, 'elecshape', 'circle', 'orientation', true, 'FaceColor','b');hold on;    %plot the electrode locations with orientation
                        ft_plot_mesh(mesh(1)); %plot the headmodel
                        title('brain')
                    %plot skull mesh
                    subplot(2,2,2)
                        ft_plot_sens(electrodeLocations, 'elec', true, 'elecshape', 'circle', 'orientation', true, 'FaceColor','b');hold on;    %plot the electrode locations with orientation
                        ft_plot_mesh(mesh(2)); %plot the headmodel
                        title('skull')
                    %plot skin mesh
                    subplot(2,2,3)
                        ft_plot_sens(electrodeLocations, 'elec', true, 'elecshape', 'circle', 'orientation', true, 'FaceColor','b');hold on;    %plot the electrode locations with orientation
                        ft_plot_mesh(mesh(3)); %plot the headmodel
                        title('skin')         
                    %plot all tissue types and electrodes together    
                    subplot(2,2,4)    
                        ft_plot_mesh(mesh(3), 'facecolor',[0.2 0.2 0.2], 'facealpha', 0.3, 'edgecolor', [1 1 1], 'edgealpha', 0.05);
                        hold on;
                        ft_plot_mesh(mesh(2),'edgecolor','none','facealpha',0.4);
                        hold on;
                        ft_plot_mesh(mesh(1),'edgecolor','none','facecolor',[0.4 0.6 0.4]);
                        ft_plot_sens(electrodeLocations, 'elec', true, 'elecshape', 'circle', 'orientation', true, 'FaceColor','b');hold on;    %plot the electrode locations with orientation
                        
                if ~isempty(titleName)
                    suptitle(titleName)
                end
                
                
            case 'fivetissues'
                
                markerIdx = find(mesh.tissue==5);
                %'downsample' so that the figure stays rotatable
                downSampledVerticesIdx = markerIdx(randperm(numel(markerIdx),1000000));
                
                
                figure;
                
                scatter3((mesh.pos(downSampledVerticesIdx,1)),(mesh.pos(downSampledVerticesIdx,2)),(mesh.pos(downSampledVerticesIdx,3)),'.','MarkerFaceAlpha',.1,'MarkerEdgeAlpha',0.1); hold on
                axis equal
                ft_plot_sens(sourceModel.alignedElectrodes, 'elec', true, 'elecshape', 'circle', 'orientation', true, 'FaceColor','b');hold on;

                if ~isempty(titleName)
                    title(titleName)
                end

        end
end
