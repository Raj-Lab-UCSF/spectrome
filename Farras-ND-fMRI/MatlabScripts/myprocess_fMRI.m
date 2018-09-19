function myprocess_fMRI(foldername)

% See CONN_fMRI_batch_manual.pdf for help.
% Note: The function assumes that the images are already preprocessed in
% SPM (realignment, coregistering, etc).
%
% For preprocessing you can use DPARSF
%
% F. Abdelnour
%
% **Modify** 
% ROILocation, 
% CSFMask,
% GreyMask,
% WhiteMask
%
% Notice:
% Since the structural and functional images are already preprocessed
% (realignment, coregistering, etc) we don't include BATCH.New.*.

% CSF/Grey/White masks locations...
CSFMask = '/home/farras/Documents/MATLAB/myToolboxes/spm8/tpm/csf.nii';
GreyMask = '/home/farras/Documents/MATLAB/myToolboxes/spm8/tpm/grey.nii';
WhiteMask = '/home/farras/Documents/MATLAB/myToolboxes/spm8/tpm/white.nii';

% Atlas file location...
ROIName = {'atlas116'};
ROILocation = '/home/farras/Documents/MATLAB/myToolboxes/spm8/conn/rois/atlas116.img';

BATCH.filename = [foldername filesep 'conn_process.mat'];
%BATCH.Setup.isnew = 1;

%BATCH.New.steps = {'realignment','coregistration','segmentation','normalization','smoothing','initialization'};

func_dir = dir([foldername filesep 'fMRI']);
func_vols = cell(1,0);
for i = 1:length(func_dir)
    if(~func_dir(i).isdir)
        [path,name,ext] = fileparts(func_dir(i).name);
        %Filter only original files - second letter is uppercase
        if(strcmp(ext,'.img') && strcmp(name(1:3),'000'))  % preprocessed images, was 'wrf' 
            func_vols = [func_vols; {[foldername filesep 'fMRI' filesep name ext]}];
        end
    end
end

BATCH.Setup.functionals = {{func_vols}};

struc_dir = dir([foldername filesep 'T1scan']);
struc_vols = cell(1,1);
for i = 1:length(struc_dir)
    if(~struc_dir(i).isdir)
        [path,name,ext] = fileparts(struc_dir(i).name);
        %Filter only original files - second letter is uppercase
        if(strcmp(ext,'.img') && strcmp(name(1),'f'))  
            struc_vols{1} = [foldername filesep 'T1scan' filesep name ext];
            break;
        end
    end
end
BATCH.Setup.structurals = struc_vols;
BATCH.Setup.RT = 2;

% Setup grey/white/CSF masks:
BATCH.Setup.masks.Grey.files{1} = GreyMask;
BATCH.Setup.masks.White.files{1} = WhiteMask;
BATCH.Setup.masks.CSF.files{1} = CSFMask;

% Choose atlas...
BATCH.Setup.rois.names = ROIName;
BATCH.Setup.rois.files{1} = ROILocation;

BATCH.Setup.nsubjects = 1; % For group analysis

%BATCH.Setup.analyses = 1; % ROI to ROI analyses, 05.10.12

BATCH.Setup.done = 1;
BATCH.Setup.conditions.names={'rest'};
BATCH.Setup.conditions.onsets{1}{1}{1} = 0;
BATCH.Setup.conditions.durations{1}{1}{1} = Inf;
BATCH.Preprocessing.done = 1;
conn_batch(BATCH);

clear BATCH;
BATCH.filename = [foldername filesep 'conn_process.mat'];
% ProcessFile = load(BATCH.filename);
% names = ProcessFile.CONN_x.Analyses.regressors.names;
% ROIName = 'atlas116';
% sources = cell(0,0);
% for i = 1:length(names)
%     if(length(names{i}) >= length(ROIName) && strcmp(names{i}(1:length(ROIName)),ROIName))
%         sources = [sources, {names{i}}];
%     end
% end
% BATCH.Analysis.sources.names = sources;
% BATCH.Analysis.sources.dimensions = repmat({1},size(sources));
% BATCH.Analysis.sources.deriv = repmat({0},size(sources));

%BATCH.Analysis.type = 1; % ROi to ROI analyses, farras, 05.10.12

BATCH.Analysis.sources = {};
BATCH.Analysis.done = 1;
BATCH.Analysis.measure = 1; % default correlation
conn_batch(BATCH);
