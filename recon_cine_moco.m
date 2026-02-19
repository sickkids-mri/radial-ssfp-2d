function [ imCine, tCine, RESULTS ] = recon_cine_moco( KSpace, Times, Noise, Header, FileName, varargin )
%RECON_CINE_MOCO  Reconstruct cine 2D SSFP image series from golden-angle radial data.
%
%   [ imCine, tCine ] = RECON_CINE_MOCO( KSpace, Times, Noise, Header, FileName ) 
%   reconstructs motion-corrected cine SSFP images from radially-sampled data. 
%   Intermediate real-time reconstructions are used to estimate motion and R-wave 
%   timing of the heart for motion-correction and cardiac synchronization.
%
%   [ ..., RESULTS ] = RECON_CINE_MOCO( ... ) returns real-time image series and 
%   intermediate motion and timing estimates.
%
%   RECON_CINE_MOCO( ..., 'name', value ) specifies additional optional input arguments:
%
%       'nCineFrame'    - Number of cardiac phases (default: 15, range: 5-35)
%       'outputDir'     - Output directory for intermediate .mat files
%       'patientName'   - Patient name for DICOM metadata
%       'studyUid'      - Study UID for DICOM metadata
%       'verbose'       - Display verbose output (default: false)
%
%   Input:
%       KSpace          - K-space data [samples x spokes x channels x slices]
%       Times           - Acquisition times [spokes x slices] in milliseconds
%       Noise           - Noise covariance data for coil selection (optional)
%       Header          - Raw data header information (TWIX header)
%       FileName        - Acquisition filename
%
%   Output:
%       imCine          - Reconstructed cine images [x x frames x slices]
%       tCine           - Timing of cine frames relative to cardiac cycle in milliseconds
%       RESULTS         - Structure with reconstruction details:
%           .PREP       - Data preprocessing results
%           .MOG        - Metric-optimized gating results
%           .RT         - Real-time image series
%           .MOCO       - Motion estimation and correction results
%           .CINE       - Cine reconstruction results
%
%   See also CINE_TOOL_READ_RAW_DATA, CINE_TOOL_TRAJECTORY, CINE_TOOL_WEIGHTING

%   Based on original implementation by Chris Roy (2018).
%   Refactored to match recon_cine_pc_moco.m architecture.


%% Parse Input Arguments

p = inputParser;
addParameter( p, 'nCineFrame', 15, @(x) isnumeric(x) && x >= 5 && x <= 35 );
addParameter( p, 'outputDir', [], @(x) isempty(x) || ischar(x) || isstring(x) );
addParameter( p, 'patientName', '', @(x) ischar(x) || isstring(x) );
addParameter( p, 'studyUid', '', @(x) ischar(x) || isstring(x) );
addParameter( p, 'verbose', false, @islogical );
parse( p, varargin{:} );

rFrames = p.Results.nCineFrame;
outputDirPath = p.Results.outputDir;
patientName = p.Results.patientName;
studyUid = p.Results.studyUid;
verbose = p.Results.verbose;


%% Setup and Preprocessing

if ~isempty(outputDirPath)
    fprintf('Output path: %s \n', outputDirPath);
end
outputFileNamePrefix = strrep( FileName, 'meas_', '' );
outputFileNamePrefix = strrep( outputFileNamePrefix, '.dat', '' );

% Calculate trajectory and ramp weighting
if ~sum(strfind(Header{end}.hdr.Config.SequenceFileName(15:end),'ty'))
    k = CINE_Tool_Trajectory(size(KSpace,1), size(KSpace,2), 1);  % standard golden angle
else
    k = CINE_Tool_Trajectory(size(KSpace,1), size(KSpace,2), 7);  % tiny golden angle
end
w = CINE_Tool_Weighting(size(KSpace,1));

% Reduce number of slices for debugging purposes
fprintf('NOTE: Reconstructing one slice only, for testing purposes.\n')
indSlice=round(size(KSpace,4)/2);
KSpace=KSpace(:,:,:,indSlice);
Times=Times(:,indSlice);

% Remove noisy channels
RESULTS.PREP.initialKSpaceSize = size(KSpace);
fprintf('Initial KSpace Size: %s\n', mat2str(size(KSpace)))

if ~isempty(Noise)
    Noise = squeeze(mean(mean(mean(abs(Noise),1),3),6));
    Coils = 1:length(Noise);
    nThresh = 2;
    while sum(Noise > nThresh*mean(Noise(~isnan(Noise))))
        Coils(Noise > nThresh*mean(Noise(~isnan(Noise)))) = nan;
        Noise(Noise > nThresh*mean(Noise(~isnan(Noise)))) = nan;
    end
    KSpace = KSpace(:,:,Coils(~isnan(Coils)),:);
end
clear Noise Coils
fprintf('KSpace Size (After Removing Noisy Channels): %s\n', mat2str(size(KSpace)))

% Perform coil compression
try
    [~,~,V] = svd(reshape(permute(KSpace,[1,2,4,3]),size(KSpace,1)*size(KSpace,2)*size(KSpace,4),size(KSpace,3)),'econ');
    KSpace = permute(reshape(reshape(permute(KSpace,[1,2,4,3]),size(KSpace,1)*size(KSpace,2)*size(KSpace,4),size(KSpace,3))*V,size(KSpace,1),size(KSpace,2),size(KSpace,4),size(KSpace,3)),[1,2,4,3]);
    KSpace = KSpace(:,:,1:8,:);
    clear V
catch
    for iSlice = 1:size(KSpace,4)
        [~,~,V] = svd(reshape(KSpace(:,:,:,iSlice),size(KSpace,1)*size(KSpace,2),size(KSpace,3)),'econ');
        KSpace(:,:,:,iSlice) = reshape(reshape(KSpace(:,:,:,iSlice),size(KSpace,1)*size(KSpace,2),size(KSpace,3))*V,size(KSpace,1),size(KSpace,2),size(KSpace,3));
    end
    KSpace = KSpace(:,:,1:8,:);
    clear V
end
fprintf('KSpace Size (After Coil Compression): %s\n', mat2str(size(KSpace)))

% Estimate coil sensitivities
FT = MCNUFFT(k,repmat(w,[1,size(k,2)]).^2,ones(0.5*size(KSpace,1),0.5*size(KSpace,1)));
STACKS = zeros(0.5*size(KSpace,1),0.5*size(KSpace,1),size(KSpace,3),size(KSpace,4),'single');
for iSlice = 1:size(KSpace,4)
    for iCoil = 1:size(KSpace,3)
        STACKS(:,:,iCoil,iSlice) = FT'*KSpace(:,:,iCoil,iSlice);
    end
end
clear FT

csm = zeros(size(STACKS),'single');
try
    parfor iSlice = 1:size(STACKS,4)
        csm(:,:,:,iSlice) = ismrm_estimate_csm_walsh(STACKS(:,:,:,iSlice));
    end
catch
    for iSlice = 1:size(STACKS,4)
        csm(:,:,:,iSlice) = ismrm_estimate_csm_walsh(STACKS(:,:,:,iSlice));
    end
end
clear STACKS


%% Real-Time Image Reconstruction (low temporal resolution, for motion correction)

fprintf('Step 2 of 6 - Reconstructing Real-Time Images \n');

[~,kRT] = CINE_Tool_RT_Sort(ones(size(k)), k, 15, 10);  % pre-compute re-sorted trajectory

param.E = MCNUFFT(kRT,repmat(w,[1,size(kRT,2),size(kRT,3)]).^2,csm(:,:,:,1));
RT = zeros(size(csm,1),size(csm,1),size(kRT,3),size(csm,4),'single');
param.lambda = 0.5;  % temporal total variation regularizer
param.nite = 4;      % number of iterations

fprintf('Step 1 of 6 - Loading Raw Data - Finished \n');

for iSlice = 1:size(csm,4)
    fprintf('    slice %i of %i\n', iSlice, size(csm,4))
    param.E = swap_b1(param.E,csm(:,:,:,iSlice));
    param.y = CINE_Tool_RT_Sort(KSpace(:,:,:,iSlice), kRT, 15, 10);
    RT(:,:,:,iSlice) = CINE_Tool_RT_Recon(param);
end
clear param kRT

if ~isempty(outputDirPath)
    matFileName = strcat('reconstructed_RT_images_', outputFileNamePrefix, '.mat');
    matFilePath = fullfile(outputDirPath, matFileName);
    fprintf('   %s\n', matFilePath);
    save(matFilePath, 'RT', '-v7.3');
end

RESULTS.PREP.kSize = size(KSpace);
RESULTS.PREP.nCoils = size(KSpace,3);
RESULTS.PREP.nSlices = size(KSpace,4);


%% Motion Estimation

fprintf('Step 3 of 6 - Calculating Translational Motion \n');

t = imresize(Times(:,1)-min(Times(:,1)), size(RT,3)/length(Times(:,1)));

[ROIy, ROIx] = CINE_Tool_AutoROI(RT, t);
[MIy, MIx] = CINE_Tool_Enlarge_ROI(ROIy, ROIx, round(size(RT,1)/2), round(size(RT,1)/2));

if max(MIy) > size(RT,1); MIy = MIy+(max(MIy)-size(RT,1)); end
if max(MIx) > size(RT,1); MIx = MIx-(max(MIx)-size(RT,1)); end

AVGFrames = round(range(t)/200);
RTAVG = CINE_Tool_Temporal_Averaging(RT, AVGFrames);
MI = CINE_Tool_Mutual_Information(Decimate_CINE(RTAVG(MIy,MIx,:,:), 4));

[MOCOy, MOCOx] = CINE_Tool_Enlarge_ROI(ROIy, ROIx, round(size(RT,1)/4), round(size(RT,1)/4));

if max(MOCOy) > size(RT,1); MOCOy = MOCOy+(max(MOCOy)-size(RT,1)); end
if max(MOCOx) > size(RT,1); MOCOx = MOCOx-(max(MOCOx)-size(RT,1)); end

myTransforms = zeros(AVGFrames, 2, size(RT,4));

for iSlice = 1:size(RT,4)
    Targets = repmat(RTAVG(:,:,MI(:,iSlice)==max(MI(:,iSlice)),iSlice), [1,1,AVGFrames]);
    myTransforms(:,:,iSlice) = CINE_Tool_MOCO(abs(RTAVG(:,:,:,iSlice)), abs(Targets), MOCOy, MOCOx);
    fprintf('    slice %i\n', iSlice)
end

RT = CINE_Tool_Correct_RT(RT, -myTransforms);

magTempMeanRT = abs(mean(RT, 3));
if ~isempty(outputDirPath)
    matFileName = strcat('calculated_translational_motion_', outputFileNamePrefix, '.mat');
    matFilePath = fullfile(outputDirPath, matFileName);
    fprintf('   %s\n', matFilePath);
    save(matFilePath, 'magTempMeanRT', 'MI', 'ROIx', 'ROIy', 'MIx', 'MIy', 'myTransforms', '-v7.3');
    matFileName = strcat('motion-corrected_RT_images_', outputFileNamePrefix, '.mat');
    matFilePath = fullfile(outputDirPath, matFileName);
    fprintf('   %s\n', matFilePath);
    save(matFilePath, 'RT', '-v7.3');
end

clear RTAVG MOCOy MOCOx AVGFrames
fprintf('Step 3 of 6 - Calculating Translational Motion - Finished \n');


%% Cardiac Synchronization (Metric-Optimized Gating)

fprintf('Step 4 of 6 - Calculating Heart Rate \n');

[ROIy, ROIx, Initial_RRs] = CINE_Tool_AutoROI(RT, t);
[MIy, MIx] = CINE_Tool_Enlarge_ROI(ROIy, ROIx, round(size(RT,1)/2), round(size(RT,1)/2));

if max(MIy) > size(RT,1); MIy = MIy+(max(MIy)-size(RT,1)); end
if max(MIx) > size(RT,1); MIx = MIx-(max(MIx)-size(RT,1)); end

MI = CINE_Tool_Mutual_Information(Decimate_CINE(RT(MIy,MIx,:,:), 4));
clear MIy MIx

MOG = struct('RWaveTimes',[],'Entropy',[]);
for iSlice = 1:size(RT,4)
    MOG(iSlice,1) = CINE_Tool_MOG(imgaussfilt(abs(RT(ROIy,ROIx,:,iSlice)), 2), t, [], 20, MI(:,iSlice), Initial_RRs(iSlice));
end
clear RT Initial_RRs

fprintf('Step 4 of 6 - Calculating Heart Rate - Finished \n');


%% Cine Image Reconstruction (high quality, motion-corrected, cardiac-gated)

fprintf('Step 5 of 6 - Reconstructing CINE Images \n');

CINEs = zeros(0.5*size(KSpace,1), 0.5*size(KSpace,1), rFrames, size(KSpace,4), 'single');

for iSlice = 1:size(KSpace,4)
    fprintf('    slice %i of %i\n', iSlice, size(KSpace,4))
    
    CP = Calculate_CardiacPhases(Times(:,iSlice)-min(Times(:,iSlice)), MOG(iSlice).RWaveTimes);
    MOCO = exp(-2*pi*1i*(-bsxfun(@times,imresize(-myTransforms(:,2,iSlice),[size(KSpace,2),1])',real(k))-...
                           bsxfun(@times,imresize(-myTransforms(:,1,iSlice),[size(KSpace,2),1])',imag(k))));
    
    W = KSpace_Interpolation_Weights(size(KSpace,1), CP, MI(:,iSlice), rFrames);
    
    [~,SpokeIndices] = sort(squeeze(W(1,:,:)), 1, 'descend');
    rSpokes = max(sum(squeeze(W(1,:,:)) > max(W(:))/1000));
    
    fprintf('    spokes per frame: rSpokes(slice=%2i) = %5i (%6i acquired spokes, %2i frames)\n', iSlice, rSpokes, size(KSpace,2), rFrames)
    
    kMOG = zeros(size(KSpace,1), rSpokes, rFrames, 'single');
    param.y = zeros(size(KSpace,1), rSpokes, size(KSpace,3), rFrames, 'single');
    
    for iFrame = 1:rFrames
        param.y(:,:,:,iFrame) = bsxfun(@times,bsxfun(@times,KSpace(:,SpokeIndices(1:rSpokes,iFrame),:,iSlice),...
                                       MOCO(:,SpokeIndices(1:rSpokes,iFrame))),W(:,SpokeIndices(1:rSpokes,iFrame),iFrame));
        kMOG(:,:,iFrame) = k(:,SpokeIndices(1:rSpokes,iFrame));
    end
    
    param.E = MCNUFFT(kMOG,repmat(w,[1,size(kMOG,2),size(kMOG,3)]).^2,csm(:,:,:,iSlice));
    param.lambda = 0.1;
    param.nite = 9;
    param.tol = 0.0025;
    
    CINEs(:,:,:,iSlice) = CINE_Tool_MOG_Recon(param);
    
    srcDelay = zeros(1,size(CINEs,4));
    if (iSlice ~= 1)
        fprintf('    Cardiac Phase Alignment\n');
        maskTgt = zeros(size(CINEs(:,:,1)));
        maskTgt(ROIy,ROIx,:) = 1;
        [CINEs(:,:,:,iSlice), srcDelay(iSlice)] = align_cardiac_phase_cine2cine(CINEs(:,:,:,iSlice-1), maskTgt, CINEs(:,:,:,iSlice), maskTgt);
    end
end

if ~isempty(outputDirPath)
    matFileName = strcat('reconstructed_CINE_images_', outputFileNamePrefix, '.mat');
    matFilePath = fullfile(outputDirPath, matFileName);
    fprintf('   %s\n', matFilePath);
    save(matFilePath, 'CINEs', 'MI', 'MOG', 'srcDelay', 'ROIx', 'ROIy', '-v7.3');
end

fprintf('Step 5 of 6 - Reconstructing CINE Images - Finished \n');


%% Save as DICOM Files

if ~isempty(outputDirPath) && ~isempty(patientName) && ~isempty(studyUid)
    fprintf('Step 6 of 6 - Saving DICOMs \n');
    save_cine_to_dicom(CINEs, MOG, Header, FileName, rFrames, outputDirPath, patientName, 'verbose', verbose);
else
    fprintf('Step 6 of 6 - Saving DICOMs - Skipped (missing parameters) \n');
end


%% Prepare Output

imCine = CINEs;

% Calculate timing for cine frames
numFrame = size(CINEs,3);
tCine = zeros(numFrame,1);
rrSlice = median(diff(double(MOG(1).RWaveTimes)));
dtFrame = rrSlice / numFrame;
tCine = dtFrame * (0:(numFrame-1))';

% Prepare results structure
RESULTS.DATA.nSample = size(KSpace,1);
RESULTS.DATA.nSpoke = size(KSpace,2);
RESULTS.DATA.nSlice = size(KSpace,4);
RESULTS.MOG.MOG = MOG;
RESULTS.MOG.MI = MI;
RESULTS.MOCO.myTransforms = myTransforms;
RESULTS.MOCO.ROIx = ROIx;
RESULTS.MOCO.ROIy = ROIy;
RESULTS.CINE.imCine = CINEs;
RESULTS.CINE.tCine = tCine;

end  % recon_cine_moco(...)

