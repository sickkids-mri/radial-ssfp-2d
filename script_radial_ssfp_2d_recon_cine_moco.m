% Example golden-angle radial 2D SSFP cine MRI reconstruction with motion-correction and cardiac metric optimized gating
%
%   Requires variable twixFilePath at minimum.
%
%   Example
%       twixFilePath    = '/absolute/path/to/raw/meas_id.dat';
%       outputDirPath   = '/absolute/path/to/output/meas_id/';
%       script_radial_ssfp_2d_recon_cine_moco

% NOTE: script is an example, however full functionality has not been tested

%% Initialize

% clear all variables except for select user-specified files, directories
% and parameters
clearvars -except debug saveLog saveResults saveFigures ...
    twixFilePath outputDirPath figDirPath dependenciesDirPath...
    nCineFrame

% close all open figures
close all

% add dependencies to Matlab path
%   e.g., 
%       dependenciesDirPath = '/path/to/dependencies';
%       origPath  = path;
%       resetPath = onCleanup( @() path( origPath ) );
%       addpath( genpath( dependenciesDirPath ) );

%% Setup

fprintf( '%s\n', mfilename )

% flags
if ~exist('debug','var'),       debug       = false;    end
if ~exist('saveLog','var'),     saveLog     = true;     end
if ~exist('saveResults','var'), saveResults = true;     end
if ~exist('saveFigures','var'), saveFigures = true;     end

% twix .dat file
assert( exist( 'twixFilePath', 'var' ) == 1 )
assert( exist( twixFilePath, 'file' ) == 2 )
fprintf( 'reconstructing %s\n', twixFilePath )
[ twixDirPath, twixFileName, twixFileExt ] = fileparts( twixFilePath );
outputFilePrefix = strrep( twixFileName, 'meas_', '' );

% output file directory
if ~exist('outputDirPath','var')
    outputDirPath = twixDirPath;
else
    fprintf( 'Using user-specified outputDirPath = %s\n', outputDirPath )
end
if ~exist(outputDirPath,'dir')
    mkdir( outputDirPath )
end

% results .mat file
resultsFilePath = fullfile( outputDirPath, strcat( outputFilePrefix, '_results.mat' ) );

% log .txt file
logFilePath = fullfile( outputDirPath, strcat( outputFilePrefix, '_log.txt' ) );

% figure directory
if ~exist('figDirPath','var')
    figDirPath = fullfile( outputDirPath, strcat( outputFilePrefix, '_figures' ) );
else
    fprintf( 'Using user-specified figDirPath = %s\n', figDirPath )
end

% number of cine frames
if ~exist('nCineFrame','var')
    nCineFrame = 15;
else
    fprintf( 'Using user-specified nCineFrame = %i\n', nCineFrame )
end


%% Start Logging

if ( saveLog )
    diary( logFilePath )
end


%% Load Raw Data

fprintf( 'Loading raw data...\n' );
[KSpace, Times, Noise, Header] = load_raw_data( twixFilePath );

% Extract basic parameters
nSample = size(KSpace, 1);
nSpoke = size(KSpace, 2);
nChannel = size(KSpace, 3);
nSlice = size(KSpace, 4);

fprintf( 'Data dimensions: %d samples x %d spokes x %d channels x %d slices\n', ...
    nSample, nSpoke, nChannel, nSlice );


%% Run Reconstruction

fprintf( 'Running reconstruction...\n' );

% Call the main reconstruction function
[ imCine, tCine, RESULTS ] = recon_cine_moco( KSpace, Times, Noise, Header, twixFileName, ...
    'nCineFrame', nCineFrame, ...
    'outputDir', outputDirPath, ...
    'verbose', true );

fprintf( 'Reconstruction completed successfully.\n' );


%% Save Results

if ( saveResults )
    fprintf( 'Saving results to %s...\n', resultsFilePath );
    save( resultsFilePath, 'imCine', 'tCine', 'RESULTS', '-v7.3' );
end


%% Generate Summary Figures

if ( saveFigures )
    fprintf( 'Generating summary figures...\n' );

    % Create figure directory
    if ~exist(figDirPath,'dir')
        mkdir( figDirPath )
    end

    % Figure 1: Final cine images
    figure('Name', 'cine_images', 'Position', [100, 100, 1200, 800]);
    nFrames = size(imCine, 3);
    nRows = ceil(sqrt(nFrames));
    nCols = ceil(nFrames / nRows);

    for iFrame = 1:nFrames
        subplot(nRows, nCols, iFrame);
        imagesc(abs(imCine(:,:,iFrame,1)));
        axis image off;
        colormap gray;
        title(sprintf('Frame %d (%.1f ms)', iFrame, tCine(iFrame)));
    end
    suptitle('Reconstructed Cine Images');

    % Figure 2: Motion correction results
    if isfield(RESULTS, 'MOCO') && isfield(RESULTS.MOCO, 'myTransforms')
        figure('Name', 'motion_correction', 'Position', [100, 100, 1000, 600]);

        subplot(2,2,1);
        imagesc(abs(RESULTS.RT.imRt(:,:,1,1)));
        axis image off;
        colormap gray;
        title('Real-time image (before MoCo)');

        subplot(2,2,2);
        imagesc(abs(RESULTS.CINE.imCine(:,:,1,1)));
        axis image off;
        colormap gray;
        title('Cine image (after MoCo)');

        subplot(2,2,3:4);
        plot(RESULTS.MOCO.myTransforms(:,1), 'b-', 'LineWidth', 2);
        hold on;
        plot(RESULTS.MOCO.myTransforms(:,2), 'r-', 'LineWidth', 2);
        grid on;
        xlabel('Frame');
        ylabel('Translation (pixels)');
        legend('X translation', 'Y translation');
        title('Motion Correction Transforms');
    end

    % Figure 3: Cardiac synchronization
    if isfield(RESULTS, 'MOG') && isfield(RESULTS.MOG, 'MOG')
        figure('Name', 'cardiac_sync', 'Position', [100, 100, 800, 600]);

        subplot(2,1,1);
        plot(RESULTS.MOG.MOG.RWaveTimes, 'ro-', 'LineWidth', 2, 'MarkerSize', 8);
        grid on;
        xlabel('R-wave number');
        ylabel('Time (ms)');
        title('Detected R-wave Times');

        subplot(2,1,2);
        plot(diff(RESULTS.MOG.MOG.RWaveTimes), 'b-', 'LineWidth', 2);
        grid on;
        xlabel('Beat number');
        ylabel('RR interval (ms)');
        title('RR Intervals');
    end

    % Save all figures
    save_figs( fullfile(figDirPath,'png'), [], fullfile(figDirPath,'fig') );
end


%% Display Summary

fprintf( '\n=== Reconstruction Summary ===\n' );
fprintf( 'Input file: %s\n', twixFilePath );
fprintf( 'Output directory: %s\n', outputDirPath );
fprintf( 'Cine dimensions: %d x %d x %d x %d (x,y,frames,slices)\n', size(imCine) );
fprintf( 'Number of cardiac phases: %d\n', nCineFrame );
fprintf( 'Temporal resolution: %.1f ms per frame\n', mean(diff(tCine)) );

if isfield(RESULTS, 'MOG') && isfield(RESULTS.MOG, 'MOG')
    rrIntervals = diff(RESULTS.MOG.MOG.RWaveTimes);
    fprintf( 'Mean heart rate: %.1f bpm\n', 60000 / mean(rrIntervals) );
    fprintf( 'Heart rate variability: %.1f ms\n', std(rrIntervals) );
end

fprintf( 'Total reconstruction time: %.1f seconds\n', toc );

if saveResults
    fprintf( 'Results saved to: %s\n', resultsFilePath );
end

if saveFigures
    fprintf( 'Figures saved to: %s\n', figDirPath );
end

fprintf( '=== End Summary ===\n\n' );


%% Stop Logging

if ( saveLog )
    diary off
end


%% Clean-Up

resetPath.task();


%% Helper Functions

function save_figs( pngDir, figHandles, figDir )
%SAVE_FIGS Save figures to PNG and FIG formats
%
%   SAVE_FIGS( pngDir, figHandles, figDir )
%
%   Inputs:
%       pngDir      - Directory to save PNG files
%       figHandles  - Figure handles to save (empty for all)
%       figDir      - Directory to save FIG files (optional)

    if nargin < 2 || isempty(figHandles)
        figHandles = get(groot,'Children');
    end

    if nargin < 3
        figDir = [];
    end

    % Create directories if they don't exist
    if ~exist(pngDir,'dir')
        mkdir(pngDir);
    end
    if ~isempty(figDir) && ~exist(figDir,'dir')
        mkdir(figDir);
    end

    % Save each figure
    for i = 1:length(figHandles)
        hFig = figHandles(i);
        figName = hFig.Name;
        if isempty(figName)
            figName = sprintf('figure_%02d', i);
        end

        % Save PNG
        pngPath = fullfile(pngDir, [figName '.png']);
        print(hFig, pngPath, '-dpng', '-r300');

        % Save FIG if requested
        if ~isempty(figDir)
            figPath = fullfile(figDir, [figName '.fig']);
            savefig(hFig, figPath);
        end

        fprintf('Saved figure: %s\n', figName);
    end
end
