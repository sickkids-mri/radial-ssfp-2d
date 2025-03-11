function [ imCine, tCine, RESULTS ] = recon_cine_moco( kData, kTraj, kTimes, kSize, DATA, RECON, MOCO, MOG, CINE )
%RECON_CINE_MOCO  Reconstruct motion-corrected cine 2D MR image series from golden-angle radial data.
%
%   [ imCine, tCine ] = RECON_CINE_MOCO( kData, kTraj, kTimes, kSize ) reconstructs
%   motion-corrected MR images from radially-sample data. Intermediate real-time
%   reconstructions are used to estimate motion and R-wave timing of the fetal heart for
%   motion-correction, data rejection and cardiac synchronization.
%
%   [ ..., RESULTS ] = recon_cine_moco( ... ) returns real-time image series and intermediate
%   motion and timing estimates in structure RESULTS.
%
%   recon_cine_moco( ..., 'name', value ) specifies additional, optional input arguments.
%
%   Input
%       kData           k-space data
%       kTraj           k-space trajectories
%       kTimes          array of acquisition time for each spoke, in milliseconds
%       kSize           dimensions of k-space data, specified as vector with each entry
%                       corresponding to dimensions as
%                           ['Col','Cha','Lin','Par','Sli','Ave','Phs','Eco',...
%                            'Rep','Set','Seg','Ida','Idb','Idc','Idd','Ide']
%
%   Output
%       imCine          cine image series after compressed sensing
%       tCine           timing of cine image frames
%
%   Refer to the arguments block in the code for an extensive list of required and optional input
%   arguments along with details of size, class, attributes, default values and validation.
%
%   Example
%
%   See also


%% Argument Validation

% Arguments Block

arguments

  % name                size        class       validation_functions                default_value   units / notes
  % -----------------   ---------   ---------   ---------------------------------   -------------   -------------

    % Positional Arguments

    kData                           double      {mustBeNonempty}
    kTraj                           double      {mustBeNonempty}
    kTimes                          double      {mustBeNonempty,mustBeReal,mustBeNonnegative}
    kSize               (1,:)       double      {mustBeNonempty,mustBePositive}

    % Name-Value Arguments

    DATA.TR             (1,1)       double      {mustBePositive}                                    % milliseconds

    RECON.sliceNo       (1,1)       double      {mustBePositive,mustBeInteger}      = 1

    MOCO.mocoRltDt      (1,1)       double      {mustBePositive}                    = 736           % milliseconds

    MOG.mogRltDt        (1,1)       double      {mustBePositive}                    = 46            % milliseconds

    CINE.nCineFrame     (1,1)       double      {mustBePositive,mustBeInteger}      = 15            % number of cine frames to reconstruct
    CINE.nCineIterCs    (1,1)       double      {mustBePositive,mustBeInteger}      = 10            % number of compressed sensing iterations

end

% Further Argument Validation

assert( ~isreal( kData ), 'kData must be complex-valued' )
assert( ~isreal( kTraj ), 'kTraj must be complex-valued' )

% Calculate Input-Depedent Optional Arguments

if ~isfield(DATA,'TR')
    DATA.TR = median(diff(kTimes));
end


%% Setup

DATA.nSample     = kSize(1);
DATA.nChannel    = kSize(2);
DATA.nSpoke      = kSize(3);
DATA.nPar        = kSize(4);
DATA.nSlice      = kSize(5);
DATA.nAve        = kSize(6);
DATA.nPhs        = kSize(7);
DATA.nEco        = kSize(8);
DATA.nRep        = kSize(9);
DATA.nSet        = kSize(10);
DATA.nSeg        = kSize(11);
DATA.nIda        = kSize(12);
DATA.nIdb        = kSize(13);
DATA.nIdc        = kSize(14);
DATA.nIdd        = kSize(15);
DATA.nIde        = kSize(16);


%% TODO: Add Reconstruction


%% Assign Output Arguments

imCine          = [];
tCine           = [];

% Remove K-space Data to Reduce Size
RESULTS.DATA    = DATA;
RESULTS.RECON   = RECON;
RESULTS.MOCO    = MOCO;     % rmfield( MOCO, { 'kData', 'kTraj' } );
RESULTS.MOG     = MOG;      % rmfield( MOG,  { 'kData', 'kTraj' } );
RESULTS.CINE    = CINE;     % rmfield( CINE, { 'kData', 'kTraj' } );


end  % recon_cine_moco(...)
