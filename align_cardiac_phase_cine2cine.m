function [ imCineSrcAlign, srcDelay ] = align_cardiac_phase_cine2cine( imCineTgt, maskTgt, imCineSrc, maskSrc, showFig );

%% Inputs

if ~exist('showFig','var')
  showFig = false;
end

%% Anon Fns

get_mask_values_vec = @( x, msk ) double( x( repmat( msk, [1,1,size(x,3)/size(msk,3)] ) ) );

get_mask_values_3d  = @( im, msk ) reshape( get_mask_values_vec( im, msk ), sum(msk(:)), 1, size(im,3) );

make_vec            = @( x ) reshape( x , [], 1 );
%% Create Common Mask
mask = maskTgt & maskSrc;

%% Get SI in Mask

siTgt = get_mask_values_3d( imCineTgt, mask );
siSrc = get_mask_values_3d( imCineSrc, mask );
nCardPhase = size(siTgt,3);
%% Find Delay

delay = 1:nCardPhase;

% cc1d = nan(1,nCardPhase);

cc3d = nan(1,nCardPhase);
for iD = 1:numel(delay),
    si2delay = circshift( siSrc, delay(iD), 3 );
%     cc1d(iD) = corr( make_vec( mean( siTgt, 1 ) ), make_vec( mean( si2delay, 1 ) ) );
    cc3d(iD) = corr( make_vec( siTgt ), make_vec( si2delay ) );
end

% [ ~, iDelay1d ] = max( cc1d );

[ ~, iDelay3d ] = max( cc3d );

% delay1d = delay(iDelay1d);

delay3d = delay(iDelay3d);

%% Offset

iDelay  = iDelay3d;
cc      = cc3d;
srcDelay = delay(iDelay);
imCineSrcAlign = circshift( imCineSrc, srcDelay, 3);

%% Visualise

if ( showFig )
  figure( 'Name', 'align_cardiac_phase' )
  yyaxis left
  plot( squeeze( abs( mean( siTgt, 1 ) ) ), '-' )
  hold on
  plot( squeeze( abs( mean( siSrc, 1 ) ) ), ':' )
  plot( squeeze( abs( mean( circshift( siSrc, delay(iDelay), 3 ), 1 ) ) ), '--' )
  ylabel( '|Y_{j \in ROI}|')
  xlabel( 'cardiac phase / offset' )
  legend('Y1','Y2','Y2 offset')
  yyaxis right
  plot( abs( cc ) )
  hold on
  plot( iDelay, abs(cc(iDelay)), 'o' )
  ylabel( 'correlation' )
end

end % align_cardiac_phase_cine2cine(...)