function w = CINE_Tool_Weighting(no_samples)
%==========================================================================
% Function which calculates analytically the sampling weights of a radial
% trajectory.
%
% Inputs:
% -------
% no_samples:       Number of samples along each projection.
% no_profiles:      Number of projections per frame and slice.
% alt_prof:         Alternating profiles flag.
% gafl:             Golden angle flag.
% normfl:           Normalization flag.
%
% Outputs:
% --------
% w:                Radial Weights
%
% Function calls:    none
% ---------------
%
% Claudio Santelli, 06/11/2012
%==========================================================================

% Check input
%--------------------------------------------------------------------------

% Initialize weight and build up kr vector (radial coordinate along spoke)
%--------------------------------------------------------------------------
w = abs(-floor(no_samples/2):ceil(no_samples/2-1));
w(w==0) = 1/4;
 w = w./max(w(:));
w=single(w');
end


