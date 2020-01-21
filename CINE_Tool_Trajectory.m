%=========================================================================
% Create radial trajectory for NUFFT
function k = CINE_Tool_Trajectory(no_samples, no_profiles, gatype)
%==========================================================================
% Returns radial trajectory in a complex valued data array. The real part
% corresponds to the x- and the imaginary part to y-component respectively.
%
% Inputs:
% -------
% no_samples:   Number of samples along each projection.
% no_rofiles:   Number of projections per frame and slice.
% gatype:         Golden angle type (i.e. 1 - standard GA, 7 - tiny golden angle
%
% Outputs:
% --------
% k:            k-space coordinates.
%
% Function calls:    none
% ---------------
%
% Claudio Santelli, 11/06/2012
% Chris Roy, 06/28/2018
%==========================================================================

% Build trajectory
%--------------------------------------------------------------------------

% Initial spoke along ky-axis
k0 = [zeros(1,no_samples); -floor(no_samples/2):ceil(no_samples/2-1)];


% Angle increment


dPhi=pi/(gatype-1+(1+sqrt(5))/2);

k=zeros(no_samples,no_profiles);
for i=1:no_profiles
    % Update rotation matrix
    rot_angle = (i-1)*dPhi;
    
    R = [cos(rot_angle), -sin(rot_angle);
        sin(rot_angle),  cos(rot_angle)];
    % Rotate k0 vector accordingly
    ktmp       = (R*k0).';
    k(:,i) = ktmp(:,1)+1i*ktmp(:,2);
end

k = k./no_samples;
k=single(k);
end



