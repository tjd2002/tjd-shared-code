function [R B] = target3dm(target, tilt, rot, tiltdir, dvintersect, quiet)
% TARGET3DM: target a location using a tilted, rotated stereotax
%
% INPUTS:
% target = [ML AP DV] atlas coordinates of target location, relative to
%        Bregma:
%          ML: animal right is +ve
%          AP: anterior is +ve
%          DV: dorsal (up) is +ve
%          (A.K.A.: R/A/D coords, same as used by Kopf digital display)
%
% tilt = degrees of CCW tilt of DV axis away from vertical, *before any
%        rotation*: +ve angles = tilt left/back, -ve angles = tilt right/fwd. 
%
% rot = degrees of CCW rotation of tilted DV axis away from animal right
%        (i.e. 'east', +ve X axis), after tilting. Must be within 0+/-30
%        (tiltLR), or 90 +/- 30 (~tiltLR)
%
% tiltdir = 'LR', or 'FB': is the stereotax arm configured so that the arm
%        tilts left-right, or front-back?
%
% dvintersect = DV heights at which to calculate intersection of injection 
%        vector (in old coordinates). Default = []. Useful for 
%        establishing craniotomy locations in unrotated coordinates before 
%        rotating the stereotax arms.
%
% quiet = true/false: silence text description of outputs, just return
%        coordinates.
%
% OUTPUTS:
% R = coordinates to use to reach the targeted location, using the 
%        rotated/tilted stereotax arm.
%
% B = location, in *unrotated* coordinates, of the intersection of the 
%        rotated injection vector with a DV plane specified by dvintersect. 
%
% EXAMPLES:
%  To target a location at ML 0.25, AP 0.7, DV -6.6, with a stereotax arm 
%  tilted 20 degrees to the right, and with no rotation:
%     [R] = target3dm([0.25 0.7 -6.6], -20, 0, 'LR');
%
%
% tjd@stanford.edu 2010-2013

% TODO:
% -make inputs/outputs optionally same order/sign as BrainNav uses
% -enforce that rot can only be [0,90]+/-30, even for left arm. 180/270 
%  would mean that sign of ML axis was reversed. test whether tiltLR is 
%  correct, given rot.

switch lower(tiltdir)
    case 'lr'
        tiltLR = true;
    case 'fb'
        tiltLR = false;
    otherwise
        error('Must provide ''tiltdir'': either ''LR'', or ''FB''');
end

if ~exist('dvintersect', 'var')
    dvintersect = [];
end


if ~exist('quiet', 'var')
    quiet = false;
end

% Convert *DV vector* to polar coordinates as used by sph2cart:
% phi = elevation of DV axis from x0y0 plane (i.e. from horizontal)
% theta = CCW rotation of DV axis in x0y0 plane away from positive 'X' 
%       (i.e. animal right=0, nose=90)
phi = deg2rad(tilt+90);
theta = deg2rad(rot);

    
% x1 (ML) vector in x0y0z0 space

if tiltLR, % Tilt Left-Right  case:
    % Relative to DV, rotation (theta) is same, tilt (phi) is 90deg CW)
    [x01(1) x01(2) x01(3)] = sph2cart(theta, phi-pi/2,1);

else % Tilt Front-Back case
    
    % x1 (ML) vector
    % ML not tilted (phi=0); rotation (theta) 90deg CW relative to DV axis
    [x01(1) x01(2) x01(3)] = sph2cart(theta-pi/2, 0, 1);

end

% y1 (AP) vector in x0y0z0 space(AP stereotax axis can't tilt or rotate)
y01 = [0 1 0];

% z1 (DV) vector in x0y0z0 space (tilted, rotated as per inputs)
[z01(1) z01(2) z01(3)] = sph2cart(theta, phi,1);

% the inverse of this matrix can be used to get coords in x1y1z1 given
% target in x0y0z0
rotmat = inv([x01' y01' z01']);
R = rotmat * target';

% output as row vector
R = R';


% get intersection of injection vector with x0y0 plane (Bregma plane), in
% x0y0z0 coordinates: (no dependence on rotated coords x1y1z1)

if ~isempty(dvintersect)
  % radial distance from target to Bregma plane
  Br = [-target(3)+dvintersect] / sin(phi);

  % x0y0 distance from target to Bregma plane
  [B(:,1) B(:,2) B(:,3)] = sph2cart(theta,phi,Br);

  % plus x0y0 distance to target
  B = bsxfun(@plus, B, target);
else
  B = [];
end
  
if ~quiet
    disp(sprintf('\nAll coords in R/A/D (as on Kopf stereotax digital display)'));
    disp(sprintf('R: %+0.2f / %+0.2f / %+0.2f (target location in new coords)\n', R));
    if ~isempty(B),
      disp(sprintf('B: %+0.2f / %+0.2f / %+0.2f (intersection with DV plane in old coords)\n', B'));
    end
end
