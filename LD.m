function [Cd,dCldmc] = LD(alpha,Re,desig,x,tmax,c,coord)

% Find the drag coefficient and lift coefficient per unit camber for a
% single hydrofoil section specified by a propeller type, at a specified
% angle of attack. The drag coefficients are, owing to a lack of other
% relevant data, calculated from data due to Burrill (1944) and Hill (1949)
% (which are experimental results pertaining to segmental and aerofoil-type
% sections). The lift-camber slope is based on data presented in O'Brien
% (1962) for Wageningen B-series sections, and includes selected effects of
% viscosity.

% TO BE ADDED:
% Dependence of drag coefficient on Reynolds number.

% INPUTS:
% - alpha : double
%           Angle of attack (rad).
% - Re    : double
%           Reynolds number.
% - desig : string
%           Propeller designation. Allowed values are: "WB" (Wageningen
%           B-series), "NACA" (NACA-series), or "C" (custom section shape).
% - x     : double
%           Radius fraction. Must be between 1.0 (obviously) and 0.2 (to
%           account for the propeller hub).
% - tmax  : double
%           Maximum thickness of section (m).
% - c     : double
%           Chord length (m).
% - coord : double, optional
%           Array of section coordinates.

% OUTPUTS:
% - Cd     : double
%            Drag coefficient.
% - dCldmc : double
%            Lift-camber slope.

% INITIAL CHECKS
% Check that the hydrofoil designation is acceptable.
if ~(desig=="WB"||desig=="NACA"||desig=="C")
    error("You have specified an invalid hydrofoil designation! Acceptable " + ...
          "designations are 'WB' (Wageningen B-series), 'NACA' (NACA series, " + ...
          "or 'C' (custom section).")
end

% Check that the radius fraction is between 0.2 and 1.0.
if x < 0.2 || x > 1.0
    error("You have specified an invalid radius fraction! Radius fractions " + ...
          "must be between 0.2 and 1.0.")
end

% DRAG COEFFICIENT
% Segmental-type sections (Burrill, 1944)
AoA = [-1 0 1 2 3 4] .* (pi/180);
tcs = [ 0.00    0.02    0.04    0.06    0.08    0.10    0.12    0.14    0.16    0.18    0.20   ];
CdS = [[0.01400 0.01179 0.01043 0.01043 0.01196 0.01434 0.01723 0.02113 0.02538 0.02996 0.03489];
       [0.01026 0.00857 0.00874 0.00992 0.01162 0.01366 0.01621 0.01926 0.02266 0.02657 0.03081];
       [0.01094 0.00942 0.00958 0.01077 0.01264 0.01485 0.01791 0.02147 0.02521 0.02928 0.03387];
       [0.01502 0.01213 0.01179 0.01281 0.01485 0.01757 0.02096 0.02470 0.02894 0.03370 0.03879];
       [0.02300 0.01842 0.01655 0.01672 0.01858 0.02113 0.02470 0.02877 0.03336 0.03811 0.04338];
       [0.03285 0.02725 0.02385 0.02266 0.02317 0.02538 0.02877 0.03302 0.03794 0.04321 0.04915]];

% Aerofoil-type sections (Hill, 1949)
CdA = [[0.01395 0.01253 0.01182 0.01210 0.01295 0.01395 0.01509 0.01651 0.01778 0.01935 0.02091];
       [0.01281 0.01168 0.01125 0.01139 0.01224 0.01352 0.01509 0.01665 0.01835 0.02006 0.02176];
       [0.01381 0.01239 0.01168 0.01196 0.01267 0.01395 0.01551 0.01722 0.01892 0.02077 0.02290];
       [0.01722 0.01452 0.01324 0.01310 0.01381 0.01509 0.01651 0.01807 0.01991 0.02162 0.02375];
       [0.02233 0.01864 0.01594 0.01466 0.01480 0.01594 0.01736 0.01920 0.02105 0.02290 0.02517];
       [0.03000 0.02489 0.02034 0.01679 0.01580 0.01707 0.01849 0.02034 0.02233 0.02432 0.02659]];

% Interpolate drag coefficient from data, depending on the section shape
if desig == "WB"
    % Wageningen B sections are aerofoil-shaped up until approximately x =
    % 0.6, after which they are round-back.
    if x < 0.6
        Cd = interp2(tcs,AoA,CdA,tmax/c,alpha,"makima");
    elseif x >= 0.6
        Cd = interp2(tcs,AoA,CdS,tmax/c,alpha,"makima");
    end
elseif desig == "NACA"
    % NACA sections are all approximately aerofoil sections
    Cd = interp2(tcs,AoA,CdA,tmax/c,alpha,"makima");
end

% LIFT-CAMBER SLOPE
% O'Brien (1962) gives values of lift-camber slopes dCl/d(m/c) for
% Wageningen B and moderate-duty standard-type sections, as functions of
% the radius fraction. These values include some effects of viscosity.
if desig == "WB"
    % Wageningen B-series sections
    rr = [0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9];
    Lm = [8.1 7.4 7.8 8.4 9.0 9.7 10.1 10.6];
    dCldmc = interp1(rr,Lm,x,"makima");
elseif desig == "NACA"
    % NACA lift-camber slopes, to be added.
end

end
