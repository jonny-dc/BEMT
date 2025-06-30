function K = Goldstein(x,phi,Z)

% Find the Goldstein correction factor based on the number of blades,
% radius fraction and inflow angle.

% INPUTS:
% - x   : double
%         Radius fraction.
% - phi : double
%         Effective inflow angle (rad).
% - Z   : integer
%         Number of blades.

% OUTPUTS:
% - K : double
%       Goldstein correction factor.

if x == 1.0
    K = 0;
else
    F = (Z/(2*x*tan(phi))) - (1/2);
    if F <= 85
        K = (2/pi) * acos(cosh(x*F)/cosh(F));
    else
        K = 1;
    end
end

end
