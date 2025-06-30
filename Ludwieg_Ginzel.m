function [k1,k2] = Ludwieg_Ginzel(x,phi,BAR)

% Function to calculate the Ludwieg-Ginzel curvature corrections for
% specified radius fraction, flow angle and BAR, based on empirical curve
% fits.

% INPUTS:
% - x   : double
%         Radius fraction.
% - phi : double
%         Flow angle (rad).
% - BAR : double
%         Blade area ratio.

% OUTPUTS:
% - k1 : double
%        First Ludwieg-Ginzel curvature correction.
% - k2 : double
%        Second Ludwieg-Ginzel curvature correction.

% CALCULATE LAMBDA
lambda = x * tan(phi);

% CALCULATE K1
u1 = (-0.65*(lambda^2)) + (1.1*lambda) + 0.664;
% u2 = -0.09 + ((lambda-0.2)*(0.85+((lambda-0.3)*(-4+((lambda-0.4)*(15.42-(47.95*(lambda-0.5))))))));
u2 = (0.85+(lambda-0.3)*(-4.0+(lambda-0.4)*(15.42-47.95*(lambda-0.5))));
u2 = -0.09+(lambda-0.2)*u2;
% u3 = -0.2 + ((lambda-0.2)*(1.375+((lambda-0.3)*(-3.75+((lambda-0.4)*(20.85-(75.7875*(lambda-0.5))))))));
u3 = (1.375+(lambda-0.3)*(-3.75+(lambda-0.4)*(20.85-75.7875*(lambda-0.5))));
u3 = -0.2+(lambda-0.2)*u3;
k1 = u1 + ((BAR-0.4)*(u2+(u3*(BAR-0.8))));

% CALCULATE K2
rr = [0.2    0.3    0.4    0.5    0.6    0.7    0.8    0.9   ];
C0 = [1.0    1.14   1.3    1.54   1.67   1.8    1.8    1.75  ];
C1 = [2.857  1.5    1.0    1.0    1.0    1.0    1.0    1.0   ];
C2 = [0.0    0.0    0.1    0.15   0.55   0.75   1.0    1.25  ];
C3 = [0.0    0.65   1.0    1.1    0.1667 0.3    1.333  1.5   ];
C4 = [1.0    0.1665 0.43   0.286  2.9465 2.835  1.905  3.55  ];
if x < 0.2
    % No data is available below x = 0.2, but the curvature correction
    % should be fairly constant in this region anyway.
    k2 = 1.0;
elseif x > 0.9
    % No data is available above x = 0.9, so just use linear extrapolation.
    % The fidelity of the curvature correction here isn't very important
    % because the tip does not generate much thrust.
    C0 = 1.75;
    C1 = 1.0;
    C2 = 1.25;
    C3 = 1.5;
    C4 = 3.55;
    % Extrapolate from x = 0.9, using an empirical gradient of 0.4.
    k2 = C0 + (C1*(BAR-0.4)*(C2+((BAR-0.6)*(C3+(C4*(BAR-0.9))))));
    k2 = k2 + (0.4*(x-0.9));
else
    % For radius fractions within the range of the data, interpolate for
    % the curvature correction factor.
    C0 = interp1(rr,C0,x,"makima");
    C1 = interp1(rr,C1,x,"makima");
    C2 = interp1(rr,C2,x,"makima");
    C3 = interp1(rr,C3,x,"makima");
    C4 = interp1(rr,C4,x,"makima");
    k2 = C0 + (C1*(BAR-0.4)*(C2+((BAR-0.6)*(C3+(C4*(BAR-0.9))))));
end

end
