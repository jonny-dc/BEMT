function [KT,KQ,eta,params] = BEMT_Solver(J,Re,D,Z,BAR,desig,x,c,tmax,m,beta)

% This function uses iteration to solve the blade element-momentum
% propeller problem, using the procedure outlined in Molland, Turnock &
% Hudson, (2011) Ship Resistance and Propulsion. Details on the
% implementation and theory can be found there. This program depends on the
% following subroutines:
% - Wageningen.m : Find and discretise the geometry of a Wageningen
%       B-series propeller.
% - Wageningen_KTKQ.m : Find the thrust and torque characteristics of a
%       Wageningen B-series propeller, to compare to the blade-element
%       momentum prediction.
% - LD.m : Find the lift and drag coefficients corresponding to a given
%       section and angle of attack.
% - Goldstein.m : Find the Goldstein correction factor for given section
%       and number of blades.
% - Ludwieg_Ginzel.m : Find the Ludwieg-Ginzel curvature correction factors
%       (as presented by Eckhardt and Morgan (1955)).
% - Req_AoA.m : Find the required angle of attack corresponding to a
%       pre-defined lift coefficient and section.

% Please note that, due to the specifics of this method, you cannot use it
% to predict propeller performance at the bollard pull condition (J = 0).

% TO DO:
% - Add drag coefficient dependence on Reynolds number, neglecting induced
%   flow but still including wake fraction distribution.
% - Add more robust way of neglecting advance ratios where the efficiency
%   is less than zero.

% INPUTS:
% - J     : double
%           Advance ratio.
% - Re    : double
%           Reynolds number.
% - D     : double
%           Propeller diameter (m).
% - Z     : integer
%           Number of blades.
% - BAR   : double
%           Blade area ratio.
% - desig : string
%           Section shape name. Valid entries are "WB" (Wageningen
%           B-series), "NACA" (NACA series), or "C" (custom sections).
% - x     : double
%           Array of radius fractions.
% - c     : double
%           Array of chord lengths at specified radius fractions (m).
% - tmax  : double
%           Array of maximum section thicknesses at specified radius
%           fractions (m).
% - m     : double
%           Array of section cambers at specified radius fractions (m).
% - beta  : double
%           Array of geometric pitch angles at specified radius fractions
%           (rad).

% OUTPUTS:
% - KT     : double
%            Thrust coefficient for the entire propeller.
% - KQ     : double
%            Torque coefficient for the entire propeller.
% - eta    : double
%            Propeller open-water efficiency.
% - params : double
%            Solution parameters at each x-value. Includes: a, a', dKT/dx,
%            dKQ/dx, alpha, Cl, Cd, n_alpha, n_eta.

% SET UP DEFINITIONS
dKTs = zeros(1,length(x));                                                  % Thrust coefficients
dKQs = zeros(1,length(x));                                                  % Torque coefficients
etas = zeros(1,length(x));                                                  % Section efficiencies
tol = 1e-5;                                                                 % Iteration tolerance
maxit = 1000;                                                               % Maximum iterations

% PARAMETER ARRAYS
param_a = zeros(1,length(x));
param_a1 = zeros(1,length(x));
param_alpha = zeros(1,length(x));
param_Cl = zeros(1,length(x));
param_Cd = zeros(1,length(x));
param_nalpha = zeros(1,length(x));
param_neta = zeros(1,length(x));

% LOOP THROUGH EACH RADIUS FRACTION
for i=1:length(x)
    X = x(i);                                                               % Current radius fraction
    phialpha = beta(i);                                                     % Current pitch angle
    psi = atan2(J,pi*X);                                                    % Undisturbed inflow angle
    alpha = 0;                                                              % Guess for angle of attack
    for j=1:maxit
        phi = phialpha - alpha;                                             % Flow angle
        etai = tan(psi) / tan(phi);                                         % Ideal efficiency
        % Initially assume efficiency is ideal
        eta = etai;
        % Assume no drag initially
        gamma = 0;
        K = Goldstein(X,phi,Z);                                             % Goldstein correction factor
        for k=1:maxit
            a = (1-etai) / (etai+((tan(psi)^2)/eta));                       % Axial induction factor
            a1 = 1 - (etai*(1+a));                                          % Tangential induction factor
            % Compute change of KT over radius
            dKTdx = pi * (J^2) * X * K * a * (1+a);
            % Compute required Cl for this dKTdx
            Cl = dKTdx / (((pi^2)/4)*((Z*c(i))/D)*(X^2)*((1-a1)^2)*sec(phi)*(1-(tan(phi)*tan(gamma))));
            % Find drag coefficient and lift-camber slope from angle of
            % attack and radius fraction
            [Cd,dCldmc] = LD(alpha,Re,desig,X,tmax(i),c(i),[]);
            % If the required lift coefficient is small enough, then Cd/Cl
            % will grow without bound.
            if Cl < 1e-12
                Cl = 1e-12;
            end
            % Calculate drag angle
            gamma = atan2(Cd,Cl);
            eta_new = tan(psi) / tan(phi+gamma);                            % New efficiency
            % If the efficiency has converged then break out of the loop
            if abs(eta_new-eta) <= tol
                eta = eta_new;
                break;
            else
                if k == maxit
                    % If efficiency fails to converge, then break out of
                    % the loop
                    error("Inner (efficiency) loop failed to converge after "+num2str(maxit)+" iterations!"+newline + ...
                          "J = "+num2str(J)+", x = "+num2str(X)+", res = "+num2str(abs(eta-eta_new)))
                end
                eta = eta_new;
                continue;
            end
        end
        alpha_req = Req_AoA(Cl,J,X,desig,m(i),c(i));                          % Required angle of attack
        % Apply Ludwieg-Ginzel curvature correction
        [k1,k2] = Ludwieg_Ginzel(X,phi,BAR);
        dClda = (2*pi) * (180/pi);                                          % Lift-curve slope (1/deg)
        % Camber required for calculated Cl at 0 angle of attack is:
        mc0 = k1 * k2 * (Cl/dCldmc);
        dmc = mc0 - (m(i)/c(i));                                            % Camber deficit
        % Update required angle of attack to take flow curvature into
        % account
        alpha_req = alpha_req + (dmc*(dCldmc/dClda));
        % If the required angle of attack (based on the converged lift
        % coefficient) is the same as the iterated angle of attack, then
        % break out of the loop
        if abs(alpha_req-alpha) <= tol
            break;
        else
            if j == maxit
                % If angle of attack fails to converge, then break out of
                % the loop
                error("Outer (angle of attack) loop failed to converge after "+num2str(maxit)+" iterations!"+newline + ...
                      "J = "+num2str(J)+", x = "+num2str(X)+", res = "+num2str(abs(alpha_req-alpha)))
            end
            alpha = alpha_req;
            continue;
        end
    end
    % Set the current values for this x
    dKTs(i) = dKTdx;
    dKQs(i) = ((pi^2)/8) * ((Z*c(i))/D) * Cl * (X^3) * ((1-a1)^2) * sec(phi) * (tan(phi)+tan(gamma));
    etas(i) = eta;
    % Collect parameter data
    param_a(i) = a;
    param_a1(i) = a1;
    param_alpha(i) = alpha;
    param_Cl(i) = Cl;
    param_Cd(i) = Cd;
    param_nalpha(i) = j;
    param_neta(i) = k;
end

% FIND TOTAL KT, KQ, ETA BY INTEGRATION
% Find thrust and torque coefficients and efficiency over the entire
% propeller by the trapezium rule.
KT = 0; KQ = 0;
for i=2:length(x)
    dx = x(i) - x(i-1);
    KT = KT + (0.5*dx*(dKTs(i)+dKTs(i-1)));
    KQ = KQ + (0.5*dx*(dKQs(i)+dKQs(i-1)));
end
eta = (KT*J) / (2*pi*KQ);
% If a range of J is specified that goes beyond the point at which
% effiency vanishes, then the method converges to extremely large thrust
% and torque coeffients. This generally results in an open-water efficiency
% greater than 1, which is obviously unphysical. Thus, we can effectively
% ignore any advance coefficients which result in efficiencies above 1.
if eta >= 1 || eta < 0
    KT = nan;
    KQ = nan;
    eta = nan;
end

params = [param_a;param_a1;dKTs;dKQs;param_alpha;param_Cl;param_Cd;param_nalpha;param_neta]';

end
