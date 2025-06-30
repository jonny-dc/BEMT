% main.m

% This program implements a blade-element momentum theory routine, as
% outlined by Molland, Turnock & Hudson (2011), to predict the performance
% of a given marine propeller geometry over a range of advance ratios.

% Depends on:
% - BEMT_Solver.m : Iteratively solve the blade-element momentum equations
%   for each blade section and advance ratio.
% - Wageningen.m : Generate the required geometric information for a
%   specific Wageningen B-series propeller.
% - Wageningen_KTKQ.m : Generate KT, KQ and efficiency curves for the
%   selected propeller, using faired experimental data.
% - Goldstein.m : Find the Goldstein correction for a given radius fraction
%   and inflow angle.
% - Ludwieg_Ginzel.m : Find the Ludwieg-Ginzel curvature corrections for a
%   given radius fraction.
% - Req_AoA.m : Find the required angle of attack corresponding to a given
%   lift coefficient.

% GENERATE PROPELLER GEOMETRY
x = [0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95]; % Radius fractions
D = 7.9;                                                                    % Diameter (m)
BAR = 0.75;                                                                 % Blade area ratio
Z = 5;                                                                      % Number of blades
P_D = 0.6;                                                                  % Pitch ratio
TE = 0.0;   LE = 0.0;                                                       % Edge thicknesses (m)
N = 200;                                                                    % Panel number
% Generate Wageningen B propeller geometry
[c,tmax,m,coord,beta] = Wageningen(x,D,Z,BAR,P_D,TE,LE,N);

% Define advance ratios, average Reynolds number, ship speed, RPM and wake
% fraction distribution
J = linspace(0.05,1.45,35);
Re = 2e6;                                                                   % Average Reynolds number
rho = 1025;                                                                 % Fluid density (kg/m3)
nu = 1.14e-6;                                                               % Fluid viscosity (m2/s)
h = 0.87 * D;                                                               % Cavitation depth (m)
p_atm = 101325;                                                             % Atmospheric pressure (Pa)
p_v = 1705.6;                                                               % Vapour pressure (Pa)
g = 9.81;                                                                   % Gravitational acceleration (m/s2)
% For the moment, assume ship speed is kept constant and RPM can therefore
% be calculated from J
Vs = 10;                                                                    % Ship speed (m/s)
% Nominal wake fraction distribution, circumferentially averaged
wf = [0.522987 0.410575 0.312838 0.229777 0.161393 0.107684 0.068652 0.044295];
wt = sum(wf) / length(wf);                                                  % Average wake fraction
Va = Vs * (1-wt);                                                           % Advance velocity (m/s)
n = Va ./ (J.*D);                                                           % Revolutions per second

% GENERATE REFERENCE PERFORMANCE CURVES
[KT_ref,KQ_ref,eta_ref] = Wageningen_KTKQ(J,P_D,BAR,Z,Re,"N");

% SOLVE BEMT METHOD
print = true;
KT = zeros(length(J),1);
KQ = zeros(length(J),1);
eta = zeros(length(J),1);
params = zeros(length(x),9,length(J));
cav = zeros(length(x),length(J));
for i=1:length(J)
    [KT(i),KQ(i),eta(i),params(:,:,i)] = BEMT_Solver(J(i),Re,D,Z,BAR,"WB",x,c,tmax,m,beta);
    Vr = sqrt((Va.^2)+((pi.*n(i).*(x.*D)).^2));                             % Section velocity (m/s)
    Re_x = (Vr.*c) / nu;                                                    % Section Reynolds number
    sigma = (p_atm+(rho.*g.*h)-p_v) ./ (0.5.*rho.*(Vr.^2));                 % Cavitation number
    Cl = params(:,7,i);
    % Perform very rough initial cavitation check using approximate Gutsche
    % bucket diagram
    for j=1:length(x)
        cav(j,i) = Cav_Bucket_Check(sigma(j),Cl(j),tmax(j),c(j));
    end
    % Display solution data for each J to the console?
    if print == true
        if ~isnan(KT(i))
            disp("J = "+num2str(J(i)))
            print = [x;params(:,:,i)']';
            disp("    x        |a        |a'       |dKT/dx   |dKQ/dx   |alpha    |Cl       |Cd       |n_alpha  |n_eta    ")
            disp("    ---------------------------------------------------------------------------------------------------")
            disp(print)
        end
    end
end

% PLOT PERFORMANCE CURVES
figure()
plot(J,KT_ref,"b-",J,10.*KQ_ref,"r-",J,eta_ref,"k-"), hold on
plot(J,KT,"bo--",J,10.*KQ,"ro--",J,eta,"ko--"), hold off
grid on, grid minor
title("Wageningen B"+num2str(Z)+"-"+num2str(BAR*100)+" Propeller")
xlabel("Advance ratio J")
ylabel("Thrust coefficient K_T, Torque coefficient 10\times K_Q, Efficiency \eta")
legend("K_T","10\times K_Q","eta")
