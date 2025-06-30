function [c,tmax,m,coord,beta] = Wageningen(x,D,Z,BAR,P_D,TE,LE,N)

% Given details about the specified propeller, this function selects the
% corresponding Wageningen B-screw series propeller and produces the chord,
% thickness, and face pitch angle distributions as functions of the
% non-dimensionalised radial coordinate x, as well as coordinate arrays of
% each hydrofoil section to analyse for cavitation.

% INPUTS:
% - x   : double
%         Array of non-dimensionalised radial coordinates. Must be lower
%         than 1.0 (obviously) and greater than 0.2 (to account for the
%         hub).
% - D   : double
%         Propeller diameter (m).
% - Z   : integer
%         Number of blades. For Wageningen propellers must be between 3 and
%         7.
% - BAR : double
%         Blade-area ratio (ratio between expanded propeller area and area
%         of equivalent circle). Refer to Carlton, J. (2019) Marine
%         Propellers and Propulsion, p.99 for limits at different blade
%         numbers.
% - P_D : double
%         Geometric pitch ratio. For Wageningen B-series propellers, this
%         must be between 0.6 and 1.4.
% - TE  : double
%         Thickness of trailing edge for manufacturing purposes. NOTE: If
%         an initial cavitation check is required, this must be set to
%         zero.
% - LE  : double
%         Thickness of leading edge for manufacturing purposes.
% - N   : integer
%         Number of panels to use for each hydrofoil discretisation. Should
%         be even.

% OUTPUTS:
% - c     : double
%           Array of chord lengths (m).
% - tmax  : double
%           Array of maximum hydrofoil thicknesses (m).
% - m     : double
%           Array of maximum cambers (m).
% - coord : double
%           Tensor of coordinate definitions of each hydrofoil section (m).
% - beta  : double
%           Array of geometric face pitch angles (rad).

% INITIAL CHECKS
%  Check that there is an acceptable number of blades.
if Z < 3 || Z > 7
    error("You have specified a propeller with " + num2str(Z) + " blades, " + ...
          "which is outside the range of validity of the Wageningen series!")
end

% Check that the BAR is acceptable.
BAR_lims = [0.35,0.80;
            0.40,1.00;
            0.45,1.05;
            0.50,0.80;
            0.55,0.85];                                                     % Series BAR limits
if BAR < BAR_lims(Z-2,1) || BAR > BAR_lims(Z-2,2)
    warning("Data for " + num2str(Z) + "-bladed Wageningen propellers is " + ...
            "only available for " + num2str(BAR_lims(Z-2,1)) + " < BAR < " + ...
            num2str(BAR_lims(Z-2,2)) + ". Proceed with caution!")
end

% Check to make sure all the x-values are acceptable.
for i=1:length(x)
    if x(i) > 1.0 || x(i) < 0.2
        error("You have specified a radial coordinate outside the range of " + ...
              "validity! Valid non-dimensional radial coordinates are " + ...
              "between 0.2 and 1.0.")
    end
end

% Check that there is an even number of panels
if mod(N,2) ~= 0
    warning("You have specified an odd number of panels for each hydrofoil " + ...
            "discretisation! Converting to even N.")
    N = N + 1;
end

% GEOMETRIC PITCH ANGLE
if Z == 4
    % 4-bladed propellers have slightly non-constant pitch distributions
    % approaching the blade root
    rR = [0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0];             % Sample radius fractions
    k1 = [0.822 0.887 0.950 0.992 1.000 1.000 1.000 1.000 1.000];
    beta = zeros(length(x),1);                                              % Face pitch angle (rad)
    for i=1:length(x)
        beta(i) = atan2(P_D*interp1(rR,k1,x(i),"makima"),pi.*x(i));
    end
else
    beta = atan2(P_D,pi.*x);                                                % Face pitch angle (rad)
end

% PRINCIPAL DIMENSIONS
if Z == 3
    % Blade outline is slightly different for 3-bladed propellers
    rR = [0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0];                             % Sample x-values
    cs = [1.633 1.832 2.000 2.120 2.186 2.168 2.127 1.657 0.000];
    as = [0.616 0.611 0.599 0.583 0.558 0.526 0.481 0.400 0.000];
    bs = [0.350 0.350 0.350 0.355 0.389 0.442 0.478 0.500 0.000];
    Ar = [0.0526 0.0464 0.0402 0.0340 0.0278 0.0216 0.0154 0.0092 0.0030];
    Br = [0.0040 0.0035 0.0030 0.0025 0.0020 0.0015 0.0010 0.0005 0.0000];
    % Interpolate based on queried x-values
    c = ((D*BAR)/Z) .* interp1(rR,cs,x,"makima");                           % Chord lengths (m)
    a = c .* interp1(rR,as,x,"makima");                                     % Directrix offset from LE (m)
    b = c .* interp1(rR,bs,x,"makima");                                     % tmax offset from LE (m)
    tmax = D .* (interp1(rR,Ar,x,"makima")-(Z.*interp1(rR,Br,x,"makima"))); % Max. thickness (m).
else
    rR = [0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0];
    cs = [1.662 1.882 2.050 2.152 2.187 2.144 1.970 1.582 0.000];
    as = [0.617 0.613 0.601 0.586 0.561 0.524 0.463 0.351 0.000];
    bs = [0.350 0.350 0.351 0.355 0.389 0.443 0.479 0.500 0.000];
    Ar = [0.0526 0.0464 0.0402 0.0340 0.0278 0.0216 0.0154 0.0092 0.0030];
    Br = [0.0040 0.0035 0.0030 0.0025 0.0020 0.0015 0.0010 0.0005 0.0000];
    % Interpolate based on queried x-values
    c = ((D*BAR)/Z) .* interp1(rR,cs,x,"makima");                           % Chord lengths (m)
    a = c .* interp1(rR,as,x,"makima");                                     % Directrix offset from LE (m)
    b = c .* interp1(rR,bs,x,"makima");                                     % tmax offset from LE (m)
    tmax = D .* (interp1(rR,Ar,x,"makima")-(Z.*interp1(rR,Br,x,"makima"))); % Max. thickness (m)
end

% HYDROFOIL COORDINATES
coord = zeros(N+1,2,length(x));                                             % Empty coordinate tensor
m = zeros(length(x),1);                                                     % Empty camber array
rR = [1.0 0.95 0.9 0.85 0.8 0.7 0.6 0.5 0.4 0.3 0.25 0.2 0.15];
Ps = [-1.0   -0.95  -0.9   -0.8   -0.7   -0.6   -0.5   -0.4   -0.2   0      0.2    0.4    0.5    0.6    0.7    0.8    0.85   0.9    0.95   1.0   ];
V1 = [0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000;
      0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000;
      0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000;
      0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000;
      0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000;
      0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000;
      0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0006 0.0022 0.0067 0.0169 0.0382;
      0.0522 0.0420 0.0330 0.0190 0.0100 0.0040 0.0012 0.0000 0.0000 0.0000 0.0000 0.0000 0.0008 0.0034 0.0085 0.0211 0.0328 0.0500 0.0778 0.1278;
      0.1467 0.1200 0.0972 0.0630 0.0395 0.0214 0.0116 0.0044 0.0000 0.0000 0.0000 0.0033 0.0090 0.0189 0.0357 0.0637 0.0833 0.1088 0.1467 0.2181;
      0.2306 0.2040 0.1790 0.1333 0.0943 0.0623 0.0376 0.0202 0.0033 0.0000 0.0027 0.0148 0.0300 0.0503 0.0790 0.1191 0.1445 0.1760 0.2186 0.2923;
      0.2598 0.2372 0.2115 0.1651 0.1246 0.0899 0.0579 0.0350 0.0084 0.0000 0.0031 0.0224 0.0417 0.0669 0.1008 0.1465 0.1747 0.2068 0.2513 0.3256;
      0.2826 0.2630 0.2400 0.1967 0.1570 0.1207 0.0880 0.0592 0.0172 0.0000 0.0049 0.0304 0.0520 0.0804 0.1180 0.1685 0.2000 0.2353 0.2821 0.3560;
      0.3000 0.2824 0.2650 0.2300 0.1950 0.1610 0.1280 0.0955 0.0365 0.0000 0.0096 0.0384 0.0615 0.0920 0.1320 0.1870 0.2230 0.2642 0.3150 0.3860];
V2 = [0.0000 0.0975 0.1900 0.3600 0.5100 0.6400 0.7500 0.8400 0.9600 1.0000 0.9600 0.8400 0.7500 0.6400 0.5100 0.3600 0.2775 0.1900 0.0975 0.0000;
      0.0000 0.0975 0.1900 0.3600 0.5100 0.6400 0.7500 0.8400 0.9600 1.0000 0.9600 0.8400 0.7500 0.6400 0.5100 0.3600 0.2775 0.1900 0.0975 0.0000;
      0.0000 0.0975 0.1900 0.3600 0.5100 0.6400 0.7500 0.8400 0.9600 1.0000 0.9600 0.8400 0.7500 0.6400 0.5100 0.3600 0.2775 0.1900 0.0975 0.0000;
      0.0000 0.0975 0.1900 0.3600 0.5100 0.6400 0.7500 0.8400 0.9600 1.0000 0.9615 0.8450 0.7550 0.6455 0.5160 0.3660 0.2830 0.1950 0.1000 0.0000;
      0.0000 0.0975 0.1900 0.3600 0.5100 0.6400 0.7500 0.8400 0.9600 1.0000 0.9635 0.8520 0.7635 0.6545 0.5265 0.3765 0.2925 0.2028 0.1050 0.0000;
      0.0000 0.0975 0.1900 0.3600 0.5100 0.6400 0.7500 0.8400 0.9600 1.0000 0.9675 0.8660 0.7850 0.6840 0.5615 0.4140 0.3300 0.2337 0.1240 0.0000;
      0.0000 0.0965 0.1885 0.3585 0.5110 0.6415 0.7530 0.8426 0.9613 1.0000 0.9690 0.8790 0.8090 0.7200 0.6060 0.4620 0.3775 0.2720 0.1485 0.0000;
      0.0000 0.0950 0.1865 0.3569 0.5140 0.6439 0.7580 0.8456 0.9639 1.0000 0.9710 0.8880 0.8275 0.7478 0.6430 0.5039 0.4135 0.3056 0.1750 0.0000;
      0.0000 0.0905 0.1810 0.3500 0.5040 0.6353 0.7525 0.8415 0.9645 1.0000 0.9725 0.8933 0.8345 0.7593 0.6590 0.5220 0.4335 0.3235 0.1935 0.0000;
      0.0000 0.0800 0.1670 0.3360 0.4885 0.6195 0.7335 0.8265 0.9583 1.0000 0.9750 0.8920 0.8315 0.7520 0.6505 0.5130 0.4265 0.3197 0.1890 0.0000;
      0.0000 0.0725 0.1567 0.3228 0.4740 0.6050 0.7184 0.8139 0.9519 1.0000 0.9751 0.8899 0.8259 0.7415 0.6359 0.4982 0.4108 0.3042 0.1758 0.0000;
      0.0000 0.0640 0.1455 0.3060 0.4535 0.5842 0.6995 0.7984 0.9446 1.0000 0.9750 0.8875 0.8170 0.7277 0.6190 0.4777 0.3905 0.2840 0.1560 0.0000;
      0.0000 0.0540 0.1325 0.2870 0.4280 0.5585 0.6770 0.7805 0.9360 1.0000 0.9760 0.8825 0.8055 0.7105 0.5995 0.4520 0.3665 0.2600 0.1300 0.0000];
% Apply cosine spacing
theta = linspace(0,pi,(N/2)+1);
for i=1:length(x)
    s = 0.5 .* c(i) .* (1-cos(theta));                                      % Chordwise ordinate
    s = s - (c(i)-a(i));
    y_face = zeros(length(s),1);
    y_back = zeros(length(s),1);
    for j=1:length(s)
        if s(j) > a(i)-b(i)
            P = (s(j)-(a(i)-b(i))) / b(i);
            v1 = interp2(Ps,rR,V1,P,x(i),"linear");
            v2 = interp2(Ps,rR,V2,P,x(i),"linear");
            y_face(j) = v1 * (tmax(i)-LE);
            y_back(j) = ((v1+v2)*(tmax(i)-LE)) + LE;
        else
            P = (s(j)-(a(i)-b(i))) / (c(i)-b(i));
            v1 = interp2(Ps,rR,V1,P,x(i),"linear");
            v2 = interp2(Ps,rR,V2,P,x(i),"linear");
            y_face(j) = v1 * (tmax(i)-TE);
            y_back(j) = ((v1+v2)*(tmax(i)-TE)) + TE;
        end
    end
    camber = 0.5 .* (y_back+y_face);                                        % Camber along hydrofoil
    m(i) = max(camber);                                                     % Maximum camber (m)
    % "Wrap" coordinates around hydrofoil
    coord(:,1,i) = [s(1:end-1) flip(s)];
    coord(:,2,i) = [y_back(1:end-1)' flip(y_face)'];
end

% PLOT PROPELLER EXPANDED OUTLINE AND PITCHED SECTIONS
figure(1)
subplot(1,2,1)
desig = "B" + num2str(Z) + "." + num2str(100*BAR);
for i=1:length(x)
    plot(coord(:,1,i),coord(:,2,i)+(0.5*x(i)*D),"k-"), hold on
end
hold off
axis equal
title("Wageningen " + desig + " Propeller Expanded Outline")
xlabel("x (m)")
ylabel("y (m)")
grid()

subplot(1,2,2)
for i=1:length(x)
    rot = [[cos(beta(i)) -sin(beta(i))];
           [sin(beta(i))  cos(beta(i))]];                                   % Rotation matrix
    xy = zeros(size(coord(:,:,1)'));                                        % Rotated section
    for j=1:length(coord(:,1,1))
        % Rotate coordinates by face pitch angle
        xy(:,j) = rot * coord(j,:,i)';
    end
    plot(xy(1,:),xy(2,:),"k-"), hold on
end
hold off
axis equal
title("Pitched Hydrofoil Sections")
xlabel("x (m)")
ylabel("y (m)")
grid()

% DRAW PROPELLER IN 3D VIEW
figure(2)
% Loop through each radius fraction, pitching each section and bending the
% coordinates around an appropriate cylinder.
xyz = zeros(length(coord(:,1,1)),3,length(x));
for i=1:length(x)
    % Define transformation matrices
    r = 0.5 * x(i) * D;                                                     % Current radius (m)
    R = [[cos(beta(i)) -sin(beta(i))];
         [sin(beta(i))  cos(beta(i))]];                                     % Pitch rotation matrix
    W = [[1/r 0];
         [0   1]];                                                          % Bending matrix
    % Loop through each point specified in coord
    for j=1:length(coord(:,1,1))
        TZ = W * R * coord(j,:,i)';                                         % Apply transformations
        theta = TZ(1);
        z = TZ(2);
        [x1,x2,x3] = pol2cart(theta,r,z);
        xyz(j,:,i) = [x1,x2,x3]';
    end
end
for J=1:length(x)-1
    for I=1:length(coord(:,1,1))-1
        ns = [xyz(I,:,J);
              xyz(I+1,:,J);
              xyz(I+1,:,J+1);
              xyz(I,:,J+1)];
        assignin("base","ns",ns);
        fill3(ns(:,1),ns(:,2),ns(:,3),"w"), hold on
    end
end
hold off
axis equal
xlabel("x")
ylabel("y")
zlabel("z")

end
