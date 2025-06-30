function cav = Cav_Bucket_Check(sigma,Cl,tmax,c)

% Function to carry out an initial cavitation check for each propeller
% section, using a Gutsche-type bucket diagram. 

% INPUTS:
% - sigma : double
%           Cavitation number.
% - Cl    : double
%           Lift coefficient.
% - tmax  : double
%           Maximum section thickness (m).
% - c     : double
%           Section chord length (m).

% OUTPUTS:
% - cav : double
%         Cavitation indicator. Is 0 for no cavitation, and 1 for a
%         cavitating section.

% DEFINITIONS
Cli = 0.25;                                                                 % Ideal lift coefficient
k = 0.15;                                                                   % Nose radius factor
t_c = tmax / c;                                                             % Thickness ratio
r_c = k * (t_c^2);                                                          % Nose radius ratio
cav = 0;

% First check back bubble cavitation
if sigma < ((2/3)*Cl) + ((5/2)*t_c)
    cav = 1;
end

% Now check sheet cavitation
if Cl < Cli
    % Face sheet cavitation, lift coefficient less than ideal
    if sigma < (0.06*((Cl-Cli)^2)) / r_c
        cav = 1;
    end
elseif Cl > Cli
    % Back sheet cavitation, lift coefficient more than ideal
    if sigma < (0.06*((Cli-Cl)^2)) / r_c
        cav = 1;
    end
end

end
