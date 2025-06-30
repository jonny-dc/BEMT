function [KT,KQ,eta] = Wageningen_KTKQ(J,P_D,BAR,Z,Re,graph)

% Finds the thrust coefficient, torque coefficient, and efficiency for a
% Wageningen B-series propeller of pitch P/D and specified blade area
% ratio, over a range of advance ratios. Based on faired data taken from
% Oosterveld & van Oosannen, (1975) Further Computer-Analysed Data of the
% Wageningen B-Screw Series.

% INPUTS:
% - J     : double
%           Array of queried advance ratios.
% - P_D   : double
%           Pitch-diameter ratio.
% - BAR   : double
%           Expanded blade area ratio.
% - Z     : integer
%           Number of blades.
% - Re    : double
%           Reynolds number.
% - graph : string
%           Whether to plot KT-KQ-eta-J graphs, valid entries are "Y" or
%           "N".

% OUTPUTS:
% - KT  : double
%         Array of thrust coefficients at the specified advance ratios.
% - KQ  : double
%         Array of torque coefficients at the specified advance ratios.
% - eta : double
%         Array of propeller efficiencies at the specified advance ratios.

% Make advance ratios into correct shape
J = J';

datT = [[ 0.00880496   0 0 0 0];
        [-0.204554     1 0 0 0];
        [ 0.166351     0 1 0 0];
        [ 0.158114     0 2 0 0];
        [-0.147581     2 0 1 0];
        [-0.481497     1 1 1 0];
        [ 0.415437     0 2 1 0];
        [ 0.0144043    0 0 0 1];
        [-0.0530054    2 0 0 1];
        [ 0.0143481    0 1 0 1];
        [ 0.0606826    1 1 0 1];
        [-0.0125894    0 0 1 1];
        [ 0.0109689    1 0 1 1];
        [-0.133698     0 3 0 0];
        [ 0.00638407   0 6 0 0];
        [-0.00132718   2 6 0 0];
        [ 0.168496     3 0 1 0];
        [-0.0507214    0 0 2 0];
        [ 0.0854559    2 0 2 0];
        [-0.0504475    3 0 2 0];
        [ 0.0104650    1 6 2 0];
        [-0.00648272   2 6 2 0];
        [-0.00841728   0 3 0 1];
        [ 0.0168424    1 3 0 1];
        [-0.00102296   3 3 0 1];
        [-0.0317791    0 3 1 1];
        [ 0.0186040    1 0 2 1];
        [-0.00410798   0 2 2 1];
        [-0.000606848  0 0 0 2];
        [-0.00498190   1 0 0 2];
        [ 0.00259830   2 0 0 2];
        [-0.000560528  3 0 0 2];
        [-0.00163652   1 2 0 2];
        [-0.000328787  1 6 0 2];
        [ 0.000116502  2 6 0 2];
        [ 0.000690904  0 0 1 2];
        [ 0.00421749   0 3 1 2];
        [ 0.0000565229 3 6 1 2];
        [-0.00146564   0 3 2 2]];                                           % KT factors

datQ = [[ 0.00379368   0 0 0 0];
        [ 0.00886523   2 0 0 0];
        [-0.0322410    1 1 0 0];
        [ 0.00344778   0 2 0 0];
        [-0.0408811    0 1 1 0];
        [-0.108009     1 1 1 0];
        [-0.0885381    2 1 1 0];
        [ 0.188561     0 2 1 0];
        [-0.00370871   1 0 0 1];
        [ 0.00513696   0 1 0 1];
        [ 0.0209449    1 1 0 1];
        [ 0.00474319   2 1 0 1];
        [-0.00723408   2 0 1 1];
        [ 0.00438388   1 1 1 1];
        [-0.0269403    0 2 1 1];
        [ 0.0558082    3 0 1 0];
        [ 0.0161886    0 3 1 0];
        [ 0.00318086   1 3 1 0];
        [ 0.0158960    0 0 2 0];
        [ 0.0471729    1 0 2 0];
        [ 0.0196283    3 0 2 0];
        [-0.0502782    0 1 2 0];
        [-0.0300550    3 1 2 0];
        [ 0.0417122    2 2 2 0];
        [-0.0397722    0 3 2 0];
        [-0.00350024   0 6 2 0];
        [-0.0106854    3 0 0 1];
        [ 0.00110903   3 3 0 1];
        [-0.000313912  0 6 0 1];
        [ 0.00359850   3 0 1 1];
        [-0.00142121   0 6 1 1];
        [-0.00383637   1 0 2 1];
        [ 0.0126803    0 2 2 1];
        [-0.00318278   2 3 2 1];
        [ 0.00334268   0 6 2 1];
        [-0.00183491   1 1 0 2];
        [ 0.000112451  3 2 0 2];
        [-0.0000297228 3 6 0 2];
        [ 0.000269551  1 0 1 2];
        [ 0.000832650  2 0 1 2];
        [ 0.00155334   0 2 1 2];
        [ 0.000302683  0 6 1 2];
        [-0.000184300  0 0 2 2];
        [-0.000425399  0 3 2 2];
        [ 0.0000869243 3 3 2 2];
        [-0.000465900  0 6 2 2];
        [ 0.0000554194 1 6 2 2]];                                           % KQ factors

% Reynolds number corrections
dKT = 0.000353485 ...
      - (0.00333758.*BAR.*(J.^2)) ...
      - (0.00478125.*BAR.*P_D.*J) ...
      + (0.000257792.*((log10(Re)-0.301).^2).*BAR.*(J.^2)) ...
      + (0.0000643192.*(log10(Re)-0.301).*(P_D.^6).*(J.^2)) ...
      - (0.0000110636.*((log10(Re)-0.301).^2).*(P_D.^6).*(J.^2)) ...
      - (0.0000276305.*(log10(Re)-0.301).*Z.*BAR.*(J.^2)) ...
      + (0.0000954000.*(log10(Re)-0.301).*Z.*BAR.*P_D.*J) ...
      + (0.0000032049.*(log10(Re)-0.301).*(Z.^2).*BAR.*(P_D.^3).*J);

dKQ = -0.000591412 ...
      + (0.00696898.*P_D) ...
      - (0.0000666654.*Z.*(P_D.^6)) ...
      + (0.0160818.*(BAR.^2)) ...
      - (0.000938091.*(log10(Re)-0.301).*P_D) ...
      - (0.000595930.*(log10(Re)-0.301).*(P_D.^2)) ...
      + (0.0000782099.*((log10(Re)-0.301).^2).*(P_D.^2)) ...
      + (0.0000052199.*(log10(Re)-0.301).*Z.*BAR.*(J.^2)) ...
      - (0.00000088528.*((log10(Re)-0.301).^2).*Z.*BAR.*P_D.*J) ...
      + (0.0000230171.*(log10(Re)-0.301).*Z.*(P_D.^6)) ...
      - (0.00000184341.*((log10(Re)-0.301).^2).*Z.*(P_D.^6)) ...
      - (0.00400252.*(log10(Re)-0.301).*(BAR.^2)) ...
      + (0.000220915.*((log10(Re)-0.301).^2).*(BAR.^2));

KT = zeros(length(J),1);                                                    % Thrust coefficient
KQ = zeros(length(J),1);                                                    % Torque coefficient
eta = zeros(length(J),1);                                                   % Efficiency

% Calculate coefficients at specified Reynolds number
for i=1:length(J)
    KT(i) = sum(datT(:,1).*(J(i).^datT(:,2)).*(P_D.^datT(:,3)).*(BAR.^datT(:,4)).*(Z.^datT(:,5)));
    KT(i) = KT(i) + dKT(i);
    KQ(i) = sum(datQ(:,1).*(J(i).^datQ(:,2)).*(P_D.^datQ(:,3)).*(BAR.^datQ(:,4)).*(Z.^datQ(:,5)));
    KQ(i) = KQ(i) + dKQ(i);
    eta(i) = (KT(i)*J(i)) / (2*pi*KQ(i));
    % If the efficiency drops below zero then the curves should stop
    if eta(i) < 0
        KT(i:end) = nan;
        KQ(i:end) = nan;
        eta(i:end) = nan;
        break;
    end
end

% Plot graphs, if requested
if graph == "Y"
    figure()
    plot(J,KT,"b-"), hold on
    plot(J,10.*KQ,"r-"), hold on
    plot(J,eta,"k-"), hold off
    grid()
    title("Propeller Performance for Wageningen B-Series")
    xlabel("Advance ratio J")
    ylabel("Thrust coeff. K_T, Torque coeff 10K_Q, Efficiency \eta")
    legend("K_T","10\times K_Q","\eta")
end

end
