% Test the accuracy of the Wageningen performance curve generator
P_D = [0.6 0.8 1.0 1.2 1.4];
J = linspace(0,1.5,300);
Z = 5;
BAR = 0.75;
Re = 2e6;
graph = "N";
figure()

for i=1:length(P_D)
    [KT,KQ,eta] = Wageningen_KTKQ(J,P_D(i),BAR,Z,Re,graph);
    plot(J,KT,"b-"), hold on
    plot(J,10.*KQ,"r-"), hold on
    plot(J,eta,"k-"), hold on
end

hold off
title("Wageningen B5-75 Propeller Curves")
grid on, grid minor
axis equal
xlabel("Advance ratio J")
ylabel("Thrust coeff. K_T, Torque coeff. 10\times K_Q, Efficiency \eta")
legend("K_T","10\times K_Q","\eta")

% Test the accuracy of the curvature correction function
BAR = linspace(0.4,1.2,100);
x1 = [0.2 0.3 0.4 0.5 0.6];
phi1 = atan(1);
x2 = [0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9];
phi2 = 0.0;
% Find and plot k1
figure()
subplot(2,1,2)
for i=1:length(x1)
    k1 = zeros(length(BAR),1);
    for j=1:length(BAR)
        [k1(j),k2] = Ludwieg_Ginzel(x1(i),phi1,BAR(j));
    end
    plot(BAR,k1,"k-"), hold on
end
hold off
grid on, grid minor
xlabel("Expanded Blade Area Ratio A_E/A_0")
ylabel("k_1")
% Find and plot k2
subplot(2,1,1)
for i=1:length(x2)
    k2 = zeros(length(BAR),1);
    for j=1:length(BAR)
        [k1,k2(j)] = Ludwieg_Ginzel(x2(i),phi2,BAR(j));
    end
    plot(BAR,k2,"k-"), hold on
end
hold off
title("Ludwieg-Ginzel Curvature Correction Factors")
grid on, grid minor
ylabel("k_2")
