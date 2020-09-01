%% INITIALIZE MATLAB
clear all
clc
close all
format long
%% Solving Y with Levenberg-Marquardt Optimization
Y0 = [0.6; 0.4; 10; 2];         % initial guess for unknowns vector
alphaGvG = 10:10:50;            % product of gas velocity and void fraction
Y = zeros(4, length(alphaGvG)); % unknowns vector

fprintf('Solving for vector of unknowns at %g points...\n', length(alphaGvG));
for i = 1:length(alphaGvG)
    [Y(:, i), fval, flag, output] = fsolve(@(Y) calcResiduals(Y, alphaGvG(i)), Y0, ...
        optimoptions('fsolve','Display','off'));
end

if flag >= 1
    fprintf('\nSimulation converged successfully.\n');
end
%% RESULTS
% INPUT DATA
pipeDiameter = .032; % m
deviationAngle = degtorad(0);
inletP = .15; % MPa
gasViscosity = 1.78e-5; % Kg/m/s
liquidViscosity = 9.98e-4; % Kg/m/s
gasDensity = 1.725; % Kg/m3
liquidDensity = 998; % Kg/m3
gravityConstant = 9.8;

% Calculating Pressure Gradient for the range of alphaG*vG
for i = 1:length(alphaGvG) 
[ RGW(i), RLW, RGL(i) ] = calcTPM( Y(1, i), Y(2, i), Y(3, i), Y(4, i), gasDensity, ...
                                            liquidDensity, gasViscosity, liquidViscosity, pipeDiameter, deviationAngle );
pressureGradient(i) = RGL(i)/Y(2, i) + RGW(i)/Y(2, i) - gasDensity*gravityConstant*sin(deviationAngle);
end

% PLOTING THE RESULTS
col = [0.850980401039124 0.325490206480026 0.0980392172932625];

figure('Position', [10 10 900 700]) 
subplot(221);
plot(alphaGvG, Y(1, :), 'Color', col ,'LineWidth', 3);
grid on;
xlabel('\alpha_G v_G, m/s');
ylabel('\alpha_G');
title('Void fraction');

subplot(222);
plot(alphaGvG, Y(2, :), 'Color', col, 'LineWidth', 3); grid on;
xlabel('\alpha_G v_G, m/s');
ylabel('\alpha_L');
title('Liquid Hold-up');

subplot(223);
plot(alphaGvG, Y(3, :), 'Color', col, 'LineWidth', 3); grid on;
xlabel('\alpha_G v_G, m/s');
ylabel('v_G, m/s');
title('Gas Velocity');

subplot(224);
plot(alphaGvG, Y(4, :), 'Color', col, 'LineWidth', 3); grid on;
xlabel('\alpha_G v_G, m/s');
ylabel('v_L, m/s');
title('Liquid Velocity');

figure
plot(alphaGvG, -pressureGradient*1e-6,'Color', col, 'LineWidth', 3); grid on;
xlabel('\alpha_G v_G, m/s');
ylabel('-dP/dx');
title('Pressure Gradient, MPa/m');
