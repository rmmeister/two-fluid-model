function [ F ] = calcResiduals(  Y, alphaGvG )
%CALCRESIDUALS Calulate the residuals vector, F.
%   Residual terms are obtained using well data and vector Y containing
%   parameters to be defined.

% INPUT DATA
pipeDiameter = .032; % m
gravityAcceleration = 9.8; % m/s2
deviationAngle = degtorad(0);
gasViscosity = 1.78e-5; % Kg/m/s
liquidViscosity = 9.98e-4; % Kg/m/s
gasDensity = 1.725; % Kg/m3
liquidDensity = 998; % Kg/m3
liquidMassRate = 80; % Liquid mass flow-rate, Kg/s/m2
gasMassRate = alphaGvG*gasDensity; % gass mass flow-rate, kg/s/m2

% Parameters to be found
alphaG = Y(1); alphaL = Y(2); vG = Y(3); vL = Y(4);

% calculate frictional CLOSURES
[ RGW, RLW, RGL ] = calcTPM( alphaG, alphaL, vG, vL, ...
gasDensity, liquidDensity, gasViscosity, liquidViscosity, pipeDiameter, deviationAngle );

% obtain RESIDUAL terms
f1 = alphaG*gasDensity*vG - gasMassRate;
f2 = alphaL*liquidDensity*vL - liquidMassRate;
f3 = RGL/alphaG - RGW/alphaG + gasDensity*gravityAcceleration*sin(deviationAngle) + RGL/alphaL + ...
RLW/alphaL - liquidDensity*gravityAcceleration*sin(deviationAngle);
f4 = alphaG + alphaL - 1;

% vector of RESIDUALS
F = [f1; f2; f3; f4];
end

