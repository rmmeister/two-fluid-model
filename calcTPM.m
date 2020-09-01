function [ RGW, RLW, RGL, hL, hG, beta ] = calcTPM( alphaG, alphaL, vG, vL, rhoG, rhoL, muG, muL, d, theta )
%CALCTPM Calculate two-phase properties
% 

e = 0.0006;
g = 9.8;
A = pi/4*d^2;
% calculating parameters
beta = 2*pi - 2*(pi*alphaL + (3*pi/2)^(1/3)*(1 - 2*alphaL + alphaL^(1/3) ...
    - (1 - alphaL)^(1/3))); % rad
hG = (1/2*cos(pi-beta/2) + 1/3/pi/alphaG*(sin(pi-beta/2))^3)*d*cos(theta);
hL = (-1/2*cos(pi-beta/2) + 1/3/pi/alphaL*(sin(pi-beta/2))^3)*d*cos(theta);
OGW = beta*d/2;
OLW = (2*pi - beta)*d/2;
OGL = d*sin(beta/2);
AG = alphaG*pi*d^2/4;
AL = alphaL*pi*d^2/4;
d_hGW = 4*AG/OGW;
d_hLW = 4*AL/OLW;
Re_hGW = vG*d_hGW*rhoG/muG;
Re_hLW = vL*d_hLW*rhoL/muL;

% determine friction factor on gas liquid interface
Re_SL = rhoL*alphaL*vL*d/muL;
FrL = vL/sqrt(g*hL);
fGL = (.004+.5e-6*Re_SL)*FrL^1.335 * rhoL*d*g/rhoG/vG^2;

% solving Colebrook and White equation for phase friction factor
fGW = fsolve(@(fGW) -1/sqrt(fGW) -2*log10(2.51/Re_hGW/sqrt(fGW) + e/3.7/d_hGW) , 0.01,...
    optimoptions('fsolve','Display','off'));
fLW = fsolve(@(fLW) -1/sqrt(fLW) -2*log10(2.51/Re_hLW/sqrt(fLW) + e/3.7/d_hLW) , 0.01,...
    optimoptions('fsolve','Display','off'));

% friction
RGW = -fGW*rhoG*OGW*vG*abs(vG)/8/A;
RLW = -fLW*rhoL*OLW*vL*abs(vL)/8/A;
RGL = fGL*rhoG*OGL*(vG-vL)*abs(vG-vL)/8/A;
end

