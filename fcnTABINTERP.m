function [ valRPM, valFx, valFy, valMx, valMy,...
    valQ, valCP, valMUinf ] ...
    = fcnTABINTERP( vecRPM, alpha, tabLOOKUP, flowq, valTHRUSTrho, flowRHO, geomDIAMETER )
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

%flow q = 2.8

% interpolate for everything based on q and T
for i = 1:length(vecRPM)
    
    tempTHRUST_rho = tabLOOKUP.Thrust_rho(tabLOOKUP.alpha == alpha & tabLOOKUP.RPM == vecRPM(i));
    tempq = tabLOOKUP.q(tabLOOKUP.alpha == alpha & tabLOOKUP.RPM == vecRPM(i));
    tempFx_rho = tabLOOKUP.Px_rho(tabLOOKUP.alpha == alpha & tabLOOKUP.RPM == vecRPM(i));
    tempFy_rho = tabLOOKUP.Py_rho(tabLOOKUP.alpha == alpha & tabLOOKUP.RPM == vecRPM(i));
    tempMx_rho = tabLOOKUP.Mx_rho(tabLOOKUP.alpha == alpha & tabLOOKUP.RPM == vecRPM(i));
    tempMy_rho = tabLOOKUP.My_rho(tabLOOKUP.alpha == alpha & tabLOOKUP.RPM == vecRPM(i));
    tempQ_rho = tabLOOKUP.Q_rho(tabLOOKUP.alpha == alpha & tabLOOKUP.RPM == vecRPM(i));
    tempCP = tabLOOKUP.CP(tabLOOKUP.alpha == alpha & tabLOOKUP.RPM == vecRPM(i));
    tempMUinf = tabLOOKUP.mu_inf(tabLOOKUP.alpha == alpha & tabLOOKUP.RPM == vecRPM(i));
    
    vecTHRUST(i) = interp1(tempq,tempTHRUST_rho,flowq)';
    vecFx_rho(i) = interp1(tempq,tempFx_rho,flowq)';
    vecFy_rho(i) = interp1(tempq,tempFy_rho,flowq)';
    vecMx_rho(i) = interp1(tempq,tempMx_rho,flowq)';
    vecMy_rho(i) = interp1(tempq,tempMy_rho,flowq)';
    vecQ_rho(i) = interp1(tempq,tempQ_rho,flowq)';
    vecCP(i) = interp1(tempq,tempCP,flowq)';
    vecMUinf(i) = interp1(tempq,tempMUinf,flowq)';


end

tempCURVE=polyfit(vecTHRUST,vecRPM',2);

valRPM=tempCURVE(1)*valTHRUSTrho^2+tempCURVE(2)*valTHRUSTrho+tempCURVE(3);


Omega = valRPM*2*pi()/60;
A = pi()*(geomDIAMETER*0.5)^2;

F = flowRHO*A*(Omega*geomDIAMETER*0.5)^2;
M = flowRHO*A*(Omega)^2*(geomDIAMETER*0.5)^3;


% interpolate rest of variables
valFx=interp1(vecRPM',vecFx_rho,valRPM,'linear','extrap')*flowRHO;
valFy=interp1(vecRPM',vecFy_rho,valRPM,'linear','extrap')*flowRHO;
valMx=interp1(vecRPM',vecMx_rho,valRPM,'linear','extrap')*flowRHO;
valMy=interp1(vecRPM',vecMy_rho,valRPM,'linear','extrap')*flowRHO;
valQ=interp1(vecRPM',vecQ_rho,valRPM,'linear','extrap')*flowRHO;
valCP=interp1(vecRPM',vecCP,valRPM,'pchip','extrap');
valMUinf=interp1(vecRPM',vecMUinf,valRPM,'linear','extrap');


end

