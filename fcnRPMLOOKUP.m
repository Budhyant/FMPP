function [ rotorRPM, rotorPx, rotorPy, rotorMx, rotorMy, rotorQ,rotorCP, ...
            rotorMUinf ] = fcnRPMLOOKUP( flowq, flowRHO, ...
            valPITCHdeg, valTHRUST, tabLOOKUP, vecANGLELST, geomDIAMETER )

%% INTERPOLATE LOOKUP TABLE BASED ON PITCH, VELOCITY, AND THRUST
% Output - rotor forces and moments
        
valTHRUSTrho = valTHRUST/flowRHO;        
        
%v2 - includes Q lookup

[ highang, lowang ] = fcnFINDANG( valPITCHdeg, vecANGLELST );

% get RPM values
vecRPMhigh = unique(tabLOOKUP.RPM(tabLOOKUP.alpha == highang));
vecRPMlow = unique(tabLOOKUP.RPM(tabLOOKUP.alpha == lowang));

[ valRPM, valPx, valPy, valMx, valMy,...
    valQ, valCP, valMUinf ] ...
    = fcnTABINTERP( vecRPMhigh, highang, tabLOOKUP, flowq, ...
    valTHRUSTrho,flowRHO, geomDIAMETER );

[ valRPM(2), valPx(2), valPy(2), valMx(2), valMy(2),...
    valQ(2), valCP(2), valMUinf(2) ] ...
    = fcnTABINTERP( vecRPMlow, lowang, tabLOOKUP, flowq, ...
    valTHRUSTrho,flowRHO, geomDIAMETER  );


tempANGS = [ highang, lowang ];
rotorRPM = interp1(tempANGS,valRPM,valPITCHdeg);

rotorPx = -1*interp1(tempANGS,valPx,valPITCHdeg); % Px is negative in lookup tables but needs to be positive for convention
rotorPy = 1*interp1(tempANGS,valPy,valPITCHdeg);
rotorMx = -1*interp1(tempANGS,valMx,valPITCHdeg);
rotorMy = -1*interp1(tempANGS,valMy,valPITCHdeg);
rotorQ = -1*interp1(tempANGS,valQ,valPITCHdeg);
rotorCP = interp1(tempANGS,valCP,valPITCHdeg);
rotorMUinf = interp1(tempANGS,valMUinf,valPITCHdeg);

end

