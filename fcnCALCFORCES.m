function [ valPITCHdeg_in, valPITCHdeg_out, rotorTHRUST, liftBODY, dragBODYinduced, ...
            rotorRPM, rotorFx, rotorFy, rotorMx, rotorMy, rotorQ, rotorCP, ...
            rotorMUinf] = fcnCALCFORCES( flowq, flowRHO, geomBODYradius, ...
            numrotors, tempPITCHrad, dragVEHICLE, valWEIGHT, tabLOOKUP, vecANGLELST,...
            geomDIAMETER, analysisBODYforces)       
%% CALCULATE FORCES IN FORCE TRIM
% Input:
%       - flowq - 1x1 - dynamic pressure
%       - flowRHO - 1x1 - density at altitude
%       - geomBODYradius - 1x1 - body radius
%       - numrotors - 1x1 - number of rotors
%       - tempPITCHrad - 1x1 - temporary pitch value in radians (last
%       iteration pitch result)
% Output:
%       - valPITCHdeg_in - 1x1 - last iteration pitch result
%       - valPITCHdeg_out - 1x1 - new pitch result
%       - rotorTHRUST - 1x1 - thrust of one rotor
%       - liftBODY - 1x1 - lift of "puck" shaped central body
%       - dragBODYinduced - 1x1 - induced drag of "puck" shaped central body
%       - rotorRPM, rotorFx, rotorFy, rotorMx, rotorMy, rotorQ, rotorCP,
%       rotorMUinf - all 1x1 - rotor variables will be updated in moment
%       trim module


    valPITCHdeg_in      = rad2deg(tempPITCHrad);
    
    if analysisBODYforces == 1
        % Calculate lift and induced drag of body
        coefLIFT            = 1.8*(tempPITCHrad); % Assume body has puck shape
        coefDRAGinduced     = 0.81*(tempPITCHrad)^2;

        liftBODY            = flowq*pi()*geomBODYradius^2*coefLIFT;
        dragBODYinduced     = flowq*pi()*geomBODYradius^2*coefDRAGinduced;
    
    elseif analysisBODYforces == 0
        
        liftBODY            = 0;
        dragBODYinduced     = 0;
                
    end
    
    % Guess thrust prediction for lookup table
    tempTHRUST          = sqrt((valWEIGHT+liftBODY)^2 + ...
                                (dragVEHICLE+dragBODYinduced)^2)/numrotors; % first guess thrust
    
    % Lookup Px for Eqns. valPITCH and rotorTHRUST   
    [ rotorRPM, rotorFx, rotorFy, rotorMx, rotorMy,rotorQ, rotorCP, rotorMUinf ] ...
        = fcnRPMLOOKUP( flowq, flowRHO, valPITCHdeg_in, tempTHRUST, ...
        tabLOOKUP, vecANGLELST, geomDIAMETER );
    
    
    valPITCHdeg_out     = asind((dragVEHICLE/numrotors*tempTHRUST+rotorFx*valWEIGHT/numrotors+...
                                rotorFx*liftBODY/numrotors+tempTHRUST*dragBODYinduced/numrotors)...
                                /(tempTHRUST^2+rotorFx^2));
    
                             
    rotorTHRUST = sqrt((valWEIGHT/numrotors+liftBODY/numrotors...
                        +rotorFx*sind(valPITCHdeg_out))^2 + ...
                        (dragVEHICLE/numrotors+dragBODYinduced/numrotors...
                        +rotorFx*cosd(valPITCHdeg_out))^2);
    


end

