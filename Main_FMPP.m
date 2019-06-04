%% FAST MULTIROTOR PERFORMANCE PREDICTION METHOD
clear,clc
close all

%% File Input
strFILE = 'inputs/SkyRangerTMotor.txt';

%% Turn analysis types ON (1) or OFF (0)

    analysisMOMENTtrim              = 1; % Turn moment trim on or off
    analysisROTORinterference       = 1; % Turn mutual wake interference velocity on or off
    analysisBODYinterference        = 1; % Turn body interference velocity on or off
    analysisBODYforces              = 1; % Turn body induced drag and lift forces on or off

%% Add point mass for mass offset test
    massOFFSET = 0; %kg
    positionOFFSET = [0, 0, 0];

    %  massOFFSET = 0.1; %kg
    %  positionOFFSET = [-0.1651, 0, 0]; 

%% Input Setup
% Load variables - atmospheric, component geometries, component masses
    [seqV, flowTEMP, flowALT, flowRHO, flowMU, flowM, flowR, flowALPHAT, ...
        angCLIMBdeg, angSIDEdeg, numLEADROTOR, geomTypeROTOR, geomNumROTORS,...
        geomDIAMETER, geomNumBLADES, geomARMlength, geomARMradius, ...
        geomBODYheight, geomBODYradius, geomLEGlength, geomLEGradius, ...
        geomLEGcentreradius, geomLEGcentreheight, geomPAYLOADlength, ...
        geomPAYLOADradius, geomPAYLOADheight, geomMOTORheight, ...
        geomMOTORradius, geomHUBheight, geomCGheight, massMOTOR, massARM, ...
        massLEG, massPAYLOAD, massBODY, massVEHICLE] ...
        = fcnMVPREAD(strFILE);

% New total vehicle mass with offset added
    massVEHICLE = massVEHICLE + massOFFSET;


% Calculate Reynolds number to CDY database for cylinder and sphere
    [cylinderRE, cylinderCDY, sphereRE, sphereCDY] = fcnRECURVE();


% Calculate wetted area of each component
    [areaLEG, areaARM, areaBODY, areaPAYLOAD, areaMOTOR]...
        = fcnCOMPONENTAREA(geomLEGlength, geomLEGradius,...
        geomARMlength, geomARMradius, geomBODYradius,...
        geomPAYLOADradius, geomPAYLOADlength,...
        geomMOTORradius, geomMOTORheight);

% Create a table with all rotor performance from database
    [tabLOOKUP, vecANGLELST] = fcnLOADTABLES(geomTypeROTOR);

% Change centre of gravity based on added point mass
    [geomCGoffset] = ...
        fcnCGOFFSET(massBODY,massOFFSET,positionOFFSET,geomCGheight);

% Setup coordinates for components
    [positionROTOR, positionMOTOR, positionARM, positionLEG,...
        positionBODY, positionPAYLOAD] ...
        = fcnCOORDSETUP(numLEADROTOR, geomNumROTORS,...
        geomARMlength, geomBODYradius, geomMOTORradius,...
        geomLEGcentreradius, geomLEGcentreheight, ...
        geomPAYLOADheight, geomHUBheight, geomCGheight, geomCGoffset);

%% Preallocate
    countV              = numel(seqV);
    
    flowV               = zeros(countV,1);
    flowq               = zeros(countV,1);
    
    % Force Trim reserved for number of velocity values under investigation
    powerPARASITIC      = zeros(countV,1);
    dragVEHICLE         = zeros(countV,1);
    dragARM             = zeros(countV,1);
    dragLEG             = zeros(countV,1);
    dragBODY            = zeros(countV,1);
    dragMOTOR           = zeros(countV,1);
    dragPAYLOAD         = zeros(countV,1);
    
    % Rotor forces and inflow angles - First Run
    rotorTHRUST         = zeros(countV,1,geomNumROTORS);
    rotorAngINFLOW      = zeros(countV,1,geomNumROTORS);
    rotorVelINFLOW      = zeros(countV,1,geomNumROTORS);
    rotorRPM            = zeros(countV,1,geomNumROTORS);
    dragBODYinduced     = zeros(countV,1); 
    liftBODY            = zeros(countV,1);
    pitchVEHICLEdeg     = zeros(countV,1);

   % Wake values 
    vi_int              = zeros(countV,3,geomNumROTORS);
    vi_self             = zeros(countV,3,geomNumROTORS);
    skewRAD             = zeros(countV,geomNumROTORS);  
    wi                  = zeros(geomNumROTORS,3, geomNumROTORS);
    rotorQ              = zeros(countV,1,geomNumROTORS);
    rotorPx             = zeros(countV,1,geomNumROTORS);
    rotorPy             = zeros(countV,1,geomNumROTORS);
    rotorMx             = zeros(countV,1,geomNumROTORS);
    rotorMy             = zeros(countV,1,geomNumROTORS);
    rotorCP             = zeros(countV,1,geomNumROTORS); 
    rotorCMx            = zeros(countV,1,geomNumROTORS);
    rotorJinf           = zeros(countV,1,geomNumROTORS);

    vi_body             = zeros(countV,3,geomNumROTORS);
    
    powerROTOR          = zeros(countV,1,geomNumROTORS);
    powerTOTAL          = zeros(countV,1,geomNumROTORS); 
    powerVEHICLE        = zeros(countV,1);

    momentROTORTHRUST   = zeros(countV,3,geomNumROTORS);
    momentROTORPx       = zeros(countV,3,geomNumROTORS); 
    momentROTORPy       = zeros(countV,3,geomNumROTORS);
    momentROTORQ        = zeros(countV,3,geomNumROTORS);
    momentROTORMx       = zeros(countV,3,geomNumROTORS);
    momentROTORMy       = zeros(countV,3,geomNumROTORS);
    momentWEIGHTMOTOR   = zeros(countV,3,geomNumROTORS);
    momentWEIGHTARM     = zeros(countV,3,geomNumROTORS);
    momentDRAGMOTOR     = zeros(countV,3,geomNumROTORS);
    momentDRAGARM       = zeros(countV,3,geomNumROTORS);
    momentWEIGHTLEG     = zeros(countV,3,geomNumROTORS);
    momentDRAGLEG       = zeros(countV,3,geomNumROTORS);
    momentWEIGHTBODY    = zeros(countV,3);
    momentWEIGHTPAYLOAD = zeros(countV,3);
    momentDRAGBODY      = zeros(countV,3);
    momentDRAGPAYLOAD   = zeros(countV,3);
    momentLIFTBODY      = zeros(countV,3);
    momentDRAGBODYinduced = zeros(countV,3);
    momentTOTAL         = zeros(countV,3);
    momentWEIGHTOFFSET  = zeros(countV,3);
    
    
    % Rotor forces and inflow angles - First Run
    rotorTHRUST1     = zeros(countV,1,geomNumROTORS);
    rotorRPM1        = zeros(countV,1,geomNumROTORS);
    
%% START VELOCITY SEQUENCE
for i = 1:size(seqV,1)
    flowV(i,1) = seqV(i)
    flowq(i,1) = 0.5*flowRHO*flowV(i,1)^2;

%% PARASITIC DRAG AND POWER PREDICTION MODEL -
    % PARASITIC DRAG OF EACH COMPONENT
    % PARASITIC POWER OF VEHICLE
    [powerPARASITIC(i,1), dragVEHICLE(i,1), dragARM(i,1), dragLEG(i,1),...
        dragBODY(i,1), dragMOTOR(i,1), dragPAYLOAD(i,1)] ...
        = fcnDRAGPREDICT(geomNumROTORS, flowV(i,1), flowRHO, flowMU,...
        cylinderRE, cylinderCDY, sphereRE, sphereCDY, areaARM, ...
        areaLEG, areaMOTOR, areaPAYLOAD, areaBODY, geomARMradius, ...
        geomLEGradius, geomMOTORradius, geomPAYLOADradius, geomBODYradius);
    
%% FORCE TRIM MODEL -    ROTOR THRUST AT FORCE TRIM
    % INDUCED DRAG AND LIFT OF BODY (FINAL)
    % PITCH OF VEHICLE (FINAL)
    % FIRST ITERATION INFLOW ANGLE, INFLOW VELOCITY
    % AND RPM
    [rotorTHRUST(i,1,:), rotorAngINFLOW(i,1,:), rotorVelINFLOW(i,1,:),...
        rotorRPM(i,1,:), dragBODYinduced(i,1), liftBODY(i,1),...
        pitchVEHICLEdeg(i,1)] ...
        = fcnFORCETRIM( flowq(i,1), flowRHO, geomNumROTORS, ...
        geomBODYradius, dragVEHICLE(i,1), massVEHICLE, ...
        tabLOOKUP, vecANGLELST,geomDIAMETER, analysisBODYforces);
    
%% BODY INTERFERENCE MODEL -
    % INTERFERENCE VELOCITY ON ROTORS IN ROTOR PLANE (FINAL)
        [vi_body(i,:,:)] ...
            = fcnBODYINTERFERENCE(flowV(i,1), geomBODYradius, ...
            pitchVEHICLEdeg(i,1), positionROTOR, geomNumROTORS,...
            analysisBODYinterference);

    
    momentTOTAL = [0, 0, 0];
    error = 1;
    count = 1;
    
while error(end) > 0.0001
%% RPM PREDICTION MODEL - AT FORCE TRIM CONDITIONS
    % FIRST ITERATION INPUTS - at force trim conditions
    % SUBSEEQUENT ITERATION INPUTS - uses updated thrust and RPM 
    [vi_int(i,:,:), vi_self(i,:,:), skewRAD(i,:,:), wi2, ...
        rotorAngINFLOW(i,1,:), rotorVelINFLOW(i,1,:), rotorvecVR(i,:,:),...
        rotorRPM(i,1,:), rotorPx(i,1,:), rotorPy(i,1,:), rotorMx(i,1,:),...
        rotorMy(i,1,:), rotorQ(i,1,:), rotorCP(i,1,:), rotorMUinf(i,1,:)]...
            = fcnPREDICTRPM(flowq(i,1),flowRHO, geomNumROTORS,geomNumBLADES,...
            geomDIAMETER, positionROTOR, rotorTHRUST(i,1,:),rotorRPM(i,1,:), ...
            rotorAngINFLOW(i,1,:),rotorVelINFLOW(i,1,:), pitchVEHICLEdeg(i,1),...
            vi_body(i,:,:),tabLOOKUP, vecANGLELST, analysisROTORinterference);
        
%% Save rotor wake interference values for plots
wi_rotor1(i,2) = wi2(1,3,2);
wi_rotor1(i,3) = wi2(1,3,3);
wi_rotor1(i,4) = wi2(1,3,4);

wi_rotor2(i,1) = wi2(2,3,1);
wi_rotor2(i,3) = wi2(2,3,3);
wi_rotor2(i,4) = wi2(2,3,4);

wi_rotor3(i,1) = wi2(3,3,1);
wi_rotor3(i,2) = wi2(3,3,2);
wi_rotor3(i,4) = wi2(3,3,4);

wi_rotor4(i,1) = wi2(4,3,1);
wi_rotor4(i,2) = wi2(4,3,2);
wi_rotor4(i,3) = wi2(4,3,3);        
        
%% MOMENT PREDICTION MODULE
    %               
    [ momentROTORTHRUST(i,:,:), momentROTORPx(i,:,:),...
        momentROTORPy(i,:,:), momentROTORMx(i,:,:),...
        momentROTORMy(i,:,:), momentROTORQ(i,:,:), ...
        momentWEIGHTMOTOR(i,:,:), momentWEIGHTARM(i,:,:),...
        momentDRAGMOTOR(i,:,:), momentDRAGARM(i,:,:),...
        momentWEIGHTLEG(i,:,:), momentDRAGLEG(i,:,:), ...
        momentWEIGHTBODY(i,:), momentWEIGHTPAYLOAD(i,:),...
        momentDRAGBODY(i,:), momentDRAGPAYLOAD(i,:),...
        momentLIFTBODY(i,:), momentDRAGBODYinduced(i,:), ...
        momentWEIGHTOFFSET(i,:), momentTOTALnew(i,:) ] ...
            = fcnCALCMOMENTS(massMOTOR, massARM, massLEG, massPAYLOAD,...
            massBODY, massOFFSET, positionROTOR, positionMOTOR,...
            positionARM, positionLEG, positionBODY, positionPAYLOAD,...
            positionOFFSET, dragVEHICLE(i,1), dragARM(i,1),...
            dragLEG(i,1), dragBODY(i,1), dragMOTOR(i,1), ...
            dragPAYLOAD(i,1), dragBODYinduced(i,1), liftBODY(i,1),...
            rotorTHRUST(i,1,:), rotorPx(i,1,:), rotorPy(i,1,:),...
            rotorMx(i,1,:), rotorMy(i,1,:),rotorQ(i,1,:),...
            pitchVEHICLEdeg(i,1), geomNumROTORS);
    
%% MOMENT TRIM MODULE 
    if analysisMOMENTtrim == 1
    [rotorTHRUST1(i,1,:), rotorRPM1(i,1,:), rotorPx(i,1,:), rotorPy(i,1,:), ...
            rotorMx(i,1,:), rotorMy(i,1,:), rotorQ(i,1,:),rotorCP(i,1,:), rotorMUinf(i,1,:)] ...
            = fcnMOMENTTRIM(rotorTHRUST(i,1,:), momentTOTALnew(i,:,:), positionROTOR,...
            flowV(i,1), flowRHO, pitchVEHICLEdeg(i,1),geomNumROTORS,tabLOOKUP, vecANGLELST,...
            numLEADROTOR, geomDIAMETER);    
        
    elseif analysisMOMENTtrim == 0
        
        rotorTHRUST1(i,1,:) = rotorTHRUST(i,1,:);
        rotorRPM1(i,1,:) = rotorRPM(i,1,:);
    
    end

% track difference between last iteration and new total moment
    error = abs(momentTOTALnew(end,2) - momentTOTAL(end,2))
    temperror(count,i) = error; % save error value based on number of iterations 

% save udpated variables
    momentTOTAL(i,:) = momentTOTALnew(i,:);    
    rotorRPM(i,1,:) = rotorRPM1(i,1,:);
    rotorTHRUST(i,1,:) = rotorTHRUST1(i,1,:);

    count = count+1 % display number of times moment loop repeats


end

%% POWER PREDICTION MODULE 
      [powerROTOR(i,1,:), powerTOTAL(i,1,:), powerVEHICLE(i,1)] ...
            = fcnROTORPOWER(flowRHO, geomDIAMETER, geomNumROTORS,...
            rotorCP(i,1,:), rotorRPM(i,1,:), powerPARASITIC(i,1));

end


 