%% TOTALS - First

momentTOTALROTORFORCES = momentROTORTHRUST(:,:,1)+momentROTORTHRUST(:,:,2)+ ...
                    momentROTORTHRUST(:,:,3)+momentROTORTHRUST(:,:,4)+ ...
                    momentROTORPx(:,:,1)+momentROTORPx(:,:,2)+ ...
                    momentROTORPx(:,:,3)+momentROTORPx(:,:,4)+ ...
                    momentROTORPy(:,:,1)+momentROTORPy(:,:,2)+ ...
                    momentROTORPy(:,:,3)+momentROTORPy(:,:,4);
                    
momentTOTALROTORMOMENTS = momentROTORMx(:,:,1)+momentROTORMx(:,:,2)+ ...
                    momentROTORMx(:,:,3)+momentROTORMx(:,:,4)+ ...
                    momentROTORMy(:,:,1)+momentROTORMy(:,:,2)+ ...
                    momentROTORMy(:,:,3)+momentROTORMy(:,:,4);

momentTOTALROTORTORQUES = momentROTORQ(:,:,1)+momentROTORQ(:,:,2)+ ...
                    momentROTORQ(:,:,3)+momentROTORQ(:,:,4);
                
momentTOTALDRAG =  momentDRAGMOTOR(:,:,1)+momentDRAGMOTOR(:,:,2)+...
                   momentDRAGMOTOR(:,:,3)+momentDRAGMOTOR(:,:,4)+...
                   momentDRAGARM(:,:,1)+momentDRAGARM(:,:,2)+...
                   momentDRAGARM(:,:,3)+momentDRAGARM(:,:,4)+...
                   momentDRAGLEG(:,:,1)+momentDRAGLEG(:,:,2)+...
                   momentDRAGLEG(:,:,3)+momentDRAGLEG(:,:,4)+...
                   momentDRAGBODY(:,:)+momentDRAGPAYLOAD(:,:)+...
                   momentDRAGBODYinduced(:,:);

momentTOTALWEIGHTS = momentWEIGHTMOTOR(:,:,1)+momentWEIGHTMOTOR(:,:,2)+...
                   momentWEIGHTMOTOR(:,:,3)+momentWEIGHTMOTOR(:,:,4)+...
                   momentWEIGHTARM(:,:,1)+momentWEIGHTARM(:,:,2)+...
                   momentWEIGHTARM(:,:,3)+momentWEIGHTARM(:,:,4)+...
                   momentWEIGHTLEG(:,:,1)+momentWEIGHTLEG(:,:,2)+... 
                   momentWEIGHTLEG(:,:,3)+momentWEIGHTLEG(:,:,4)+...
                   momentWEIGHTBODY(:,:)+momentWEIGHTPAYLOAD(:,:)+momentWEIGHTOFFSET(:,:);

momentTOTALLIFTBODY = momentLIFTBODY(:,:);

%% PITCH
 figure;
 plot(flowV, pitchVEHICLEdeg)
title('pitch vs airspeed') 

%% THRUST
% Rotor thrust (force trim only) vs airspeed
figure;
plot(flowV, rotorTHRUST(:,:,1),'-ro', flowV, rotorTHRUST(:,:,2),'-b+',...
    flowV, rotorTHRUST(:,:,3),'-ks',flowV, rotorTHRUST(:,:,4),'-gx')
title('Rotor Thrust vs Airspeed Force Trim')
xlabel('Airspeed [m/s]')
ylabel('Rotor Thrust [N]')
legend('Rotor 1 - Lead','Rotor 2 - Left','Rotor 3 - Rear','Rotor 4 - Right','Location','northwest')

%% RPM
% Rotor RPM (after moment trim) vs airspeed
figure;
hold on
plot(flowV, rotorRPM(:,:,1),'-ro', flowV, rotorRPM(:,:,2),'-b+',...
    flowV, rotorRPM(:,:,3),'-ks',flowV, rotorRPM(:,:,4),'-gx')
title('Rotor RPM vs Airspeed')
xlabel('Airspeed [m/s]')
ylabel('Rotor RPM')
legend('Rotor 1 - Lead','Rotor 2 - Left','Rotor 3 - Rear','Rotor 4 - Right','Location','northwest')

%% POWER
% Rotor Power (force trim only) vs airspeed
figure;
plot(flowV, powerROTOR(:,:,1),'-ro', flowV, powerROTOR(:,:,2),'-b+',...
    flowV, powerROTOR(:,:,3),'-ks',flowV, powerROTOR(:,:,4),'-gx')
title('Rotor Power vs Airspeed')
xlabel('Airspeed [m/s]')
ylabel('Rotor Power [W]')
legend('Rotor 1 - Lead','Rotor 2 - Left','Rotor 3 - Rear','Rotor 4 - Right','Location','northwest')

%% POWER DECOMPOSITION
% Rotor Power (force trim only) vs airspeed
figure;
plot(flowV, sum(powerROTOR,3),'-.r', ...
    flowV, powerPARASITIC, '--g', flowV, sum(powerROTOR,3)+powerPARASITIC, '-k')
title('Rotor Power vs Airspeed')
xlabel('Airspeed [m/s]')
ylabel('Power [W]')
legend('total','parasitic','total rotor','Location','northwest')

figure(21);
hold on
plot(flowV, sum(powerROTOR,3)+powerPARASITIC, '-k')
xlabel('Airspeed [m/s]')
ylabel('Power [W]')


%% Residual moments (force and moment trim) vs airspeed
% About y-axis only (pitch)
% Rotor Forces - Thrust, Px, Py
% Rotor Moments - Mx, My
% Drag - Dpar, Dind
% Lift - from body
% Weights - from components
figure;
plot(flowV, momentTOTALROTORFORCES(:,1),'--b+',...
    flowV, momentTOTALROTORMOMENTS(:,1),'--mo',...
    flowV, momentTOTALROTORTORQUES(:,1),'--co',...
    flowV, momentTOTALDRAG(:,1),':rs',...
    flowV, momentTOTALLIFTBODY(:,1),'-.kh',...
    flowV, momentTOTALWEIGHTS(:,1), '-gd')
hold on
plot (flowV, momentTOTALnew(:,1), 'LineWidth',2)
title('Roll')
xlabel('Airspeed [m/s]')
ylabel('Residual Rolling Moment [Nm]')
legend('Rotor Forces','Rotor Moments','Rotor Torques','Total Drag','Lift','Weights','Total','Location','southwest')

figure;
plot(flowV, momentTOTALROTORFORCES(:,2),'--b+',...
    flowV, momentTOTALROTORMOMENTS(:,2),'--mo',...
    flowV, momentTOTALROTORTORQUES(:,2),'--co',...
    flowV, momentTOTALDRAG(:,2),':rs',...
    flowV, momentTOTALLIFTBODY(:,2),'-.kh',...
    flowV, momentTOTALWEIGHTS(:,2), '-gd')
hold on
plot (flowV, momentTOTALnew(:,2), 'LineWidth',2)
title('Pitch diamond Force Trim')
xlabel('Airspeed [m/s]')
ylabel('Pitching Moment [Nm]')
legend('Rotor Forces','Rotor Moments','Rotor Torques','Total Drag','Lift','Weights','Total','Location','southwest')


figure;
plot(flowV, momentTOTALROTORFORCES(:,3),'--b+',...
    flowV, momentTOTALROTORMOMENTS(:,3),'--mo',...
    flowV, momentTOTALROTORTORQUES(:,3),'--co',...
    flowV, momentTOTALDRAG(:,3),':rs',...
    flowV, momentTOTALLIFTBODY(:,3),'-.kh',...
    flowV, momentTOTALWEIGHTS(:,3), '-gd')
hold on
plot (flowV, momentTOTALnew(:,3), 'LineWidth',2)
title('Yaw')
xlabel('Airspeed [m/s]')
ylabel('Residual Yawing Moment [Nm]')
legend('Rotor Forces','Rotor Moments','Rotor Torques','Total Drag','Lift','Weights','Total','Location','southwest')

%% Rotor Interference
figure(61);
plot(flowV, wi_rotor1(:,2), '-b+',...
    flowV, wi_rotor1(:,3), '-ks',...
    flowV,wi_rotor1(:,4), '-gx')
xlabel('Airspeed [m/s]')
ylabel('Induced Velocity Through Rotor Plane [m/s]')
legend('Rotor 2','Rotor 3','Rotor 4')

figure(62);
plot(flowV, wi_rotor2(:,1), '-ro',...
    flowV, wi_rotor2(:,3), '-ks',...
    flowV,wi_rotor2(:,4), '-gx')
xlabel('Airspeed [m/s]')
ylabel('Induced Velocity Through Rotor Plane [m/s]')
legend('Rotor 1','Rotor 3','Rotor 4')

figure(63);
plot(flowV, wi_rotor3(:,1), '-ro',...
    flowV, wi_rotor3(:,2), '-b+',...
    flowV,wi_rotor3(:,4), '-gx')
xlabel('Airspeed [m/s]')
ylabel('Induced Velocity Through Rotor Plane [m/s]')
legend('Rotor 1','Rotor 2','Rotor 4')

figure(64);
plot(flowV, wi_rotor4(:,1), '-ro',...
    flowV, wi_rotor4(:,2), '-b+',...
    flowV,wi_rotor4(:,3), '-ks')
legend('Rotor 1','Rotor 2','Rotor 3') 

%% BODY INTERFERENCE
vecBODY = vi_body;

B=repmat(pitchVEHICLEdeg, 1,1,4);
Vbxx = vecBODY(:,1,:).*cosd(B(:,1,:));
Vbxz = vecBODY(:,1,:).*sind(B(:,1,:));
Vbzx = vecBODY(:,3,:).*sind(B(:,1,:));
Vbzz = vecBODY(:,3,:).*cosd(B(:,1,:));

Vbrx = Vbxx + Vbxz;
Vbrz = Vbzx + Vbzz;
Vbry = zeros(length(vecBODY),1,geomNumROTORS);
vecBODYrotor = [Vbrx, Vbry , Vbrz];

xINT = zeros(length(vecBODY),1,geomNumROTORS);
yINT = zeros(length(vecBODY),1,geomNumROTORS);
zINT = vi_int(:,3,:);
vecINTrotor = [xINT, yINT,zINT];

Vlist = repmat(seqV,1,1,4);
xFREE = Vlist.*cosd(B);
yFREE = zeros(length(vecBODY),1,4);
zFREE = -Vlist.*sind(B);
vecFREErotor = [xFREE, yFREE,zFREE]; 

vecVRrotor = vecFREErotor + vecINTrotor + vecBODYrotor;
magVRrotor = sqrt(sum(abs(vecVRrotor).^2,2));

%% Body Interference
figure;
plot(flowV, Vbrx(:,:,1),'-ro',...
    flowV, Vbrx(:,:,2),'-b+',...
    flowV, Vbrx(:,:,3),'-ks',...
    flowV, Vbrx(:,:,4),'-gx')
title('Normal Component of Induced Velocity')
xlabel('Airspeed [m/s]')
ylabel('Body Interference Velocity Parallel to Rotor Plane [m/s]')
legend('Rotor 1 - Lead Left','Rotor 2 - Rear Left','Rotor 3 - Rear Right','Rotor 4 - Lead Right')

figure;
plot(flowV, Vbrz(:,:,1),'-ro',...
    flowV, Vbrz(:,:,2),'-b+',...
    flowV, Vbrz(:,:,3),'-ks',...
    flowV, Vbrz(:,:,4),'-gx')
title('Normal Component of Induced Velocity')
xlabel('Airspeed [m/s]')
ylabel('Body Interference Velocity Through to Rotor Plane [m/s]')
legend('Rotor 1 - Lead Left','Rotor 2 - Rear Left','Rotor 3 - Rear Right','Rotor 4 - Lead Right')




% plot rotor 1 combined velocities in x
plot(flowV, vecBODYrotor(:,1,1),'.-g',...
    flowV, vecINTrotor(:,1,1),'.-r',...
    flowV, vecFREErotor(:,1,1),'.-k',...
    flowV, vecVRrotor(:,1,1),'b')
title('rotor 1 in x')
xlabel('Vehicle Airspeed [m/s]')
ylabel('Inflow Velocity [m/s]')
legend('Body Interference','Rotor Interference','Freestream','Total')
% plot rotor 1 combined velocities in z
figure;
plot(flowV, vecBODYrotor(:,3,1),'.-g',...
    flowV, vecINTrotor(:,3,1),'.-r',...
    flowV, vecFREErotor(:,3,1),'.-k',...
    flowV, vecVRrotor(:,3,1),'b')
title('rotor 1 in z')
xlabel('Vehicle Airspeed [m/s]')
ylabel('Inflow Velocity [m/s]')
legend('Body Interference','Rotor Interference','Freestream','Total')

% plot rotor 2 combined velocities in x
figure;
plot(flowV, vecBODYrotor(:,1,2),'.-g',...
    flowV, vecINTrotor(:,1,2),'.-r',...
    flowV, vecFREErotor(:,1,2),'.-k',...
    flowV, vecVRrotor(:,1,2),'b')
title('rotor 2 in x')
xlabel('Vehicle Airspeed [m/s]')
ylabel('Inflow Velocity [m/s]')
legend('Body Interference','Rotor Interference','Freestream','Total')
% plot rotor 2 combined velocities in z
figure;
plot(flowV, vecBODYrotor(:,3,2),'.-g',...
    flowV, vecINTrotor(:,3,2),'.-r',...
    flowV, vecFREErotor(:,3,2),'.-k',...
    flowV, vecVRrotor(:,3,2),'b')
title('rotor 2 in z')
xlabel('Vehicle Airspeed [m/s]')
ylabel('Inflow Velocity [m/s]')
legend('Body Interference','Rotor Interference','Freestream','Total')
% plot rotor 3 combined velocities in x
figure;
plot(flowV, vecBODYrotor(:,1,3),'.-g',...
    flowV, vecINTrotor(:,1,3),'.-r',...
    flowV, vecFREErotor(:,1,3),'.-k',...
    flowV, vecVRrotor(:,1,3),'b')
title('rotor 3 in x')
xlabel('Vehicle Airspeed [m/s]')
ylabel('Inflow Velocity [m/s]')
legend('Body Interference','Rotor Interference','Freestream','Total')
% plot rotor 3 combined velocities in z
plot(flowV, vecBODYrotor(:,3,3),'.-g',...
    flowV, vecINTrotor(:,3,3),'.-r',...
    flowV, vecFREErotor(:,3,3),'.-k',...
    flowV, vecVRrotor(:,3,3),'b')
title('rotor 3 in z')
xlabel('Vehicle Airspeed [m/s]')
ylabel('Inflow Velocity [m/s]')
legend('Body Interference','Rotor Interference','Freestream','Total')

% plot rotor 4 combined velocities in x
figure;
plot(flowV, vecBODYrotor(:,1,4),'.-g',...
    flowV, vecINTrotor(:,1,4),'.-r',...
    flowV, vecFREErotor(:,1,4),'.-k',...figure;
    flowV, vecVRrotor(:,1,4),'b')
title('rotor 4 in x')
xlabel('Vehicle Airspeed [m/s]')
ylabel('Inflow Velocity [m/s]')
legend('Body Interference','Rotor Interference','Freestream','Total')
% plot rotor 4 combined velocities in z
figure;
plot(flowV, vecBODYrotor(:,3,4),'.-g',...
    flowV, vecINTrotor(:,3,4),'.-r',...
    flowV, vecFREErotor(:,3,4),'.-k',...
    flowV, vecVRrotor(:,3,4),'b')
title('rotor 4 in z')
xlabel('Vehicle Airspeed [m/s]')
ylabel('Inflow Velocity [m/s]')
legend('Body Interference','Rotor Interference','Freestream','Total')


