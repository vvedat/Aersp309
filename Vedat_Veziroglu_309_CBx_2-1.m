%%% Vedat Veziroglu
%%% I have completed this work with integrity
clc,clear
%% create rotation angle and vector in A
beta = 1.2; %rad
rA = [4.2;7.3;-10.8]; %km

%% Call rotation DCM function and translate A vector into B form
C=DCM_3(beta);
rB=C*rA;

%% print out values
fprintf('Rotation angle beta is %g [rad]',beta)
fprintf('\nRotation DCM around axis-3 by beta:\n')
disp(C)
fprintf('R vector in frame A is [%g,%g,%g](km):',rA)
fprintf('\nR vector in frame B is [%g,%g,%g](km):',rB)



%% returns a rotation DCM around axis-3 by alpha degrees
function DCM = DCM_3(alpha)
    DCM =[cos(alpha), sin(alpha), 0
         -sin(alpha), cos(alpha), 0
         0          ,          0, 1];
end