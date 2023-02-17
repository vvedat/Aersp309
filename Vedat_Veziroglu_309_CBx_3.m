%%% Vedat Veziroglu
%%% I have completed this work with integrity
clc;clear; close all

%% create rotation angles beta, alpha, theta
beta = 1.7; %rad
alpha = 4.2; %rad
theta = 2.5; %rad


%% Call rotation DCM function and translate A vector into B form
C_AB = DCM_3(theta)*DCM_1(alpha)*DCM_3(beta);

%% print out values
fprintf('Rotation angle around axis-3 beta is %g [rad]',beta)
fprintf('\nRotation angle around axis-1 alpha is %g [rad]',alpha)
fprintf('\nRotation angle around axis-3 theta is %g [rad]',theta)
fprintf('\nRotation DCM with 3-1-3 sequence by beta, alpha, and theta:\n')
disp(C_AB)

%% returns a rotation DCM around axis-3 by alpha degrees
function DCM_3 = DCM_3(alpha)
    DCM_3 =[cos(alpha), sin(alpha), 0
         -sin(alpha), cos(alpha), 0
         0          ,          0, 1];
end

%% returns a rotation DCM around axis-1 by alpha degrees
function DCM_1 = DCM_1(alpha)
    DCM_1 =[1, 0, 0
         0, cos(alpha), sin(alpha)
         0 , -sin(alpha), cos(alpha)];
end