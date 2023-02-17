%%% Vedat Veziroglu
%%% I have completed this work with integrity
clc ; clear; close all

%% position and velocity in b-i vectors
rB = [7145 ; 1564 ; 2832]; % (km)
vB = [1.4 ; -3.7 ; 2.9]; % (km/s)

%% define rotation angles alpha and lambda
alpha = 124; %deg
lambda= 41 ; %deg

%% Calculate rotation DCM 3-alpha and 2-lambda
C_NB = DCM_2(-lambda*pi/180)*DCM_3(alpha*pi/180); % convert degrees to radians for calling DCM function
C_BN = C_NB'; % tranpose C_NB to get rotation matrix from B to N

%% convert r and v into N basis vectors using C_BN
rN = C_BN*rB;
vN = C_BN*vB;

%% print values
fprintf('alpha = %g deg, lambda =  %g deg\n',[alpha, lambda])
fprintf('rB = [%g  %g  %g] (km)\n', rB)
fprintf('vB = [%g  %g  %g] (km/s)\n', vB)
disp("C_NB = ")
fprintf('[%g  %g  %g \n', C_NB(1,:))
fprintf(' %g  %g  %g \n', C_NB(2,:))
fprintf(' %g  %g  %g] \n', C_NB(3,:))
fprintf('rN = [%g  %g  %g] (km)\n', rN)
fprintf('vN = [%g  %g  %g] (km/s)\n', vN)

%% returns a rotation DCM around axis-3 by alpha rad
function DCM_3 = DCM_3(alpha)
    DCM_3 =[cos(alpha) ,    sin(alpha)  , 0
            -sin(alpha),    cos(alpha)  , 0
            0          ,    0           , 1];
end

%% returns a rotation DCM around axis-2 by alpha rad
function DCM_2 = DCM_2(alpha)
    DCM_2 =[cos(alpha)  , 0 ,   -sin(alpha)
            0           , 1 ,   0     
            sin(alpha)  , 0 ,   cos(alpha)];
end