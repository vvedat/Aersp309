%%% Vedat Veziroglu
%%% I have completed this work with integrity

%% create rotation angle variable
alpha = 0.6;

%% Call rotation DCM function
fprintf('\nRotation DCM around axis-3 by %g radians:\n\n', alpha)
disp(DCM_3(0.6))

%% returns a rotation DCM around axis-3 by alpha degrees
function DCM = DCM_3(alpha)
    DCM =[cos(alpha), cos(alpha-pi/2), 0
        cos(alpha+pi/2), cos(alpha), 0
        0, 0, 1];
end