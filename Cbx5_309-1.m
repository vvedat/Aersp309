%%% Vedat Veziroglu
%%% I have completed this work with integrity
clc ; clear; close all

%% known quantitie of s/c
rE = [7350 ; -2264 ; 2938];     % (km) position in ECI coordinates
vE = [-3.426 ; 5.994 ; -2.328]; % (km/s) velocity in ECI coordinates
t=1564;         % time in seconds
mu=398600;      %gravitional parameter km^3/s^2

%% define I,J,K unit vectors for ECI
I = [1 0 0];
J = [0 1 0];
K = [0 0 1];

%% find omega, inclination, and OMEGA
h = cross(rE,vE);               % angular velocity
i=acos(K*h/norm(h));            % inclination
n = cross(K,h)/norm(cross(K,h)); % node unit vector
e=cross(vE,h)/mu - rE./norm(rE); % eccentricity vector

% calculate omega based on K and e and n
if K*e >= 0
    w=acos(n*e/norm(e));
else
    w=-acos(n*e/norm(e));
end

% calculate OMEGA based on n and I and J
if n*I' >= 0
    om = atan((n*J')/(n*I'));
else
    om = atan((n*J')/(n*I'))+pi;
end


%% Calculate rotation 3-1-3 euler sequence 3-OMEGA -> 1-inclination -> 3-omega
C_EP = DCM_3(w)*DCM_1(i)*DCM_3(om); 

%% convert r and v into perifocal coordinates vectors using C_EP
rP = C_EP*rE;
vP = C_EP*vE;

%% calculate classical orbit elements
p=(norm(h))^2/mu;                    % p
Ep = (norm(vP))^2/2 - mu/norm(rP);   % epsilon for energy
a=-mu/Ep/2;                          % semi-major axis
e=norm(e);                           %scalar of eccentricity vector

% find correct theta based on rP
if rP(2)>= 0
    theta=acosd((p/norm(rP))/e-1/e);
else
    theta=-acosd((p/norm(rP))/e-1/e);
end

E=2*atan(tand(theta/2)*sqrt((1-e)/(1+e)));  % find eccentric anomoly
M = E-e*sin(E);                             % mean anomoly
tp=t-sqrt(a^3/mu)*M;                        %find time of periapsis

% print out 6 classical orbital elements
A = {'a', 'e','i','OMEGA','w','tp'};        % A nicer notation
B = [a e i*180/pi om*180/pi w*180/pi tp];
C = {'km', '','deg','deg','deg','seconds'}; % A nicer notation

C = cat(1, A, num2cell(B),C);
fprintf('%s = %g %s\n', C{:});

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

%% returns a rotation DCM around axis-1 by alpha degrees
function DCM_1 = DCM_1(alpha)
    DCM_1 =[1, 0            , 0
            0, cos(alpha)   , sin(alpha)
            0, -sin(alpha)  , cos(alpha)];
end