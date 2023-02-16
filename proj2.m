%%% Vedat Veziroglu
%%% I have completed this work with integrity
clc ; clear; close all

%% known quantitie of s/c
rE = [-2195 ; 8153 ; -5700];        % [km] position in ECI coordinates
vE = [-5.080 ; -5.889 ; 0.239];     % [km/s] velocity in ECI coordinates
t1=7180;            % time in seconds
t2=24680;           % time at t2 [s]
mu=398600;      %gravitional parameter [km^3/s^]

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
e=norm(e);                           % scalar of eccentricity vector

% find correct theta based on rP
if rP(2)>= 0
    theta=acosd((p/norm(rP))/e-1/e);
else
    theta=-acosd((p/norm(rP))/e-1/e);
end

E=2*atan(tand(theta/2)*sqrt((1-e)/(1+e)));  % find eccentric anomoly
M = E-e*sin(E);                             % mean anomoly
tp=t1-sqrt(a^3/mu)*M;                        %find time of periapsis

% find rP2 and vP2 from t2
% find E for t2 using newton iteration (kepler function)
E2 = kepler(a,e,mu,t2,tp,0.00001);

% convert E to theta
theta2 = 2*atand(tan(E2/2)*sqrt((1+e)/(1-e)));

% find rp2 and vp2 from theta2
rP2 = p/(1+e*cosd(theta2));
vP2 = sqrt(2*(Ep+mu/rP2));
RP2 = [rP2*cosd(theta2), rP2*sind(theta2),0];

% find phi from h vector
h_P = C_EP*h;
phi2 = acosd(h_P(3)/rP2/vP2);
ang = theta2+90-phi2;
VP2 = [vP2*cosd(ang), vP2*sind(ang),0];

rE2 = C_EP'*RP2';
vE2 = C_EP'*VP2';

% print out 6 classical orbital elements
A = {'a', 'e','i','OMEGA','w','tp'};        % A nicer notation
B = [a e i*180/pi om*180/pi w*180/pi tp];
C = {'km', '','deg','deg','deg','seconds'}; % A nicer notation


C = cat(1, A, num2cell(B),C);
fprintf('%s = %g %s\n', C{:});
disp('DCM C^EP=')

fprintf('[%2.3g %2.3g %2.3g]\n ',C_EP)

fprintf('\nr_2E = [%g %g %g] (km)\n ',rE2)
fprintf('v_2E = [%g %g %g] (km/s)',vE2)

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

%% NR iteration Function to solve keplers equation
function E = kepler(a,e,mu,t,tp,delta)
    error = 2*delta;  % initialize error
    M = (t-tp)*sqrt(mu/a^3); % find M for 
    xold = M ; % define first guess as M
    while error>=delta
        f = M+e*sin(xold)-xold; % f(x)
        df = e*cos(xold)-1; % f'(x)
        xnew = xold - f/df; % find new x
        error =abs(xnew-xold); % find error
        xold=xnew; % define new x as old to be used in the next iteration
    end
    E=xold; % define E to be returned
end