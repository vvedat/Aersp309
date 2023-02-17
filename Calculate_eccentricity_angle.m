% Vedat Veziroglu
% I have completed this work with integrity.

a = 20540; %km
e = 0.7; %eccentricity
mu = 3.986e5;  %km^3/s^2
t = 12480; %secondd
tp = 1342; %second in periapsis
delta = 1e-4; %rad

E=kepler(a,e,mu,t,tp,delta); % calls Kepler function to find E through NR iteration

A = {'a', 'e','mu','t','tp','Delta','E'};  % A nicer notation
B = [a e mu t tp delta E];
C = {'km', '','km^3/s^2','s','s','rad','rad'};  % A nicer notation

% for k = 1:length(A)
%   fprintf('%s  %g\n', A{k}, B(k));
% end

C = cat(1, A, num2cell(B),C);
fprintf('%s = %g %s\n', C{:});



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

