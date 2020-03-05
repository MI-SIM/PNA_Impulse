function out = pna_sys

% Mathematical model of PN/A hybrid system
% Code developed by Matthew J. Wade, 2019

% Last Updated: 05-03-2020

out{1} = [];%@init;
out{2} = @fun_eval;
out{3} = [];
out{4} = [];
out{5} = [];
out{6} = [];
out{7} = [];
out{8} = [];
out{9} = [];

% --------------------------------------------------------------------------
function unew = fun_eval(~,u,m1,ka1,ko1,y1,in1,m2,kn2,ko2,y2,in2,ka3,kn3,y3,in3,aIN,oIN,R3,s1bar,r,fwas)

rmin = .01;
rmax = .99;
r = min(r,rmax);
r = max(r,rmin);

% Calculate starting point of impulse
sp = r*aIN + (1-r)*s1bar;

s2n = u(1);
x1n = u(2);
x2n = u(3);
so2 = oIN;

%ODE Kinetics
syms s1 s2 x1 x2

% Calculation of Anammox process rate (\rho in literature)
p3 = R3/12.5;

% Growth functions (Monod form)
mua = m1*(s1./(ka1+s1));
mun = m2*(s2./(kn2+s2));
muo1 = (so2./(ko1+so2));
muo2 = (so2./(ko2+so2));
muaa = p3*(s1./(ka3+s1));
muna = (s2./(kn3+s2));

F1 = mua*muo1;
F2 = mun*muo2;
F3 = muaa*muna;

% Differential equations
f1 = -(1/y1 + in1)*F1*x1 - in2*F2*x2 - (1/y3 + in3)*F3;  %dS1/dt
f2 = (1/y1)*F1*x1 - (1/y2)*F2*x2 - ((1/y3)+(1/1.14))*F3; %dS2/dt
f3 = F1*x1;                                              %dX1/dt
f4 = F2*x2;                                              %dX2/dt

deq = matlabFunction([f2./f1;f3./f1;f4./f1]);
deq1 = @(s,w) deq(s,w(1),w(2),w(3));

de_opt = odeset('RelTol',1e-6,'AbsTol',1e-6);
sspan = [sp s1bar];

w0 = [(1-r)*s2n;(1-fwas)*x1n;(1-fwas)*x2n];

% Integrate system of ODEs between impulse times
[~,w] = ode23s(deq1,sspan,w0,de_opt);

s2m = w(end,1); % Nitrite concentration
x1m = w(end,2); % AOB concentration
x2m = w(end,3); % NOB concentration

unew = [s2m;x1m;x2m];


