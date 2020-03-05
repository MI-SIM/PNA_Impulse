
%% Impulsive system: MBBR with SBR cycles - Partial Nitritation/Anammox
% Matlab code using MatContM 5p4, command-line
% Developed by Matthew. J. Wade, 2019

% Original 0D hybrid model by Laureni et al. (2019). Water Research.
% https://www.sciencedirect.com/science/article/pii/S0043135419300016
% MatContM implementation by Ting-Hao Hsu (Uni. Miami)

% Last Updated: 05-03-2020

%% Preparation
try
    run('init')
catch
    run('./matcontm5p4/init'); 
end
global opt cds fpmds
opt = contset();

%% Set Numerics
sys = @pna_sys; % The model system
ap = 16;        % Index of active parameter for bifurcation

% Guide (see 'par' below): 16 = Dissolved Oxygen, oIN, 17 = AMX activity,
% R3, 20 = solids wastage, fwas

% These were the operating parameters used for this study

%% Set Parameters (See Table A2 in mansucript for units)
m1 = 0.3;       % Maximum specific growth rate, AOB
y1 = 0.18;      % Growth yield, AOB
ka1 = 2.4;      % Ammonium half-saturation constant, AOB
ko1 = 0.6;      % Oxygen half-saturation constant, AOB
in1 = 0.083;    % Nitrogen required for cell synthesis, AOB
m2 = 0.34;      % Maximum specific growth rate, NOB
y2 = 0.08;      % Growth yield, NOB
kn2 = 0.5;      % Nitrite half-saturation constant, NOB
ko2 = 0.4;      % Oxygen half-saturation constant, NOB
in2 = 0.083;    % Nitrogen required for cell synthesis, NOB
y3 = 0.17;      % Growth yield, AMX
ka3 = 0.03;     % Ammonium half-saturation constant, AMX
kn3 = 0.005;    % Nitrite half-saturation constant, AMX
in3 = 0.058;    % Nitrogen required for cell synthesis, AMX

aIN = 20;       % Ammonium influent concentration
oIN = 0.3;      % Default oxygen influent concentration
s1bar = 2;      % Threshold ammonium concentration to trigger SBR impulse
fwas = 0.005;   % Default fWAS (percentage / 100)
R3 = 86;        % Default maximum volumetric Anammox activity

r = .5;         % Volumetric SBR exchange fraction per cycle

par = [m1,ka1,ko1,y1,in1,m2,kn2,ko2,y2,in2,ka3,kn3,y3,in3,aIN,oIN,R3,s1bar,r,fwas].';

%% Set Functions
sp = r*aIN + (1-r)*s1bar;

% Initial state variables, obtained by running XPP (Impulse.ode)
xini = [0.2717899; 215.5964; 40.76058]; % [S2, X1, X2]

i = 1; % i iteration#

% Converge to initial fixed-point
[u0,v0] = init_FPm_FPm(sys,xini,par,ap,i);

%% Follow Periodic Orbit

% Options for continuation from fixed-point (impulse point)
opt = contset;
opt = contset(opt,'MaxNumPoints',200);
opt = contset(opt,'Singularities',1);
opt = contset(opt,'Multipliers',1);
opt = contset(opt,'backward',0);
opt = contset(opt,'Increment',1e-2);
opt = contset(opt, 'MaxNewtonIters', 3);
opt = contset(opt, 'MaxCorrIters', 10);
opt = contset(opt, 'MaxTestIters', 10);
opt = contset(opt, 'FunTolerance', 1e-6);
opt = contset(opt, 'VarTolerance', 1e-6);
opt = contset(opt, 'TestTolerance', 1e-6);
opt = contset(opt, 'InitStepSize', 0.01);

[u1,v1,s1,~,~] = cont(@fixedpointmap,u0,v0,opt);
opt = contset(opt,'backward',1);

u0n = u0;
v0n = v0;
niter = 1;
s2 = 0;
uz=[];

% Continue from last point until branch point is found
% niter should be increased if this loop fails
while length(s2)<=2 && niter<20
    
    [u2,v2,s2,~,~] = cont(@fixedpointmap,u0n,v0n,opt);
    uz =[uz,u2];
    [valc,dirc] = min([u2(4,1),u2(4,end)]);
    if dirc == 2
        indx = length(u2(4,:));
    else
        indx = 1;
    end
    u0n = u2(:,indx);
    v0n = v2(:,indx);
    niter = niter+1;
    
end

%Print eigenvalues
for k = 2:length(s1)-1
    fprintf('\n')
    fprintf(['Eigenvalues for ',s1(k).label,'\n'])
    fprintf('%f \n',s1(k).data.eval')
end

for k = 2:length(s2)-1
    fprintf('\n')
    fprintf(['Eigenvalues for ',s2(k).label,'\n'])
    fprintf('%f \n',s2(k).data.eval')
end

% Now continue to follow the branch curve
% Find branch point (BP) index on previous bifurcation curve

for ii = 1:length(s2)
    if strmatch(s2(ii).label,'BP  ');
        bpindex = ii;
    end
end

ap2 = [16 17]; %Bifurcate in two parameters [DO, Ramx], indices

u2i = u2(1:3,s2(bpindex).index); %Starting point
p1 = par; p1(fpmds.ActiveParams) = u2(4,s2(bpindex).index);

[ubp,vbp] = init_LPm_LPm(sys,u2i,p1,ap2,2); % Using continuation from limit
% point (LP) as conditions are similar to branch point

% Adjust continuation options if continuation stalls or fails to complete
% curve
opt = contset(opt,'backward',0);
opt = contset(opt,'singularities',0);
opt = contset(opt,'Increment',1e-2);
opt = contset(opt, 'MaxNewtonIters', 3);
opt = contset(opt, 'MaxCorrIters', 10);
opt = contset(opt, 'MaxTestIters', 10);
opt = contset(opt, 'FunTolerance', 1e-6);
opt = contset(opt, 'VarTolerance', 1e-6);
opt = contset(opt, 'TestTolerance', 1e-5);
opt = contset(opt, 'InitStepSize', 0.01);
opt = contset(opt, 'Adapt', 3);
opt = contset(opt,'MaxNumPoints',5);

[xbp1,vbp1,sbp1,hbp1,fbp1]=cont(@limitpointmap,ubp,vbp,opt);

opt = contset(opt,'backward',1);
[xbp2,vbp2,sbp2,hbp2,fbp2]=cont(@limitpointmap,ubp,vbp,opt);

% Plotting branch point curve in two parameter space
% Move figure to new directory to avoid overwriting if you run this code
% again and wish to retain these figures

figure(1)
cpl(xbp1,vbp1,sbp1,[4 5])
hold on
cpl(xbp2,vbp2,sbp2,[4 5])
xlabel('O_2')
ylabel('r_{AMX}')
savefig('Par_o2_rAMX.fig')

% Change indices and redo analysis
ap3 = [16 20]; %Two parameter [DO, fwas]

[ubpa,vbpa] = init_LPm_LPm(sys,u2i,p1,ap3,2);

opt = contset(opt,'backward',0);
[xbp3,vbp3,sbp3,hbp3,fbp3]=cont(@limitpointmap,ubpa,vbpa,opt);

opt = contset(opt,'backward',1);
[xbp4,vbp4,sbp4,hbp4,fbp4]=cont(@limitpointmap,ubpa,vbpa,opt);

figure(2)
cpl(xbp3,vbp3,sbp3,[4 5])
hold on
cpl(xbp4,vbp4,sbp4,[4 5])
xlabel('O_2')
ylabel('f_{WAS}')
savefig('Par_o2_fwas.fig')

% Change indices and redo analysis
ap4 = [17 20]; %Two parameter [RAMX, fwas]

[ubpb,vbpb] = init_LPm_LPm(sys,u2i,p1,ap4,2);

opt = contset(opt,'backward',0);
[xbp5,vbp5,sbp5,hbp5,fbp5]=cont(@limitpointmap,ubpb,vbpb,opt);

opt = contset(opt,'backward',1);
[xbp6,vbp6,sbp6,hbp6,fbp6]=cont(@limitpointmap,ubpb,vbpb,opt);

figure(3)
cpl(xbp5,vbp5,sbp5,[4 5])
hold on
cpl(xbp6,vbp6,sbp6,[4 5])
xlabel('r_{AMX}')
ylabel('f_{WAS}')
savefig('Par_rAMX_fwas.fig')

% Save data from analysis
save(['Data_Run_',date])

return
