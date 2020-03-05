function [ outS ] = pna_sys_perf(t,init,params)
%DM_5: Functions describing growth of nitrifiers on multiple substrates

%   Complete model of hybrid MBBR partial nitritation and Anammox;
%   Original model developed by Laureni et al. (2019). 
%   M.J. Wade (Newcastle University, McMaster University)

% Last Updated: 05-03-2020

%Parameters
muA = params(1); %Max Specific Growth Rate (AOB) 1/d
muN = params(2); %Max Specific Growth Rate (NOB) 1/d
muX = params(3); %Max Specific Growth Rate (AMX) 1/d
K1 = params(4); %Half-saturation coefficient (AOB on NH4) g/m3
K2 = params(5); %Half-saturation coefficient (NOB on NO2) g/m3
K0a = params(6); %Half-saturation coefficient (AOB on O2) g/m3
K0n = params(7); %Half-saturation coefficient (NOB on O2) g/m3
K1x = params(8); %Half-saturation coefficient (AMX on NH4) g/m3
K2x = params(9); %Half-saturation coefficient (AMX on NO2) g/m3
YA = params(10); %Growth yield constant (AOB) gCOD/gN
YN = params(11); %Growth yield constant (NOB) gCOD/gN
YX = params(12); %Growth yield constant (AMX) gCOD/gN
iA = params(13); %Fraction of N in AOB gN/gCOD
iN = params(14); %Fraction of N in NOB gN/gCOD
iX = params(15); %Fraction of N in AMX gN/gCOD
bA = params(16); %Decay constant for AOB 1/d
bN = params(17); %Decay constant for NOB 1/d

s0 = params(18); %Fixed O2 concentration mg/L
x3 = params(19); %Fixed AMX biomass concentration mg/L

%Initial state conditions
s1 = init(1); s2 = init(2); s3 = init(3);
s4 = init(4); x1 = init(5); x2 = init(6);

%Monod Growth functions

mu1 = muA*(s1/(K1+s1))*(s0/(K0a+s0)); %AOB growth on ammonium and oxygen

mu2 = muN*(s2/(K2+s2))*(s0/(K0n+s0)); %NOB growth on nitrite and oxygen

mu3 = (s1/(K1x+s1))*(s2/(K2x+s2)); %AMX growth on ammonium and nitrite

%ODEs

ds1 = -(1/YA+iA)*mu1*x1 - iN*mu2*x2 - (1/YX + iX)*mu3*x3; %+ iA*bA*x1 + iN*bN*x2; %NH4+

ds2 = 1/YA*mu1*x1 - 1/YN*mu2*x2 - (1/YX+1/1.14)*mu3*x3; %NO2-

ds3 = 1/YN*mu2*x2 + (1/1.14)*mu3*x3; %NO3-

ds4 = 2/YX*mu3*x3; %N2

dx1 = mu1*x1 - bA*x1; %AOB

dx2 = mu2*x2 - bN*x2; %NOB

outS = [ds1;ds2;ds3;ds4;dx1;dx2];

end

