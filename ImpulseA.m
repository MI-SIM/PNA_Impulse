function [t,x,ts,xs,ximp]=ImpulseA(f,g,imp,x0,ti)

% Calculation of Impulsive Map
% Matthew J. Wade (Newcastle University), Tyler Meadows (University of
% Idaho) and Ting-Hao Tsu (University of Miami)

%Last Updated: 05-03-2020

xtmp =@(x) imp(x(1))*x; % Rturns 0 if imp(x) if false
dx =@(t,x) imp(x(1))*f(t,xtmp(x));

t0 = ti(1);
Tend = ti(end);
Ttol = .001;
tolr = 1e-6;
tola = 1e-6;

Ntol = 0.0002; %Norm tolerance
NrmXs = 1; %Initial Norm

options = rdpset('RelTol',tolr,'AbsTol',tola,'MaxNbrStep',100000);
x0 = reshape(x0,1,[]);
ts = t0;
xs = x0;
ximp = x0;
t = t0;
x = x0;
%while (t0<Tend-Ttol) && NrmXs > Ntol 
while (t0<Tend-Ttol)
    [tnew,xnew]=RADAUsolver(dx,[t0,Tend],x0,options);
    x1new = xnew(:,1);
    ind1 = find([~imp(x1new);1],1)-1;
    t = [t; NaN; tnew(1:ind1)];
    x = [x; NaN*x0; xnew(1:ind1,:)];
    ts = [ts;t(end)];
    xs = [xs;g(x(end,:))];
    ximp = [ximp;x(end,:)];
    t0 = t(end);
    x0 = g(x(end,:));
    
    %Convergence check
    if size(xs,1)>2
        NrmXs = norm(xs(end-1,:)-xs(end-2,:));
    end
end

end