%Performance metrics (HRT and N-removal efficiency) in two parameter space

%Select points along BP curve

%Get grid dimensions
[len,wid] = size(gridNS);

%Decay rates (not used in this model)
bA = 0;
bN = 0;

%Nitrogen incorporation
iNa = 0.083;
iNn = 0.083;
iNx = 0.058;
kq = 0;
for k = 1:len
    for kk = 1:wid
        
        switch flag3
            
            case 1
                WAS = 0.005;
                r = 0.5;
                DO = gr_limits1(k);
                rAmx = gr_limits2(kk);

            case 2
                rAmx = 86;
                r = 0.5;
                DO = gr_limits1(k);
                WAS = gr_limits2(kk);
                
            case 3
                DO = 1.5;
                r = 0.5;
                rAmx = gr_limits1(k);
                WAS = gr_limits2(kk);
                
        end
        
        %AMX concentration dependent on rAMX
        XAmx = rAmx/12.5; %mg/L

        %Biomass retention
        Bret = 1-WAS;
         
        params = [0.297,0.337,0.017,2.4,0.5,0.6,0.4,0.03,0.005,0.18,0.08,0.17,iNa,iNn,iNx,bA,bN,DO,XAmx]; %Growth parameters and fixed state values (O2 and AMX)
        tspan = linspace(0,50000,1000); % t large enough to converge to quasi-steady state impulse point
        
        s1in = 20; % NH4 in
        
        % Determine initial conditions (X2 > 0 or X2 = 0)
        % Load fittedmodel curve - Data generated using Matlab curve
        % fitting toolbox against branch point curve data (branch_curves.m)
        load(['fit_model',num2str(flag3)])
        
        switch flag3
            case 1
                
                % Determine location of point
                p1 = DO;
                p2 = rAmx;
                
                % Determine location of point on curve (@X)
                poc1 = fittedmodel1(p1);
                
                % Determine region of point [P(1,1) or P(1,0)]
                if p2 >= poc1
                    x2_init = 0;
                else
                    x2_init = 1;
                end
                
            case 2
                
                % Determine location of point
                p1 = DO;
                p2 = WAS;
                
                % Determine location of point on curve (@X)
                poc2 = fittedmodel2(p1);
                
                % Determine region of point
                if p2 >= poc2
                    x2_init = 0;
                else
                    x2_init = 1;
                end
                
            case 3
                
                % Determine location of point
                p1 = rAmx;
                p2 = WAS;
                
                % Determine location of point on curve (@X)
                poc3 = fittedmodel3(p1);
                
                % Determine region of point
                if p2 >= poc3
                    x2_init = 0;
                else
                    x2_init = 1;
                end
                
             case 4
                
                % Determine location of point
                p1 = DO;
                p2 = r;
                
                % Determine location of point on curve (@X)
                
                poc4 = fittedmodel4(p1);
                
                % Determine region of point
                
                if p2 >= poc4
                    x2_init = 0;
                else
                    x2_init = 1;
                end
        end
        
        init = [20;0;1;0;1;x2_init];
        
        s1bar = 2;
        
        imp = @(s1) s1 > s1bar;
        xin = [s1in];
        
        % Impulse points for each state variable [S1, S2, S3, S4, X1, X2]
        g = @(x) [(x(1)*(1-r)+xin*r) x(2)*(1-r) x(3)*(1-r) x(4)*(1-r) x(5:6)*Bret];
        
        %% Run solver
        [t,x,t0,x0,ximp] = ImpulseA(@(t,x)pna_sys_perf(t,x,params),g,imp,init,tspan);

        % Hydraulic Retention Time
        grad_time = gradient(t0);
        Tcycle = grad_time(end-1);
        
        HRT = Tcycle*24*2; % HRT (hours)
        
        HRTlist(k,kk) = HRT;
        
        % Nitrogen removal efficiency
        Neff_out(k,kk) = ximp(end,4)*100/s1in; %N-removal eff. (%)
        
        % Time to steady-state
        SS_time(k,kk) = t0(end); % Time (days)
        
        fprintf('Calculation for point %1.0f of %1.0f\n',kq,len*wid)
        kq=kq+1;
    end
end

% Plotting contour maps for performance metrics
% Change the colour scheme using:
%  colormap(...), use help colormap for more information and names of
%  colormaps

mat = dir('Data_contour/*.mat');
for q = 1:length(mat)-1
    load(['Data_contour/',mat(q).name]);
    
figure(q)
set(gca,'xlim',[min(par1) max(par1)])
hold on
set(gca,'ylim',[min(par2) max(par2)])
contourf(gr_limits1,gr_limits2,Neff_out',100,'linecolor','none');
hold on
plot(par1,par2,'r','linewidth',3)
xlabel(xparam,'fontsize',14)
ylabel(yparam,'fontsize',14)
colormap(brewermap([],'Greys'))
hcb=colorbar;
title(hcb,'N-removal (%)','fontsize',9)

figure(q+3)
set(gca,'xlim',[min(par1) max(par1)])
hold on
set(gca,'ylim',[min(par2) max(par2)])
contourf(gr_limits1(2:25),gr_limits2(1:25),(HRTlist(2:25,1:25)'),100,'linecolor','none');
hold on
plot(par1,par2,'r','linewidth',3)
xlabel(xparam,'fontsize',14)
ylabel(yparam,'fontsize',14)
colormap(brewermap([],'Greys'))
hcb=colorbar;
title(hcb,'HRT (d)','fontsize',9)
end

