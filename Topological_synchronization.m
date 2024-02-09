%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Script that numerically integrates the dynamics of Dirac synchronization
%on fully connected networks as described in the paper
%
%
% This code can be redistributed and/or modified
% under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or (at
% your option) any later version.
%  
% This code is distributed by the authors in the hope that it will be 
% useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
%
% If you use this code please cite the following paper:
%
%”Dirac synchronization is rhythmic and explosive”
%L. Calmon, J. G. Restrepo, J. Torres, G. Bianconi, Communications Physics 5.1 (2022): 1-17.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialisation of the system
clear all

%%%%%%Parameters for numerical integration in RK4%%%%%%%%%
Tmax = 10; % max time to integrate
dt = 0.005; %time step
sigma_max=9; %max coupling strength
dsigma=0.03;
func = @ddt_fc; %dynamical equations to be loaded from ddt_fc.m

%%%%%%%%Frequencies and phases initialisation%%%%%%%%

N=100; %number of nodes
Omega0 = 0; %nodes mean frequencies
Omega1=0; %edges mean frequencies

w = Omega0 + randn(N,1); %Gaussian node frequencies
w = w - mean(w); %enforces mean 0

omegaedge = Omega1/N+ 1*randn(N,N)/(sqrt(N));
omegaedge=triu(omegaedge,1);
omegaedge=omegaedge-sum(sum(omegaedge))/nnz(omegaedge);
omegaedge=triu(omegaedge,1); %Appropriately scaled Gaussian links frequencies

what=sum(omegaedge,2)-(sum(omegaedge,1))'; %Normally distributed projected edges frequencies

theta0 = 2*pi*rand(N,1); % random initial conditions on nodes
psi0 = 2*pi*rand(N,1); % random initial conditions for link phases projected to nodes

%initialise the phases
theta = theta0;
psi = psi0;

%%%%%%%%%%%%%variables that will store the time series%%%%%%%%%%
Xats_up = cell(length(0:dsigma:sigma_max),1); %time series of X_\alpha in forward direction
Xbts_up = cell(length(0:dsigma:sigma_max),1); %time series of X_\beta in forward direction

Xthetats_up = cell(length(0:dsigma:sigma_max),1); %time series of X_\theta in forward direction
Xpsits_up = cell(length(0:dsigma:sigma_max),1); %time series of X_\beta in forward direction

Xats_down = cell(length(sigma_max:-dsigma:0),1); %same in backward direction.
Xbts_down = cell(length(sigma_max:-dsigma:0),1);

Xthetats_down = cell(length(0:dsigma:sigma_max),1);
Xpsits_down = cell(length(0:dsigma:sigma_max),1);

%% Numerical integration and construction of the phase diagram: first in forward direction
kcount = 1; %coupling strength step counter

for sigma = 0:dsigma:sigma_max %sweep sigma up
    
    Normcount = 0; %time series averages counter
    
    %variables to store the averaged real order parameters
    Rbetaup_ave = 0;
    Ralphaup_ave = 0;
    Rthetaup_ave = 0;
    Rpsiup_ave = 0;
    
    dangle=0; %phase of Xalpha
    D_angle_ave=0; %average of emergent frequency
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    te = 1; %time step counter
    
    for t = 0:dt:Tmax %time iterations for numerical integration
        
        %Runge-Kutta 4
        [k1{1},k1{2}] = func(theta, psi, sigma,w,what);
        [k2{1},k2{2}] = func(theta+0.5*dt*k1{1}, psi+0.5*dt*k1{2},sigma,w,what);
        [k3{1},k3{2}] = func(theta+0.5*dt*k2{1}, psi+0.5*dt*k2{2},sigma,w,what);
        [k4{1},k4{2}] = func(theta+dt*k3{1}, psi+dt*k3{2},sigma,w,what);
        
        %Final update to the phases
        theta = theta + (dt/6)*(k1{1}+2*k2{1}+2*k3{1}+k4{1});
        psi = psi + (dt/6)*(k1{2}+2*k2{2}+2*k3{2}+k4{2});
        
        %alpha and beta from theta and psi at the time step
        alpha = theta + psi/2;
        beta = (N/(N-1))*(theta-sum(theta)/N) - psi;
        
        %Complex order parameters at the time step
        Xbeta = sum(exp(1i*beta))/N;
        Xalpha = sum(exp(1i*alpha))/N;
        Xtheta = sum(exp(1i*theta))/N;
        Xpsi = sum(exp(1i*psi))/N;
        
        %Calculate the phase of Xalpha to compute the emergent frequency
        diffangle=(cos(angle(Xalpha))-cos(dangle))/sin(dangle);
        dangle=angle(Xalpha);
        
        %Store the order parameters as time series
        Xats_up{kcount}(te) = Xalpha;
        Xbts_up{kcount}(te) = Xbeta;
        Xthetats_up{kcount}(te) = Xtheta;
        Xpsits_up{kcount}(te) = Xpsi;
        
        if(te > (4*Tmax/dt)/5) %time average in final fifth of time series
            Rbetaup_ave = Rbetaup_ave + abs(Xbeta);
            Ralphaup_ave = Ralphaup_ave + abs(Xalpha);
            Rpsiup_ave = Rpsiup_ave + abs(Xpsi);
            Rthetaup_ave = Rthetaup_ave + abs(Xtheta);
            
            %Emergent frequency as average phase of X_\alpha
            D_angle_ave=D_angle_ave+diffangle/dt;
            
            Normcount = Normcount + 1; %update counter
            
        end
        te = te+1; %update time step counter
        
    end
    
    %Normalise appropriately and store the time averages of the order
    %parameters for each coupling strength
    Rbetaup(kcount) = Rbetaup_ave/(Normcount);
    Ralphaup(kcount) = Ralphaup_ave/(Normcount);
    Rthetaup(kcount) = Rthetaup_ave/(Normcount);
    Rpsiup(kcount) = Rpsiup_ave/(Normcount);
    
    Dangleup(kcount)=D_angle_ave/(Normcount);
    
    
    %update counter and store sigma
    kup(kcount) = sigma; 
    kcount = kcount+1; 
    sigma %print coupling strength
    
end

%% Downwards transition operates as upwards, reversing sigma
kcount = 1;
for sigma = sigma_max:-dsigma:0 %sweep sigma down, starting from previous final state
    
    Normcount = 0;
    Rbetadown_ave = 0;
    Ralphadown_ave = 0;
    Rthetadown_ave = 0;
    Rpsidown_ave = 0;
    
    dangle = 0;
    D_angle_ave=0;
    
    
    te = 1;
    for t = 0:dt:Tmax
        
        [k1{1},k1{2}] = func(theta, psi, sigma,w,what);
        [k2{1},k2{2}] = func(theta+0.5*dt*k1{1}, psi+0.5*dt*k1{2},sigma,w,what);
        [k3{1},k3{2}] = func(theta+0.5*dt*k2{1}, psi+0.5*dt*k2{2},sigma,w,what);
        [k4{1},k4{2}] = func(theta+dt*k3{1}, psi+dt*k3{2},sigma,w,what);
        
        theta = theta + (dt/6)*(k1{1}+2*k2{1}+2*k3{1}+k4{1});
        psi = psi + (dt/6)*(k1{2}+2*k2{2}+2*k3{2}+k4{2});
        
        alpha = theta + psi/2;
        beta = (N/(N-1))*(theta-sum(theta)/N) - psi;
        
        Xbeta = sum(exp(1i*beta))/N;
        Xalpha = sum(exp(1i*alpha))/N;
        Xtheta = sum(exp(1i*theta))/N;
        Xpsi = sum(exp(1i*psi))/N;
        
        diffangle=(cos(angle(Xalpha))-cos(dangle))/sin(dangle);
        dangle=angle(Xalpha);
        
        Xats_down{kcount}(te) = Xalpha;
        Xbts_down{kcount}(te) = Xbeta;
        Xthetats_down{kcount}(te) = Xtheta;
        Xpsits_down{kcount}(te) = Xpsi;
        
        
        if(te > 4*(Tmax/dt)/5)
            
            Rbetadown_ave = Rbetadown_ave + abs(Xbeta);
            Ralphadown_ave = Ralphadown_ave + abs(Xalpha);
            Rpsidown_ave = Rpsidown_ave + abs(Xpsi);
            Rthetadown_ave = Rthetadown_ave + abs(Xtheta);
            D_angle_ave=D_angle_ave+diffangle/dt;
            
            Normcount = Normcount +1;
        end
        te = te+1;
        
    end
    
    Rbetadown(kcount) = Rbetadown_ave/(Normcount);
    Ralphadown(kcount) = Ralphadown_ave/(Normcount);
    Rthetadown(kcount) = Rthetadown_ave/(Normcount);
    Rpsidown(kcount) = Rpsidown_ave/(Normcount);
    
    Dangledown(kcount)=D_angle_ave/(Normcount);
    
    kdown(kcount) = sigma;
    kcount = kcount+1;
    sigma
end
%% Plot the phase diagram
figure()

subplot(2,2,1);
plot(kup,Ralphaup,'-','linewidth',1, 'marker','o','markersize', 5, 'displayname', 'R_\alpha forward'); hold on;
plot(kdown, Ralphadown,'-','linewidth',1,'marker','d','markersize', 5,'displayname', 'R_\alpha backward');
set(gca,'Xlim', [0,5],'Ylim',[0,1])
yticks([0:0.2:1])
yticklabels({'0.0','0.2','0.4','0.6','0.8','1.0'})

xlabel('\sigma','Interpreter', 'tex','fontsize', 24, 'fontweight','bold')
ylabel('R_\alpha','Interpreter', 'tex','fontsize', 24, 'fontweight','bold')
legend('location',' best','Interpreter', 'tex','box','off','fontsize', 16,'fontweight','bold')

subplot(2,2,2); plot(kup,Rbetaup,'-','linewidth',1, 'marker','o','markersize', 5, 'displayname', 'R_\beta forward'); hold on;
plot(kdown, Rbetadown,'-','linewidth',1, 'marker','d','markersize', 5,'displayname', 'R_\beta backward');
set(gca,'Xlim', [0,5],'Ylim',[0,1])
yticks([0:0.2:1])
yticklabels({'0.0','0.2','0.4','0.6','0.8','1.0'})


xlabel('\sigma','Interpreter', 'tex','fontsize', 24, 'fontweight','bold')
ylabel('R_{\beta}','Interpreter', 'tex','fontsize', 24, 'fontweight','bold')
legend('location',' best','Interpreter', 'tex','box','off','fontsize', 16,'fontweight','bold')


subplot(2,2,3)
plot(kup,Rthetaup,'-','linewidth',1, 'marker','o','markersize', 5, 'displayname', 'R_\theta forward'); hold on;
plot(kdown, Rthetadown,'-','linewidth',1, 'marker','d','markersize', 5,'displayname', 'R_\theta backward');
set(gca,'Xlim', [0,5],'Ylim',[0,1])
yticks([0:0.2:1])
yticklabels({'0.0','0.2','0.4','0.6','0.8','1.0'})

xlabel('\sigma','Interpreter', 'tex','fontsize', 24, 'fontweight','bold')
ylabel('R_{\theta}','Interpreter', 'tex','fontsize', 24, 'fontweight','bold')
legend('location',' best','Interpreter', 'tex','box','off','fontsize', 16,'fontweight','bold')

subplot(2,2,4)
plot(kup,Rpsiup,'-','linewidth',1, 'marker','o','markersize', 5, 'displayname', 'R_\psi forward'); hold on;
plot(kdown, Rpsidown,'-','linewidth',1, 'marker','d','markersize', 5,'displayname', 'R_\psi backward');
set(gca,'Xlim', [0,5],'Ylim',[0,1])
yticks([0:0.2:1])
yticklabels({'0.0','0.2','0.4','0.6','0.8','1.0'})

xlabel('\sigma','Interpreter', 'tex','fontsize', 24, 'fontweight','bold')
ylabel('R_{\psi}','Interpreter', 'tex','fontsize', 24, 'fontweight','bold','fontweight','bold')
legend('location',' best','Interpreter', 'tex','box','off','fontsize', 16,'fontweight','bold')

%% Store the results and initial conditions
save('Results.mat','Dangledown','Dangleup','Ralphaup','Rbetaup','Ralphadown','Rbetadown', 'Rthetaup','Rthetadown','Rpsiup','Rpsidown','kdown', 'kup','Xats_down','Xbts_down','Xats_up','Xbts_up','Tmax','dt','-v7.3','N');
save('Initial_Conditions.mat','w', 'what', 'theta0','psi0');