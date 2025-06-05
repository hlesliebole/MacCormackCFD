
clear all
close all 


stepx = 20;
stepy = 20; %specify location of the step edge

[xx,yy] = ndgrid (linspace(0,1*10^(-5),75),linspace(0,8*10^(-6),80)); %generate grid
dt = 2.35*10^(-11);

dx = xx(2,1) - xx(1,1);
dy = yy(1,2) - yy(1,1);

% define physical parameters 
M =4;
rho_0 = 1.225; %kg/m^3
R = 287; % J/kg.K
cp = 1005; % J/kg.K
cv = 718; % J/kg.K
gamma = 1.4;
Pr = 0.71; 
Nt = 1500; % number of time steps
Nx = 75;
Ny = 80;

%Boundary condition
T_inf = 288.15; % let the initial condition of the temperature be 288.15 k
u_inf = M*(gamma*R*T_inf)^0.5; 
p_inf = rho_0*R*T_inf;

u = u_inf*ones(Nx,Ny);
v = zeros(Nx,Ny);
p = p_inf*ones(Nx,Ny);
T = T_inf*ones(Nx,Ny);
rho = rho_0*ones(Nx,Ny);

% boundary conditions (rewrite it for clarity)
% wall BC
u(:,1) = 0;
v(:,1) = 0;
T(:,1) = T_inf;
p(:,1) = 2*p(:,2) - p(:,3);
% inlet&far-field BC
u(1,:) = u_inf;
u(:,end) = u_inf;
v(1,:) = 0;
v(:,end) = 0;
p(1,:) = p_inf;
p(:,end) = p_inf;
T(1,:) = T_inf;
T(:,end) = T_inf;
% outlet BC
u(end,:) = 2*u(end-2,:)-u(end-1,:);
v(end,:) = 2*v(end-2,:)-v(end-1,:);
p(end,:) = 2*p(end-2,:)-p(end-1,:);
T(end,:) = 2*T(end-2,:)-T(end-1,:);

% leading edge
% since leading edge cell does not enter compuational domain, the BC is omitted for
% computational efficiency

%obtain conservative var
U_cons = prim2cons(rho,u,v,T,cv);

% allocate cons with zeros
U_cons_tot = zeros(4,Nx,Ny,Nt);
U_cons_pred = zeros(4,Nx,Ny);
E = zeros(4,Nx,Ny);
F = zeros(4,Nx,Ny);


%initialize U, E and F
U_cons_tot(:,:,:,1) = U_cons;

% Compute initial physical parameters
a = (gamma*R*T_inf)^0.5; % speed of sound
mu= sutherland(T_inf); %viscosity 
k = cp*mu/Pr; %thermal conductivity
t_tot = 0;



for tind = 1:Nt


    % =============================> Predictor step (foward in dx and dy) 
    %update physical parameters
    a = (gamma*R*T).^0.5; % speed of sound
    mu= sutherland(T); %viscosity 
    k = cp*mu/Pr; %thermal conductivity

    % ingredients for E(use bwd in dx as fwd for dE/dx)
    dudx_bwd = ddx_bwd(u,dx);
    dudy_central = ddy_central(u,dy);
    dvdx_bwd = ddx_bwd(v,dx);
    dvdy_central = ddy_central(v,dy);
    dTdx_bwd = ddx_bwd(T,dx);
    %dTdy_bwd = ddy_bwd(T,dy);

    tau_xx_pred_E = 2*mu.*(dudx_bwd-(1/3)*(dudx_bwd+dvdy_central));
    %tau_yy_pred = 2*mu(dvdy_central-(1/3)*(dudx_bwd+dvdy_central));
    tau_xy_pred_E = mu.*(dudy_central+dvdx_bwd);
    qx_pred_E = -k.*dTdx_bwd;
    %qy_pred = -k*dTdy_bwd;

    % assemble E
    E(1,:,:) = rho.*u;
    E(2,:,:) = rho.*u.^2+p - tau_xx_pred_E;
    E(3,:,:) = rho.*u.*v - tau_xy_pred_E;
    E(4,:,:) = (squeeze(U_cons(4,:,:))+p).*u - u.*tau_xx_pred_E - v.*tau_xy_pred_E + qx_pred_E;

    % ingredients for F(use bwd in dy as fwd for dF/dy)
    dudx_central = ddx_central(u,dx);
    dudy_bwd = ddy_bwd(u,dy);
    dvdx_central = ddx_central(v,dx);
    dvdy_bwd = ddy_bwd(v,dy);
    dTdy_bwd = ddy_bwd(T,dy);

    tau_yy_pred_F = 2*mu.*(dvdy_bwd-(1/3)*(dudx_central+dvdy_bwd));
    tau_xy_pred_F = mu.*(dudy_bwd+dvdx_central);
    qy_pred_F = -k.*dTdy_bwd;

    % assemble F
    F(1,:,:) = rho.*v;
    F(2,:,:) = rho.*u.*v-tau_xy_pred_F;
    F(3,:,:) = rho.*v.^2+p-tau_yy_pred_F;
    F(4,:,:) = (squeeze(U_cons(4,:,:))+p).*v - u.*tau_xy_pred_F - v.*tau_yy_pred_F + qy_pred_F;

    % compute dE/dx and dF/dy in predictor step (using FORWARD diff)
    for i  = 1:4
        dEdx_pred(i,:,:) = ddx_fwd(squeeze(E(i,:,:)),dx);
        dFdy_pred(i,:,:) = ddy_fwd(squeeze(F(i,:,:)),dy);
    end 

    % time marching U in predictor step
    U_cons_pred = U_cons+dt*(-dEdx_pred-dFdy_pred);

    % update primitive var
    [rho,u,v,T,p,e,Et] = cons2prim(U_cons_pred,R,cv);

    % apply BC
    [u,v,p,T,rho] = applyBC_AdbWall_step(u,v,p,T,u_inf,p_inf,T_inf,stepx,stepy);    

    % update conserv var U_cons_pred
    U_cons_pred = prim2cons(rho,u,v,T,cv);



    %===========================> Corrector step (backward in dx and dy)
    %update physical parameters
    mu= sutherland(T); %viscosity 
    k = cp*mu/Pr; %thermal conductivity


    % ingredients for E(use FORWARD in dx as bwd for dE/dx, CENTRAL in dy)
    dudx_fwd = ddx_fwd(u,dx);
    dudy_central = ddy_central(u,dy);
    dvdx_fwd = ddx_fwd(v,dx);
    dvdy_central = ddy_central(v,dy);
    dTdx_fwd = ddx_fwd(T,dx);

    tau_xx_corr_E = 2*mu.*(dudx_fwd-(1/3)*(dudx_fwd+dvdy_central));
    tau_xy_corr_E = mu.*(dudy_central+dvdx_fwd);
    qx_corr_E = -k.*dTdx_fwd;

    % assemble E (using updated prim var from pred with BC)
    E(1,:,:) = rho.*u;
    E(2,:,:) = rho.*u.^2+p - tau_xx_corr_E;
    E(3,:,:) = rho.*u.*v - tau_xy_corr_E;
    E(4,:,:) = (squeeze(U_cons_pred(4,:,:))+p).*u - u.*tau_xx_corr_E - v.*tau_xy_corr_E + qx_corr_E;

    % ingredients for F(use FORWARD in dy as bwd for dF/dy, CENTRAL for dx)
    dudx_central = ddx_central(u,dx);
    dudy_fwd = ddy_fwd(u,dy);
    dvdx_central = ddx_central(v,dx);
    dvdy_fwd = ddy_fwd(v,dy);
    dTdy_fwd = ddy_fwd(T,dy);

    tau_yy_corr_F = 2*mu.*(dvdy_fwd-(1/3)*(dudx_central+dvdy_fwd));
    tau_xy_corr_F = mu.*(dudy_fwd+dvdx_central);
    qy_corr_F = -k.*dTdy_fwd;

    % assemble F
    F(1,:,:) = rho.*v;
    F(2,:,:) = rho.*u.*v-tau_xy_corr_F;
    F(3,:,:) = rho.*v.^2+p-tau_yy_corr_F;
    F(4,:,:) = (squeeze(U_cons_pred(4,:,:))+p).*v - u.*tau_xy_corr_F - v.*tau_yy_corr_F + qy_corr_F;


    % compute dE/dx and dF/dy in corrector step (using BACKWARD diff)
    for i  = 1:4
        dEdx_corr(i,:,:) = ddx_bwd(squeeze(E(i,:,:)),dx);
        dFdy_corr(i,:,:) = ddy_bwd(squeeze(F(i,:,:)),dy);
    end 


    % time marching U in corrector step
    U_cons = 0.5*(U_cons+U_cons_pred+dt*(-dEdx_corr-dFdy_corr));
    
    % update primitive var
    [rho,u,v,T,p,e,Et] = cons2prim(U_cons,R,cv);

   % apply BC
    [u,v,p,T,rho] = applyBC_AdbWall_step(u,v,p,T,u_inf,p_inf,T_inf,stepx,stepy); 


    % update U_cons
    U_cons = prim2cons(rho,u,v,T,cv);


    % save snapshot to total data 
    U_cons_tot(:,:,:,tind+1) = U_cons;
    t_tot = t_tot+dt;

    %check convergence
    [~,~,~,T_past,~,~,~] = cons2prim(squeeze(U_cons_tot(:,:,:,tind)),R,cv);
    T_change_sum(tind) = mean(T-T_past,"all"); % use average change of temperature 
    % from the current time step to the previous time steps as a metric to check convergence

    % ================================> Visualization
    if mod(tind,50) == 0 
        plot_field(rho,u,v,e,p,T,T_change_sum(1:tind),tind,xx,yy,dt)

    end 


end 

U_cons_final = squeeze(U_cons_tot(:,:,:,end));
[rho,u,v,T,p,e,Et] = cons2prim(U_cons_final,R,cv);
[S] = get_Schlieren(rho,dx,dy);

%% final plot at 1500 time steps


figure()
subplot(321)
pcolor(xx,yy,rho(:,:))
axis equal tight
shading interp
cb = colorbar; 
ylabel(cb,'\rho [kg/m^3]')
title('\rho [kg/m^3]')
xlabel('x(m)','FontSize',16)
ylabel('y(m)','FontSize',16)
%caxis([0,1])

subplot(322)
pcolor(xx,yy,u(:,:))
axis equal tight
shading interp
cb = colorbar; 
ylabel(cb,'u [m/s]')
title('u [m/s]')
xlabel('x(m)','FontSize',16)
ylabel('y(m)','FontSize',16)

subplot(323)
pcolor(xx,yy,v(:,:))
axis equal tight
shading interp
cb = colorbar; 
ylabel(cb,'v [m/s]')
title('v [m/s]')
xlabel('x(m)','FontSize',16)
ylabel('y(m)','FontSize',16)

subplot(324)
pcolor(xx,yy,e(:,:))
axis equal tight
shading interp
cb = colorbar; 
ylabel(cb,'e [m^2/s^2]')
title('e [m^2/s^2]')
xlabel('x(m)','FontSize',16)
ylabel('y(m)','FontSize',16)


subplot(325)
pcolor(xx,yy,p(:,:))
axis equal tight
shading interp
cb = colorbar; 
ylabel(cb,'p [kg/s^2/m]')
title('p [kg/s^2/m]')
xlabel('x(m)','FontSize',16)
ylabel('y(m)','FontSize',16)

subplot(326)
pcolor(xx,yy,T(:,:))
axis equal tight
shading interp
cb = colorbar; 
ylabel(cb,'T [K]')
title('T [K]')
xlabel('x(m)','FontSize',16)
ylabel('y(m)','FontSize',16)

sgtitle('1500 Time Steps','FontSize',20)



%% Functions (new)

function [u,v,p,T,rho] = applyBC_AdbWall_step(u,v,p,T,u_inf,p_inf,T_inf,steploc_x,steploc_y)
    % boundary condition with adiabatic wall 
    R = 287; % J/kg.K

    % enforce BCs on primitive var
    % wall BC
    u(:,1) = 0;
    v(:,1) = 0;
    T(:,1) = T(:,2); % the adiabatic wall BC via dT/dy = 0 at the wall 
    p(:,1) = 2*p(:,2) - p(:,3);
    % inlet&far-field BC
    u(1,:) = u_inf;
    u(:,end) = u_inf;
    v(1,:) = 0;
    v(:,end) = 0;
    p(1,:) = p_inf;
    p(:,end) = p_inf;
    T(1,:) = T_inf;
    T(:,end) = T_inf;
    % outlet BC
    u(end,:) = 2*u(end-1,:)-u(end-2,:);
    v(end,:) = 2*v(end-1,:)-v(end-2,:);
    p(end,:) = 2*p(end-1,:)-p(end-2,:);
    T(end,:) = 2*T(end-1,:)-T(end-2,:);

    % BC for step
    u(steploc_x:end,1:steploc_y) = 0;
    v(steploc_x:end,1:steploc_y) = 0;
    T(steploc_x:end,steploc_y+1) = T(steploc_x:end,steploc_y); % the adiabatic wall BC via dT/dy = 0 at the wall 
    T(steploc_x-1,1:steploc_y) = T(steploc_x,1:steploc_y);
    T(steploc_x-1,steploc_y+1) = T(steploc_x,steploc_y);
    p(steploc_x+1:end,1:steploc_y-1) = 0;
    p(steploc_x:end,steploc_y) = 2*p(steploc_x:end,steploc_y+1) - p(steploc_x:end,steploc_y+2);
    p(steploc_x,1:steploc_y) = 2*p(steploc_x-1,1:steploc_y) - p (steploc_x-2,1:steploc_y);

    

    % use updated uvpT with BC to update rho to make sure consistency of BC
    rho = p./(R*T);
end 


function [u,v,p,T,rho] = applyBC(u,v,p,T,u_inf,p_inf,T_inf)
    R = 287; % J/kg.K

    % enforce BCs on primitive var
    % wall BC
    u(:,1) = 0;
    v(:,1) = 0;
    T(:,1) = T_inf;
    p(:,1) = 2*p(:,2) - p(:,3);
    % inlet&far-field BC
    u(1,:) = u_inf;
    u(:,end) = u_inf;
    v(1,:) = 0;
    v(:,end) = 0;
    p(1,:) = p_inf;
    p(:,end) = p_inf;
    T(1,:) = T_inf;
    T(:,end) = T_inf;
    % outlet BC
    u(end,:) = 2*u(end-1,:)-u(end-2,:);
    v(end,:) = 2*v(end-1,:)-v(end-2,:);
    p(end,:) = 2*p(end-1,:)-p(end-2,:);
    T(end,:) = 2*T(end-1,:)-T(end-2,:);

    % use updated uvpT with BC to update rho to make sure consistency of BC
    rho = p./(R*T);
end 

function [u,v,p,T,rho] = applyBC_AdbWall(u,v,p,T,u_inf,p_inf,T_inf)
    % boundary condition with adiabatic wall 
    R = 287; % J/kg.K

    % enforce BCs on primitive var
    % wall BC
    u(:,1) = 0;
    v(:,1) = 0;
    T(:,1) = T(:,2); % the adiabatic wall BC via dT/dy = 0 at the wall 
    p(:,1) = 2*p(:,2) - p(:,3);
    % inlet&far-field BC
    u(1,:) = u_inf;
    u(:,end) = u_inf;
    v(1,:) = 0;
    v(:,end) = 0;
    p(1,:) = p_inf;
    p(:,end) = p_inf;
    T(1,:) = T_inf;
    T(:,end) = T_inf;
    % outlet BC
    u(end,:) = 2*u(end-1,:)-u(end-2,:);
    v(end,:) = 2*v(end-1,:)-v(end-2,:);
    p(end,:) = 2*p(end-1,:)-p(end-2,:);
    T(end,:) = 2*T(end-1,:)-T(end-2,:);

    % use updated uvpT with BC to update rho to make sure consistency of BC
    rho = p./(R*T);
end 

function [] = plot_field(rho,u,v,e,p,T,T_past,tind,xx,yy,dt)
figure(1)
subplot(331)
pcolor(xx,yy,rho(:,:))
axis equal tight
shading interp
cb = colorbar; 
ylabel(cb,'\rho [kg/m^3]')
title('\rho [kg/m^3]')
xlabel('x(m)','FontSize',16)
ylabel('y(m)','FontSize',16)

subplot(332)
pcolor(xx,yy,u(:,:))
axis equal tight
shading interp
cb = colorbar; 
ylabel(cb,'u [m/s]')
title('u [m/s]')
xlabel('x(m)','FontSize',16)
ylabel('y(m)','FontSize',16)

subplot(333)
pcolor(xx,yy,v(:,:))
axis equal tight
shading interp
cb = colorbar; 
ylabel(cb,'v [m/s]')
title('v [m/s]')
xlabel('x(m)','FontSize',16)
ylabel('y(m)','FontSize',16)

subplot(334)
pcolor(xx,yy,e(:,:))
axis equal tight
shading interp
cb = colorbar; 
ylabel(cb,'e [m^2/s^2]')
title('e [m^2/s^2]')
xlabel('x(m)','FontSize',16)
ylabel('y(m)','FontSize',16)


subplot(335)
pcolor(xx,yy,p(:,:))
axis equal tight
shading interp
cb = colorbar; 
ylabel(cb,'p [kg/s^2/m]')
title('p [kg/s^2/m]')
xlabel('x(m)','FontSize',16)
ylabel('y(m)','FontSize',16)

subplot(336)
pcolor(xx,yy,T(:,:))
axis equal tight
shading interp
cb = colorbar; 
ylabel(cb,'T [K]')
title('T [K]')
xlabel('x(m)','FontSize',16)
ylabel('y(m)','FontSize',16)

subplot(3,3,[7 9])
plot((1:tind)*dt,T_past)
xlabel('t(s)','FontSize',16)
ylabel('$\Delta \bar{T}$ (K)','Interpreter','latex','FontSize',16)
title('Convergence','FontSize',16)

drawnow
end 


function [S] = get_Schlieren(rho,dx,dy)
beta = 0.8;
kappa = 10;
grad_rho = (ddx_central(rho,dx).^2+ddy_central(rho,dy).^2).^(1/2);

S = beta*exp(-kappa*(grad_rho)/(max(grad_rho,[],'all')));
end 


function [T_3loc] = get_3loc(T)
% This function select 3 different x location x/L = 0.25, 0.5, 0.75 of the
% varialbe T
dim = size(T);
loc_ind = [19,38,57]; % the index of x such that xx(loc_ind,1)/L = 0.25, 0.5, 0.75 
T_3loc = zeros(length(loc_ind),dim(2));

for i = 1:length(loc_ind)
    T_3loc(i,:) = T(loc_ind(i),:);
end 

end 

function [T_normalize] = normalize3loc1(T)
% this function normalize the variable T by the maximum T and repeat that
% for all 3 x locaitons
    for i = 1:3
        T_normalize(i,:) = T(i,:) / max(T(i,:));
    end 

end 

function [T_normalize] = normalize3loc2(T)
% this function normalize the variable T by the T_inf and repeat that
% for all 3 x locaitons
    for i = 1:3
        T_normalize(i,:) = T(i,:) / (T(i,end));
    end 

end 

function U_cons_final_adb=runcase_adbwall
% this function essentially run the simulation but change the wall BC to
% Adiabatic wall. The output is the 1500 time step conservative var U

clearvars -except U_cons_final


[xx,yy] = ndgrid (linspace(0,1*10^(-5),75),linspace(0,8*10^(-6),80)); %generate grid
dt = 2.35*10^(-11);

dx = xx(2,1) - xx(1,1);
dy = yy(1,2) - yy(1,1);

% define physical parameters 
M =4;
rho_0 = 1.225; %kg/m^3
R = 287; % J/kg.K
cp = 1005; % J/kg.K
cv = 718; % J/kg.K
gamma = 1.4;
Pr = 0.71; 
Nt = 1500; % number of time steps
Nx = 75;
Ny = 80;

%Boundary condition
T_inf = 288.15; % let the initial condition of the temperature be 288.15 k
u_inf = M*(gamma*R*T_inf)^0.5; 
p_inf = rho_0*R*T_inf;

u = u_inf*ones(Nx,Ny);
v = zeros(Nx,Ny);
p = p_inf*ones(Nx,Ny);
T = T_inf*ones(Nx,Ny);
rho = rho_0*ones(Nx,Ny);

% boundary conditions (rewrite it for clarity)
% wall BC
u(:,1) = 0;
v(:,1) = 0;
T(:,1) = T_inf;
p(:,1) = 2*p(:,2) - p(:,3);
% inlet&far-field BC
u(1,:) = u_inf;
u(:,end) = u_inf;
v(1,:) = 0;
v(:,end) = 0;
p(1,:) = p_inf;
p(:,end) = p_inf;
T(1,:) = T_inf;
T(:,end) = T_inf;
% outlet BC
u(end,:) = 2*u(end-2,:)-u(end-1,:);
v(end,:) = 2*v(end-2,:)-v(end-1,:);
p(end,:) = 2*p(end-2,:)-p(end-1,:);
T(end,:) = 2*T(end-2,:)-T(end-1,:);

% leading edge
% since leading edge cell does not enter compuational domain, the BC is omitted for
% computational efficiency

%obtain conservative var
U_cons = prim2cons(rho,u,v,T,cv);

% allocate cons with zeros
U_cons_tot = zeros(4,Nx,Ny,Nt);
U_cons_pred = zeros(4,Nx,Ny);
E = zeros(4,Nx,Ny);
F = zeros(4,Nx,Ny);


%initialize U, E and F
U_cons_tot(:,:,:,1) = U_cons;

% Compute initial physical parameters
a = (gamma*R*T_inf)^0.5; % speed of sound
mu= sutherland(T_inf); %viscosity 
k = cp*mu/Pr; %thermal conductivity
t_tot = 0;



for tind = 1:Nt


    % =============================> Predictor step (foward in dx and dy) 
    %update physical parameters
    a = (gamma*R*T).^0.5; % speed of sound
    mu= sutherland(T); %viscosity 
    k = cp*mu/Pr; %thermal conductivity

    % ingredients for E(use bwd in dx as fwd for dE/dx)
    dudx_bwd = ddx_bwd(u,dx);
    dudy_central = ddy_central(u,dy);
    dvdx_bwd = ddx_bwd(v,dx);
    dvdy_central = ddy_central(v,dy);
    dTdx_bwd = ddx_bwd(T,dx);
    %dTdy_bwd = ddy_bwd(T,dy);

    tau_xx_pred_E = 2*mu.*(dudx_bwd-(1/3)*(dudx_bwd+dvdy_central));
    %tau_yy_pred = 2*mu(dvdy_central-(1/3)*(dudx_bwd+dvdy_central));
    tau_xy_pred_E = mu.*(dudy_central+dvdx_bwd);
    qx_pred_E = -k.*dTdx_bwd;
    %qy_pred = -k*dTdy_bwd;

    % assemble E
    E(1,:,:) = rho.*u;
    E(2,:,:) = rho.*u.^2+p - tau_xx_pred_E;
    E(3,:,:) = rho.*u.*v - tau_xy_pred_E;
    E(4,:,:) = (squeeze(U_cons(4,:,:))+p).*u - u.*tau_xx_pred_E - v.*tau_xy_pred_E + qx_pred_E;

    % ingredients for F(use bwd in dy as fwd for dF/dy)
    dudx_central = ddx_central(u,dx);
    dudy_bwd = ddy_bwd(u,dy);
    dvdx_central = ddx_central(v,dx);
    dvdy_bwd = ddy_bwd(v,dy);
    dTdy_bwd = ddy_bwd(T,dy);

    tau_yy_pred_F = 2*mu.*(dvdy_bwd-(1/3)*(dudx_central+dvdy_bwd));
    tau_xy_pred_F = mu.*(dudy_bwd+dvdx_central);
    qy_pred_F = -k.*dTdy_bwd;

    % assemble F
    F(1,:,:) = rho.*v;
    F(2,:,:) = rho.*u.*v-tau_xy_pred_F;
    F(3,:,:) = rho.*v.^2+p-tau_yy_pred_F;
    F(4,:,:) = (squeeze(U_cons(4,:,:))+p).*v - u.*tau_xy_pred_F - v.*tau_yy_pred_F + qy_pred_F;

    % compute dE/dx and dF/dy in predictor step (using FORWARD diff)
    for i  = 1:4
        dEdx_pred(i,:,:) = ddx_fwd(squeeze(E(i,:,:)),dx);
        dFdy_pred(i,:,:) = ddy_fwd(squeeze(F(i,:,:)),dy);
    end 

    % time marching U in predictor step
    U_cons_pred = U_cons+dt*(-dEdx_pred-dFdy_pred);

    % update primitive var
    [rho,u,v,T,p,e,Et] = cons2prim(U_cons_pred,R,cv);

    % apply BC
    [u,v,p,T,rho] = applyBC_AdbWall(u,v,p,T,u_inf,p_inf,T_inf);    

    % update conserv var U_cons_pred
    U_cons_pred = prim2cons(rho,u,v,T,cv);



    %===========================> Corrector step (backward in dx and dy)
    %update physical parameters
    mu= sutherland(T); %viscosity 
    k = cp*mu/Pr; %thermal conductivity


    % ingredients for E(use FORWARD in dx as bwd for dE/dx, CENTRAL in dy)
    dudx_fwd = ddx_fwd(u,dx);
    dudy_central = ddy_central(u,dy);
    dvdx_fwd = ddx_fwd(v,dx);
    dvdy_central = ddy_central(v,dy);
    dTdx_fwd = ddx_fwd(T,dx);

    tau_xx_corr_E = 2*mu.*(dudx_fwd-(1/3)*(dudx_fwd+dvdy_central));
    tau_xy_corr_E = mu.*(dudy_central+dvdx_fwd);
    qx_corr_E = -k.*dTdx_fwd;

    % assemble E (using updated prim var from pred with BC)
    E(1,:,:) = rho.*u;
    E(2,:,:) = rho.*u.^2+p - tau_xx_corr_E;
    E(3,:,:) = rho.*u.*v - tau_xy_corr_E;
    E(4,:,:) = (squeeze(U_cons_pred(4,:,:))+p).*u - u.*tau_xx_corr_E - v.*tau_xy_corr_E + qx_corr_E;

    % ingredients for F(use FORWARD in dy as bwd for dF/dy, CENTRAL for dx)
    dudx_central = ddx_central(u,dx);
    dudy_fwd = ddy_fwd(u,dy);
    dvdx_central = ddx_central(v,dx);
    dvdy_fwd = ddy_fwd(v,dy);
    dTdy_fwd = ddy_fwd(T,dy);

    tau_yy_corr_F = 2*mu.*(dvdy_fwd-(1/3)*(dudx_central+dvdy_fwd));
    tau_xy_corr_F = mu.*(dudy_fwd+dvdx_central);
    qy_corr_F = -k.*dTdy_fwd;

    % assemble F
    F(1,:,:) = rho.*v;
    F(2,:,:) = rho.*u.*v-tau_xy_corr_F;
    F(3,:,:) = rho.*v.^2+p-tau_yy_corr_F;
    F(4,:,:) = (squeeze(U_cons_pred(4,:,:))+p).*v - u.*tau_xy_corr_F - v.*tau_yy_corr_F + qy_corr_F;


    % compute dE/dx and dF/dy in corrector step (using BACKWARD diff)
    for i  = 1:4
        dEdx_corr(i,:,:) = ddx_bwd(squeeze(E(i,:,:)),dx);
        dFdy_corr(i,:,:) = ddy_bwd(squeeze(F(i,:,:)),dy);
    end 


    % time marching U in corrector step
    U_cons = 0.5*(U_cons+U_cons_pred+dt*(-dEdx_corr-dFdy_corr));
    
    % update primitive var
    [rho,u,v,T,p,e,Et] = cons2prim(U_cons,R,cv);

   % apply BC
    [u,v,p,T,rho] = applyBC_AdbWall(u,v,p,T,u_inf,p_inf,T_inf); 


    % update U_cons
    U_cons = prim2cons(rho,u,v,T,cv);


    % save snapshot to total data 
    U_cons_tot(:,:,:,tind+1) = U_cons;
    t_tot = t_tot+dt;


end 


% grab var at 1500 time step
U_cons_final_adb = squeeze(U_cons_tot(:,:,:,end));

clearvars -except U_cons_final_adb

end 


%% function (previous homeworks/attached here in case the files in Zip file doesn't work)
function [mu] = sutherland(T)
mu0 = 1.735*10^(-5); %N.S/m^2
T0 = 288.15; %K
S1 = 110.4; % K
mu = mu0.*(T/T0).^(3/2).*(T0+S1)./(T+S1);
end 

function [rho,u,v,T,p,e,Et] = cons2prim(U,R,cv)
rho = squeeze(U(1,:,:));
u = squeeze(U(2,:,:))./rho;
v = squeeze(U(3,:,:))./rho;
Et = squeeze(U(4,:,:));
e = Et./rho - 0.5*(u.^2+v.^2);
T = e/cv;
p = rho.*R.*T;
end 

function U = prim2cons(rho,u,v,T,cv)
U(1,:,:) = rho;
U(2,:,:) = rho.*u;
U(3,:,:) = rho.*v;
Et = rho.*(cv*T+0.5*((u.^2+v.^2)));
U(4,:,:) = Et;
end 

function df = ddx_bwd(f,dx,BC)
if nargin<3
    BC= 'one-sided';
end 
switch BC
    case 'periodic'
    dim = size(f);
    Nx= dim(1); %dimensional in x direction --> using ndgrid
    A = -1*diag(ones(1,Nx-1),-1)+diag(ones(1,Nx));
    A(1,end) = -1;
    df = (1/dx)*A*f;
    otherwise
    dim = size(f);
    Nx= dim(1); %dimensional in x direction --> using ndgrid
    A = -1*diag(ones(1,Nx-1),-1)+diag(ones(1,Nx));
    A(1,1) = -1; %forward diff for first pt 
    A(1,2) = 1;
    df = (1/dx)*A*f;
end 
end 

function df = ddx_central(f,dx,BC)
if nargin<3
    BC= 'one-sided';
end 
switch BC
    case 'periodic'

    dim = size(f);
    Nx= dim(1); %dimensional in x direction --> using ndgrid
    A = -1*diag(ones(1,Nx-1),-1)+diag(ones(1,Nx-1),1);
    A(1,end)=-1;%periodic BC
    A(end,1) = 1;
    df = (1/(2*dx))*A*f;       
    otherwise
    dim = size(f);
    Nx= dim(1); %dimensional in x direction --> using ndgrid
    A = -1*diag(ones(1,Nx-1),-1)+diag(ones(1,Nx-1),1);
    A(1,1:3)=[-3 4 -1]; %use forward diff (2nd order) for first pt
    df = (1/(2*dx))*A*f;   
end 
end 


function df = ddx_fwd(f,dx,BC)
if nargin<3
    BC= 'one-sided';
end 
switch BC
    case 'periodic'
    dim = size(f);
    Nx= dim(1); %dimensional in x direction --> using ndgrid
    A = diag(ones(1,Nx-1),1)-diag(ones(1,Nx));
    A(end,1) = 1; % periodic BC
    df = (1/dx)*A*f;
    otherwise
    dim = size(f);
    Nx= dim(1); %dimensional in x direction --> using ndgrid
    A = diag(ones(1,Nx-1),1)-diag(ones(1,Nx));
    A(end,end-1) = -1; %use backward for last pt
    A(end,end) = 1;%use backward for last pt
    df = (1/dx)*A*f;

end 
end 


function df = ddy_bwd(f,dx,BC)
if nargin<3
    BC= 'one-sided';
end 
switch BC
    case 'periodic'
    f_T = f';
    dim = size(f_T);
    Nx= dim(1); %dimensional in x direction --> using ndgrid
    A = -1*diag(ones(1,Nx-1),-1)+diag(ones(1,Nx));
    A(1,end) = -1;
    df = ((1/dx)*A*f_T)';
    otherwise
    f_T = f';
    dim = size(f_T);
    Nx= dim(1); %dimensional in x direction --> using ndgrid
    A = -1*diag(ones(1,Nx-1),-1)+diag(ones(1,Nx));
    A(1,1) = -1; %forward diff for first pt 
    A(1,2) = 1;
    df = ((1/dx)*A*f_T)';

end 
end 


function df = ddy_central(f,dx,BC)
if nargin<3
    BC= 'one-sided';
end 
switch BC
    case 'periodic'
    f_T = f';
    dim = size(f_T);
    Nx= dim(1); %dimensional in x direction --> using ndgrid
    A = -1*diag(ones(1,Nx-1),-1)+diag(ones(1,Nx-1),1);
    A(1,end)=-1; %use forward diff (2nd order) for first pt
    A(end,1) = 1;%use backward diff (2nd order) for last pt
    df = ((1/(2*dx))*A*f_T)';       
    otherwise
    f_T = f';
    dim = size(f_T);
    Nx= dim(1); %dimensional in x direction --> using ndgrid
    A = -1*diag(ones(1,Nx-1),-1)+diag(ones(1,Nx-1),1);
    A(1,1:3)=[-3 4 -1]; %use forward diff (2nd order) for first pt
    df = ((1/(2*dx))*A*f_T)'; 

end
end 

function df = ddy_fwd(f,dx,BC)
if nargin<3
    BC= 'one-sided';
end 
switch BC
    case 'periodic'
    f_T = f';
    dim = size(f_T);
    Nx= dim(1); %dimensional in x direction --> using ndgrid
    A = diag(ones(1,Nx-1),1)-diag(ones(1,Nx));
    A(end,1) = 1; % periodic BC
    df = ((1/dx)*A*f_T)';
    otherwise
    f_T = f';
    dim = size(f_T);
    Nx= dim(1); %dimensional in x direction --> using ndgrid
    A = diag(ones(1,Nx-1),1)-diag(ones(1,Nx));
    A(end,end-1) = -1; %use backward for last pt
    A(end,end) = 1;%use backward for last pt
    df = ((1/dx)*A*f_T)';

end 
end 



