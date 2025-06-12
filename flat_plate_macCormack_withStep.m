%% flat_plate_macCormack_withStep.m – Mach-4 flat-plate solver (compressible Navier-Stokes)
% ----------------------------------------------------------------------------------
% Numerical method  : 2-D MacCormack (predictor / corrector)
% Governing eqns    : Compressible Navier-Stokes with Sutherland viscosity
% Grid              : Uniform x–grid, uniform y–grid (stretching template left in)
% Temporal scheme   : Explicit, variable Δt from CFL + viscous constraint
% Boundary-conds    : Bottom wall (no-slip, isothermal), top/left freestream,
%                     right non-reflecting outflow (see applyBCs.m)
% Visualization     : 6 scalar fields every 10 steps + running residual

clear;  clc;  close all     % fresh workspace

%% 1. Physical / simulation parameters ------------------------------------
thermalRegime = 1; %1 for isothermal, 2 for adiabatic
gamma  = 1.40;                       % ratio of specific heats
cp     = 1005;                       % J/(kg·K)
cv     = cp/gamma;                   % from cp = cv + R
R      = cp - cv;                    % specific gas constant (air)
Pr     = 0.72;                       % Prandtl number (air, ≈ constant)
nu_ref = 1.8e-5;                     % reference kinematic viscosity (unused)
M      = 4;                          % freestream Mach number
T_inf  = 300;                        % K (freestream temperature)
p_inf  = 101300;                     % Pa (freestream pressure)

rho_inf = 1.225;                     % kg/m³  (freestream density)
a_inf   = sqrt(gamma*R*T_inf);       % speed of sound
u_inf   = M * a_inf;                 % freestream x-velocity
v_inf   = 0;                         % freestream y-velocity

%% 2. Grid generation ------------------------------------------------------
L  = 1e-5;    H  = 8e-6;             % domain length / height
nx = 75;      ny = 80;               % grid nodes

% --- y-stretching --------------------
% beta  = 1.2;
% y_u   = linspace(0,1,ny);
% y_s   = (exp(beta*y_u) - 1) / (exp(beta) - 1);
%y_coords = y_s * H;
beta   = -2.0;               % stretching strength (tune between ~1–5)
eta   = linspace(0,1,ny);  % uniform parameter
y_coords = H * tanh(beta*eta)/tanh(beta);

%y_coords = linspace(0, H, ny);       % uniform in y
x_coords = linspace(0, L, nx);       % uniform in x
[x, y]   = ndgrid(x_coords, y_coords);

dx = L /(nx-1);
dy = H /(ny-1);
dy_avg = dy;                         % used only for reference

% --- forward‐facing step geometry ---
x_step = 0.5 * L;
h_step = 0.2 * H;
[~, i_step] = min(abs(x_coords - x_step));
[~, j_step] = min(abs(y_coords - h_step));
solidMask = false(nx, ny);
for ii = i_step:nx
    for jj = 1:(j_step-1)
        solidMask(ii,jj) = true;
    end
end
%% 3. Time-stepping controls ----------------------------------------------
CFL     = 0.3;                       % Courant number (reduced for margin)
dt_init = 2.35e-11/2;                  % s, first step guess
dt = dt_init;
t_end   = 1e-6;                      % total simulated time
nsteps  = 1500; %ceil(t_end/dt_init);%       % upper-bound on iterations
t_tot   = 0;                        % initialize total time counter
%res_rho = zeros(1, nsteps);          % L2 residual history

%% 4. Initial flow field (uniform freestream) -----------------------------
rho = rho_inf * ones(nx,ny);
u   = u_inf  * ones(nx,ny);
v   = zeros(nx,ny);
T   = T_inf  * ones(nx,ny);
p   = p_inf  * ones(nx,ny);

T_prev = T;          % stores ρ from the previous step for convergence


T(T < 1e-3) = 1e-3;   % temperature floor
mu = sutherland(T);                  % dynamic viscosity [Pa·s]
k  = cp .* mu / Pr;                  % thermal conductivity [W/(m·K)]

U  = prim2cons(rho,u,v,T,cv);        % conservative vector [4×nx×ny]

%% 5. Figure style defaults (LaTeX, CMU Serif etc.) -----------------------
set(groot,'DefaultFigureUnits','centimeters','DefaultFigurePosition',[0 0 16 12]);
set(groot,'DefaultAxesFontSize',12,'DefaultAxesLineWidth',1.5,'DefaultLineLineWidth',1.5);
set(groot,'DefaultAxesTickLabelInterpreter','latex','DefaultTextInterpreter','latex');
set(groot,'DefaultAxesFontName','CMU Serif','DefaultLegendInterpreter','latex');
set(groot,'DefaultColorbarTickLabelInterpreter','latex','DefaultFigureColormap',jet);
set(groot,'DefaultFigureRenderer','painters');    % vector graphics
fig1 = figure('Units','normalized','OuterPosition',[0 0 1 1]); % full-screen
gifFileName = ['flat_plate_iso.gif' 'flat_plate_adi.gif'];
%% 6. MAIN TIME LOOP ------------------------------------------------------
%for thermalRegime=1:2
% gifFile = gifFileName(thermalRegime);
% if exist(gifFile,'file'), delete(gifFile); end
for n = 1:nsteps

    % ---------------------------------------------------------------------
    % 6.1  Variable timesteps (explicit CFL + viscous criterion)
    % ---------------------------------------------------------------------
    a  = sqrt(gamma*R.*T);           % local sound speed
    T(T < 1e-3) = 1e-3;   % temperature floor
    mu = sutherland(T);              % update μ & k with current T
    k  = cp .* mu / Pr;

    dt_conv = CFL ./ ( abs(u)/dx + abs(v)/dy + a .* sqrt(1/dx^2 + 1/dy^2));
    dt_visc = 0.25 * min(dx,dy)^2 / ( max(mu(:)) / min(rho(:)) );
    %dt      = min( min(dt_conv(:)), dt_visc );   % global Δt
    t_tot   = t_tot+dt;

    % --- console progress ------------------------------------------------
    fprintf(['step %4d | dt=%8.2e | ρ[min/max]=[%8.2e %8.2e]  ' ...
        'p[min/max]=[%8.2e %8.2e]  T[min/max]=[%8.2e %8.2e]\n'],...
        n, dt, min(rho(:)), max(rho(:)), min(p(:)), max(p(:)), ...
        min(T(:)), max(T(:)));

    % ---------------------------------------------------------------------
    % 6.2  Predictor stage  (forward differences in both directions)
    % ---------------------------------------------------------------------
    % --- velocity / temperature gradients --------------------------------
    dudx_bwd    = ddx_bwd   (u, dx);
    dudy_central= ddy_central(u, dy);
    dvdx_bwd    = ddx_bwd   (v, dx);
    dvdy_central= ddy_central(v, dy);
    dTdx_bwd    = ddx_bwd   (T, dx);

    % --- viscous stresses & heat flux for E-flux -------------------------
    tau_xx_E = 2*mu .* (dudx_bwd - (dudx_bwd+dvdy_central)/3);
    tau_xy_E =     mu .* (dudy_central + dvdx_bwd);
    qx_E     = -k  .*  dTdx_bwd;

    % --- assemble E (convective + viscous) -------------------------------
    E = zeros(4,nx,ny);
    E(1,:,:) = rho.*u;
    E(2,:,:) = rho.*u.^2 + p - tau_xx_E;
    E(3,:,:) = rho.*u.*v       - tau_xy_E;
    E(4,:,:) = (squeeze(U(4,:,:))+p).*u - u.*tau_xx_E - v.*tau_xy_E + qx_E;

    % --- analogous build for F ------------------------------------------
    dudx_central= ddx_central(u, dx);
    dudy_bwd    = ddy_bwd   (u, dy);
    dvdx_central= ddx_central(v, dx);
    dvdy_bwd    = ddy_bwd   (v, dy);
    dTdy_bwd    = ddy_bwd   (T, dy);

    tau_yy_F = 2*mu .* (dvdy_bwd - (dudx_central+dvdy_bwd)/3);
    tau_xy_F =     mu .* (dudy_bwd + dvdx_central);
    qy_F     = -k  .*  dTdy_bwd;

    F = zeros(4,nx,ny);
    F(1,:,:) = rho.*v;
    F(2,:,:) = rho.*u.*v       - tau_xy_F;
    F(3,:,:) = rho.*v.^2 + p   - tau_yy_F;
    F(4,:,:) = (squeeze(U(4,:,:))+p).*v - u.*tau_xy_F - v.*tau_yy_F + qy_F;

    % --- spatial derivatives (forward) -----------------------------------
    dEdx_fwd = zeros(size(E));   dFdy_fwd = zeros(size(F));
    for iVar = 1:4
        dEdx_fwd(iVar,:,:) = ddx_fwd(squeeze(E(iVar,:,:)), dx);
        dFdy_fwd(iVar,:,:) = ddy_fwd(squeeze(F(iVar,:,:)), dy);
    end

    % --- predictor update ------------------------------------------------
    %U_star = U - dt * solidMask.*(dEdx_fwd + dFdy_fwd);
    % Form the “fluid” update only where solidMask == false
    U_star = U;   % start by copying everything
    for iv = 1:4
        for i = 1:nx
            for j = 1:ny
                if ~solidMask(i,j)
                    U_star(iv,i,j) = U(iv,i,j) - dt*( dEdx_fwd(iv,i,j) + dFdy_fwd(iv,i,j) );
                else
                    % leave U_star(iv,i,j) = U(iv,i,j)  (so solid never changes)
                end
            end
        end
    end

    % --- convert to primitives, apply BCs, rebuild U* ---------------------
    [rho,u,v,T,p,e,Et] = cons2prim(U_star, R, cv);
    [p,rho,u,v,T]     = applyBCstep(rho,u,v,T,R,T_inf,p,p_inf,u_inf,i_step,j_step,solidMask);
    T(T < 1e-3) = 1e-3;   % temperature floor
    mu = sutherland(T);  k = cp.*mu/Pr;         % update transport props
    U_star = prim2cons(rho,u,v,T,cv);

    % ---------------------------------------------------------------------
    % 6.3  Corrector stage  (backward differences)
    % ---------------------------------------------------------------------
    % --- gradients for corrected fluxes ----------------------------------
    dudx_fwd    = ddx_fwd   (u, dx);
    dudy_central= ddy_central(u, dy);
    dvdx_fwd    = ddx_fwd   (v, dx);
    dvdy_central= ddy_central(v, dy);
    dTdx_fwd    = ddx_fwd   (T, dx);

    tau_xx_E = 2*mu .* (dudx_fwd - (dudx_fwd+dvdy_central)/3);
    tau_xy_E =     mu .* (dudy_central + dvdx_fwd);
    qx_E     = -k  .*  dTdx_fwd;

    E(1,:,:) = rho.*u;
    E(2,:,:) = rho.*u.^2 + p - tau_xx_E;
    E(3,:,:) = rho.*u.*v       - tau_xy_E;
    E(4,:,:) = (squeeze(U_star(4,:,:))+p).*u - u.*tau_xx_E - v.*tau_xy_E + qx_E;

    dudx_central= ddx_central(u, dx);
    dudy_fwd    = ddy_fwd   (u, dy);
    dvdx_central= ddx_central(v, dx);
    dvdy_fwd    = ddy_fwd   (v, dy);
    dTdy_fwd    = ddy_fwd   (T, dy);

    tau_yy_F = 2*mu .* (dvdy_fwd - (dudx_central+dvdy_fwd)/3);
    tau_xy_F =     mu .* (dudy_fwd + dvdx_central);
    qy_F     = -k  .*  dTdy_fwd;

    F(1,:,:) = rho.*v;
    F(2,:,:) = rho.*u.*v       - tau_xy_F;
    F(3,:,:) = rho.*v.^2 + p   - tau_yy_F;
    F(4,:,:) = (squeeze(U_star(4,:,:))+p).*v - u.*tau_xy_F - v.*tau_yy_F + qy_F;

    % --- spatial derivatives (backward) ----------------------------------
    dEdx_bwd = zeros(size(E));   dFdy_bwd = zeros(size(F));
    for iVar = 1:4
        dEdx_bwd(iVar,:,:) = ddx_bwd(squeeze(E(iVar,:,:)), dx);
        dFdy_bwd(iVar,:,:) = ddy_bwd(squeeze(F(iVar,:,:)), dy);
    end

    % --- corrector update (MacCormack) -----------------------------------
    %U_new = 0.5 * ( U + U_star - dt*solidMask.*(dEdx_bwd + dFdy_bwd) );
    U_new = U;   % again copy everything
    for iv = 1:4
        for i = 1:nx
            for j = 1:ny
                if ~solidMask(i,j)
                    U_new(iv,i,j) = 0.5*( U(iv,i,j) + U_star(iv,i,j) ...
                        - dt*( dEdx_bwd(iv,i,j) + dFdy_bwd(iv,i,j) ) );
                else
                    % solid cell unchanged
                end
            end
        end
    end

    % --- final primitive conversion + BCs + safety check -----------------
    [rho,u,v,T,p,e,Et] = cons2prim(U_new, R, cv);
    [p,rho,u,v,T]     = applyBCstep(rho,u,v,T,R,T_inf,p,p_inf,u_inf,i_step,j_step,solidMask);

    if any(T(:) < 0) || any(~isfinite(T(:)))
        error('Temperature blew up at step %d', n);
    end

    U_new = prim2cons(rho,u,v,T,cv);

    % ---------------------------------------------------------------------
    % 6.4  Convergence monitoring & graphics ------------------------------

    [~,~,~,T_past,~,~,~] = cons2prim(U,R,cv);
    T_change_sum(n) = mean(T-T_past,"all"); % use average change of
    % temperature from the current time step to the previous time steps as
    % a metric to check convergence

    U = U_new;                       % advance to next time level

    if mod(n,50)==0 || n==1
        % dvdx = ddx_central(v, dx);                  % ∂v/∂x
        % dudy = ddy_central(u, dy);                  % ∂u/∂y
        % omega    = dvdx - dudy;                         % vorticity

        figure(fig1); clf;
        tl = tiledlayout(fig1,3,3,'TileSpacing','compact','Padding','compact');

        % ---- hard-wired limits so color ranges never shift -------------
        cLims = [ ...
            0.998,       3.54    ;   % rho
            0,           1.3891e3;   % u
            0,           180.137 ;   % v
            2.0381e5,    3.6489e5;   % e
            1.01e5,      3.05e5 ;    % p
            300,         600     ;
            -1e+10, 1e+10];  % T

        data  = {rho, u, v, e, p, T};
        names = {'$\rho$ ($kg/m^3$)','u ($m/s$)','v ($m/s$)','e ($m^2/s^2$)','p ($Pa$)','T ($K$)'};

        % ---- first six tiles: field contours ----------------------------
        for k = 1:6
            ax = nexttile(tl,k);
            data{k}(solidMask) = NaN;
            pcolor(ax,x,y,data{k});
            shading(ax,'interp');
            axis(ax,'equal','tight');
            clim(ax,cLims(k,:));
            colorbar(ax);
            title(ax,names{k},'Interpreter','latex');
        end

        % ---- bottom-row (tiles 7-9 merged): convergence history ---------
        axConv = nexttile(tl,7,[1 3]);   % span entire bottom row
        semilogy(axConv,1:n,T_change_sum(1:n),'LineWidth',1.2);
        ylim([1e-4 0.6])
        xlim([0 nsteps])
        xlabel(axConv,'Time step');
        ylabel(axConv,'$\Delta \bar{T}$ (K)','Interpreter','latex');
        grid(axConv,'on');
        title(axConv,'Temperature convergence','Interpreter','latex');

        % --- LaTeX-compatible time string --------------------------------
        expT  = floor(log10(t_tot));        % exponent (…-8, -5,  0, 3, …)
        coefT = t_tot/10^expT;              % coefficient in [1,10)

        timeStr = sprintf('$t = %.3g\\,\\times\\,10^{%d}\\;\\mathrm{s}$', ...
            coefT, expT);

        % --- global title ------------------------------------------------
        sgtitle(tl, { ...
            '\bf Mach 3 Compressible Navier--Stokes MacCormack Solver', ...
            timeStr}, ...
            'Interpreter','latex');

        drawnow

        % --- GIF capture & write ---
        gifFile = 'step05x02ymach3.gif';
        frame   = getframe(fig1);     % grab the full window
        im      = frame2im(frame);
        [A,map] = rgb2ind(im,256);

        if n == 1                    % first frame: create / overwrite
            imwrite(A,map,gifFile,'gif','LoopCount',Inf,'DelayTime',0.10);
        else                         % subsequent frames: append
            imwrite(A,map,gifFile,'gif','WriteMode','append','DelayTime',0.10);
        end
    end
    % if thermalRegime ==1
    %     rho_iso=rho;
    %     u_iso=u;
    %     v_iso=v;
    %     T_iso=T;
    %     p_iso=p;
    % elseif thermalRegime==2
    %     rho_adi=rho;
    %     u_adi=u;
    %     v_adi=v;
    %     T_adi=T;
    %     p_adi=p;
    % end
end
%end

%% ------------------------------------------------------------------------
%  Numerical Schlieren of the converged solution
% -------------------------------------------------------------------------
% 1.  ∂ρ/∂x and ∂ρ/∂y with existing central-difference helpers
drhodx = ddx_central(rho,dx);
drhody = ddy_central(rho,dy);

% 2.  Gradient magnitude |∇ρ|
absGrad = hypot(drhodx,drhody);          % hypot = sqrt(a.^2 + b.^2)

% 3.  Scaled Schlieren field  S(x,y)  ∈ (0,β]
beta   = 0.8;
kappa  = 10;
S      = beta * exp( -kappa * absGrad / max(absGrad(:)) );

% 4.  Plot (monochrome schlieren-like image)
figure('Name','Numerical Schlieren','Units','normalized','OuterPosition',[0 0 1 1])

p = pcolor(x,y,S);         shading interp
axis equal tight
set(gca,'YDir','normal')    % keep y increasing upward
colormap(flipud(gray))      % white = weak gradient, black = strong
clim([0 1])                % lock color scale
colorbar
title('Numerical Schlieren (density gradient)','Interpreter','latex')
xlabel('$x$ [m]','Interpreter','latex')
ylabel('$y$ [m]','Interpreter','latex')

%% ------------------------------------------------------------------------
% plot velocity sections
y_max= 10;
x_u(1:3) = x([25,30,35]) ;
u_vert(:,1) = u(25,1:y_max);
u_vert(:,2) = u(30,1:y_max);
u_vert(:,3) = u(35,1:y_max);
y_u = y(1,1:y_max);
figure()
plot(u_vert,y_u,'LineWidth',1.5);
legend('$x=0.32*10^{-5}$','$x=0.39*10^{-5}$','$x=0.46*10^{-5}$','interpreter','latex','Fontsize',20)
xlabel('u(m/s)','FontSize',26)
ylabel('y(m)','FontSize',26)

figure()
pcolor(x,y,u(:,:))
hold on 
plot([x_u(1),x_u(1)],[y_u(1),y_u(end)],'LineWidth',1.5,'Color','k')
hold on 
plot([x_u(2),x_u(2)],[y_u(1),y_u(end)],'LineWidth',1.5,'Color','k')
hold on 
plot([x_u(3),x_u(3)],[y_u(1),y_u(end)],'LineWidth',1.5,'Color','k')
axis equal tight
shading interp
cb = colorbar; 
ylabel(cb,'u [m/s]','FontSize',20)
title('u [m/s]','FontSize',30)
xlabel('x(m)','FontSize',26)

%% --- Mach Angle ---

% --- 1. find theoretical Mach angle ---
thetaTh = asin(1/M);               % radians
thetaThDeg = rad2deg(thetaTh);

% --- 2. regenerate schlieren image ---
figure; clf;
pcolor(x, y, S); shading interp; axis equal tight;
colormap(flipud(gray))      % white = weak gradient, black = strong
clim([0 1])                % lock color scale
colorbar
hold on;
xlabel('$x$ [m]','Interpreter','latex')
ylabel('$y$ [m]','Interpreter','latex')

% --- 3. find the numerical Mach angle
bufX = round(0.05*nx);                 % ignore 5 % columns on each side
bufY = round(0.05*ny);                 % ignore 5 % rows top & bottom
cols  = (bufX+1):(nx-bufX);            % allowed i-indices
rows  = (bufY+1):(ny-bufY);            % allowed j-indices

% pick shock location (min-S) in every interior column
xShock = x(cols,1);                    % physical x for those columns
yShock = zeros(size(xShock));          % to be filled

for k = 1:numel(cols)
    i = cols(k);
    [~,relJ] = min( S(i,rows) );       % row index *within* rows-subarray
    j        = rows(relJ);             % absolute j-index
    yShock(k)= y(i,j);                 % store physical y
end

% discard first/last few points since shock still slightly curved
xShock = xShock(3:end-2);
yShock = yShock(3:end-2);

% least-squares line:  y = m·x + b
p        = polyfit(xShock, yShock, 1); % p(1)=slope=m
thetaNum = atan(abs(p(1)));            % radians
thetaDeg = rad2deg(thetaNum);

% --- 3. draw the Mach angle rays ---
xMax = max(x(:));
yRayTh = tan(thetaTh)  * xMax;
yRayNu = tan(thetaNum) * xMax;

% LaTeX-ready labels
lblTh = sprintf('Theoretical θ =  %.2fº', thetaThDeg);
lblNu = sprintf('Numerical θ = %.2fº', thetaDeg);

plot([0, xMax], [0, yRayTh],  'r-',  'LineWidth',1.8, ...
     'DisplayName', lblTh);
plot([0, xMax], [0, yRayNu], 'b--', 'LineWidth',1.8, ...
     'DisplayName', lblNu);

title(sprintf('Numerical Schlieren with Mach-θ Overlay at Mach %.d',M),'Interpreter','latex');
lgd = legend('Location','best');
set(lgd, 'Interpreter','latex');
hold off;

%% -------------------------- Helper functions ----------------------------
function mu = sutherland(T)
% Dynamic viscosity of air via Sutherland (all units SI)
mu0 = 1.735e-5;   T0 = 288.15;   S = 110.4;
mu  = mu0 * (T/T0).^(3/2) .* (T0+S) ./ (T+S);
mu  = max(mu, 1e-10);            % numerical floor
end