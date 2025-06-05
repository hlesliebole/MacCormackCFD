function [p, rho, u, v, T] = applyBCs(rho, u, v, T, R, T_inf, p, p_inf, u_inf,thermalRegime)
% applyBCs.m - Apply BCs to primitive variables with enhanced numerical stability
% Enforces no-slip wall (bottom), freestream (top/left), outflow (right)

[nx, ny] = size(rho);

% --- Bottom wall (j = 1): no-slip, isothermal ---
u(:,1) = 0;
v(:,1) = 0;
if thermalRegime == 1 % iso
    T(:,1) = T_inf;
elseif thermalRegime == 2 % adi
    T(:,1)=T(:,2);
else
    error('Mismatched thermal regime')
end
p(:,1) = 2*p(:,2)-p(:,3);

% --- Top boundary (j = ny): freestream ---
u(:,ny)   = u_inf;
v(:,ny)   = 0;
T(:,ny)   = T_inf;
p(:,ny)   = p_inf;

% --- Left boundary (i = 1): freestream ---
u(1,:)   = u_inf;          
v(1,:)   = 0;
T(1,:)   = T_inf;          
p(1,:)   = p_inf;

% --- Right boundary (i = nx): sub-sonic non-reflective ---
u(nx,:)   = 2*u(nx-1,:) - u(nx-2,:);
v(nx,:)   = 2*v(nx-1,:) - v(nx-2,:);
T(nx,:)   = 2*T(nx-1,:) - T(nx-2,:);
p(nx,:)   = 2*p(nx-1,:) - p(nx-2,:);

rho = p./R./T;
end