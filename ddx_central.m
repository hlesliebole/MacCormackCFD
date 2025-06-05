function dfdx = ddx_central(f, dx)
    [nx, ny] = size(f);
    dfdx = zeros(nx, ny);
    
    % apply central difference for i = 2 to nx-1
    dfdx(2:nx-1,:) = (f(3:nx,:) - f(1:nx-2,:))/(2*dx);
    
    % at i=1, use second order forward difference
    dfdx(1,:) = (-3*f(1,:) + 4*f(2,:) - f(3,:))/(2*dx);

    % at i=nx, use second order backward difference
    dfdx(nx,:) = (3*f(nx,:) - 4*f(nx-1,:) + f(nx-2,:))/(2*dx);
end