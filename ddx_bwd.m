function dfdx = ddx_bwd(f, dx, ~)
    % dfdx = (f - circshift(f, [1, 0, 0])) / dx;
    [nx, ny, nv] = size(f);
    dfdx = zeros(nx, ny, nv);

    % apply backward difference for i = 2 to nx
    dfdx(2:nx,:,:) = (f(2:nx,:,:) - f(1:nx-1,:,:))/dx;

    % at i=1, use forward difference
    %dfdx(1,:,:) = (f(2,:,:) - f(1,:,:))/dx;
    dfdx(1,:,:)    = dfdx(2,:,:);               % copy first interior gradient
end