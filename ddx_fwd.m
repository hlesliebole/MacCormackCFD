function dfdx = ddx_fwd(f, dx , ~)
% F is [nx, ny, nvar]
    % dfdx = (circshift(f, [-1, 0, 0]) - f) / dx;    

    [nx, ny, nv] = size(f);
    dfdx = zeros(nx, ny, nv);

    % apply forward difference for i = 1 to nx-1
    dfdx(1:nx-1,:,:) = (f(2:nx,:,:) - f(1:nx-1,:,:))/dx;

    % at i=nx, use backward difference
    %dfdx(nx,:,:) = (f(nx,:,:) - f(nx-1,:,:))/dx;
    % or copy last interior gradient?
    dfdx(nx,:,:)     = dfdx(nx-1,:,:);          % copy last interior gradient

end