function dfdy = ddy_central(f, dy)
    [nx, ny] = size(f);
    dfdy = zeros(nx, ny);
    
    % Central difference: second-order accurate interior (j = 2 to ny-1)
    dfdy(:, 2:ny-1) = (f(:, 3:ny) - f(:, 1:ny-2)) / (2 * dy);

    % Left boundary (j = 1): second-order forward difference
    dfdy(:, 1) = (-3*f(:,1) + 4*f(:,2) - f(:,3)) / (2 * dy);

    % Right boundary (j = ny): second-order backward difference
    dfdy(:, ny) = (3*f(:,ny) - 4*f(:,ny-1) + f(:,ny-2)) / (2 * dy);
end