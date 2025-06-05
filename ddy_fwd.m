function dfdy = ddy_fwd(f, dy, ~)
% F is [nx, ny, nvar]
    % dfdy = (circshift(f, [0, -1, 0]) - f) / dy;
    [nx, ny, nv] = size(f);
    dfdy = zeros(nx, ny, nv);

    % Forward difference in y: j = 1 to ny-1
    dfdy(:, 1:ny-1, :) = (f(:, 2:ny, :) - f(:, 1:ny-1, :)) / dy;

    % At the last column (j = ny), use backward difference
    %dfdy(:, ny, :) = (f(:, ny, :) - f(:, ny-1, :)) / dy;
    dfdy(:,ny,:)     = dfdy(:,ny-1,:);          % copy last interior gradient
end