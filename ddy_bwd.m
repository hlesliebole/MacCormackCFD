function dfdy = ddy_bwd(f, dy, ~)
    %dfdy = (f - circshift(f, [0, 1, 0])) / dy;
    [nx, ny, nv] = size(f);
    dfdy = zeros(nx, ny, nv);

    % Backward difference in y: j = 2 to ny
    dfdy(:, 2:ny, :) = (f(:, 2:ny, :) - f(:, 1:ny-1, :)) / dy;

    % At the first column (j = 1), use forward difference
    %dfdy(:, 1, :) = (f(:, 2, :) - f(:, 1, :)) / dy;
    dfdy(:,1,:)    = dfdy(:,2,:);               % copy first interior gradient
end