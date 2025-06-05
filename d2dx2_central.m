function d2fdx2 = d2dx2_central(f, dx)
% D2DX2_CENTRAL Computes the second derivative of f in x using central differences
% with periodic boundary conditions.
%
%   d2fdx2 = d2dx2_central(f, dx)
%
%   Input:
%       f  - input field [nx x ny]
%       dx - grid spacing in x-direction
%
%   Output:
%       d2fdx2 - second derivative ∂²f/∂x² [nx x ny]

    % Periodic central difference in x
    f_plus  = circshift(f, [-1, 0]);  % f(i+1, j)
    f_minus = circshift(f, [ 1, 0]);  % f(i-1, j)

    d2fdx2 = (f_plus - 2*f + f_minus) / dx^2;
end