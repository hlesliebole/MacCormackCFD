function d2fdy2 = d2dy2_central(f, dy)
% D2DY2_CENTRAL Computes the second derivative of f in y using central differences
% with periodic boundary conditions.
%
%   d2fdy2 = d2dy2_central(f, dy)
%
%   Input:
%       f  - input field [nx x ny]
%       dy - grid spacing in y-direction
%
%   Output:
%       d2fdy2 - second derivative ∂²f/∂y² [nx x ny]

    % Periodic central difference in y
    f_plus  = circshift(f, [0, -1]);  % f(i, j+1)
    f_minus = circshift(f, [0,  1]);  % f(i, j-1)

    d2fdy2 = (f_plus - 2*f + f_minus) / dy^2;
end