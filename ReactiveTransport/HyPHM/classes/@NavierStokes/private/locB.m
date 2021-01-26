% locB performs the assembly of a local assembly matrix
% under consideration of the global orientation.
%
% Input
%   d        [ NavierStokes ]  data required by solver
%   kT       [ scalar ]        number of current element
%
% Output
%   stiff    [6 x 6]             local assembly matrix
%
%   stiff(N', N) = <dxphi_N, dxphi_N'> ,  N, N' = 1,...,6.
%
%
% See also NavierStokes.computeLevel
%
% Copyright 2009, 2010 F. Frank, Chair of Applied Mathematics I,
%           Department of Mathematics, University of Erlangen-Nuremberg,
%           G E R M A N Y

function stiff = locB(d, kT)

B1 = ...
    [1 / 2, 1 / 6, 0, 0, 0, -2 / 3; ...
    1 / 6, 1 / 2, 0, 0, 0, -2 / 3; ...
    0, 0, 0, 0, 0, 0; ...
    0, 0, 0, 4 / 3, -4 / 3, 0; ...
    0, 0, 0, -4 / 3, 4 / 3, 0; ...
    -2 / 3, -2 / 3, 0, 0, 0, 4 / 3];

B2 = ...
    [1 / 2, 0, 1 / 6, 0, -2 / 3, 0; ...
    0, 0, 0, 0, 0, 0; ...
    1 / 6, 0, 1 / 2, 0, -2 / 3, 0; ...
    0, 0, 0, 4 / 3, 0, -4 / 3; ...
    -2 / 3, 0, -2 / 3, 0, 4 / 3, 0; ...
    0, 0, 0, -4 / 3, 0, 4 / 3];

B3 = ...
    [1, 1 / 6, 1 / 6, 0, -2 / 3, -2 / 3; ...
    1 / 6, 0, -1 / 6, 2 / 3, 0, -2 / 3; ...
    1 / 6, -1 / 6, 0, 2 / 3, -2 / 3, 0; ...
    0, 2 / 3, 2 / 3, 4 / 3, -4 / 3, -4 / 3; ...
    -2 / 3, 0, -2 / 3, -4 / 3, 4 / 3, 4 / 3; ...
    -2 / 3, -2 / 3, 0, -4 / 3, 4 / 3, 4 / 3];


A = squeeze(d.grid.A(kT, :, :));

stiff = (A(2, 2)^2 * B1 - A(2, 1) * A(2, 2) * B3 + A(2, 1)^2 * B2) / abs(det(A));

end

% Remark: d.isEvolution is the flag which decides whether the evolution
% term should be taken into account.