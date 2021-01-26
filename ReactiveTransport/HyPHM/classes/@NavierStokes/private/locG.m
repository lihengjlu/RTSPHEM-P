% locG performs the assembly of a local assembly matrix
% under consideration of the global orientation.
%
% Input
%   d        [ NavierStokes ]  data required by solver
%   kT       [ scalar ]        number of current element
%
% Output
%   stiff    [3 x 6]             local assembly matrix
%
%   stiff(N', N) = <lambda_N', dyphi_N> ,   N' = 1,...,3, N = 1,...,6.
%
%
% See also NavierStokes.computeLevel
%
% Copyright 2009, 2010 F. Frank, Chair of Applied Mathematics I,
%           Department of Mathematics, University of Erlangen-Nuremberg,
%           G E R M A N Y

function stiff = locG(d, kT)

G1 = [-1 / 6, 0, 0, 1 / 6, -1 / 6, 1 / 6; ...
    0, 1 / 6, 0, 1 / 6, -1 / 6, -1 / 6; ...
    0, 0, 0, 1 / 3, -1 / 3, 0];

G2 = [-1 / 6, 0, 0, 1 / 6, 1 / 6, -1 / 6; ...
    0, 0, 0, 1 / 3, 0, -1 / 3; ...
    0, 0, 1 / 6, 1 / 6, -1 / 6, -1 / 6];

A = squeeze(d.grid.A(kT, :, :));

stiff = -A(1, 2) * G1 + A(1, 1) * G2;

end

% Remark: d.isEvolution is the flag which decides whether the evolution
% term should be taken into account.