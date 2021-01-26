%> @file locMatAlg.m Assembly of a local assembly matrix under consideration of the global orientation.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%> @param  this  problem class [ Transport ]
%> @param  kT    number of element  [ scalar ]
%> @param  di  the inverse diffusion coefficient [ 2 x 2]
%>
%> @retval stiff  local assembly matrix @f$[\int_T D^-1 \vec{phi}_k \cdot\vec{phi}_j]_{j,k}@f$ [3 x 3]
%>

function stiff = locMatAlgVector(this, di)


g = this.grid;
s = g.sigE0T(:, :); % signi of local edges
J = g.A; %improved SG
arE = g.areaE(g.E0T(:, :));

stiff = zeros(g.numT, 3, 3);
stiff(:, 1, 1) = (arE(:, 1).^2 .* s(:, 1).^2 .* (2 * J(:, 1).^2 .* di(:, 1) + 2 * J(:, 3).^2 .* di(:, 1) + 2 * J(:, 2).^2 .* di(:, 4) + 2 * J(:, 4).^2 .* di(:, 4) + 2 * J(:, 1) .* J(:, 3) .* di(:, 1) + 2 * J(:, 1) .* J(:, 2) .* di(:, 3) + J(:, 1) .* J(:, 4) .* di(:, 3) + J(:, 3) .* J(:, 2) .* di(:, 3) + 2 * J(:, 3) .* J(:, 4) .* di(:, 3) + 2 * J(:, 1) .* J(:, 2) .* di(:, 2) + J(:, 1) .* J(:, 4) .* di(:, 2) + J(:, 3) .* J(:, 2) .* di(:, 2) + 2 .* J(:, 3) .* J(:, 4) .* di(:, 2) + 2 * J(:, 2) .* J(:, 4) .* di(:, 4))) ./ (12 * sqrt(2)^2) ./ (g.areaT(:) * 2);
stiff(:, 1, 2) = -(2^(1 / 2) * arE(:, 1) .* arE(:, 2) .* s(:, 1) .* s(:, 2) .* (2 * J(:, 1).^2 .* di(:, 1) - 2 * J(:, 3).^2 .* di(:, 1) + 2 * J(:, 2).^2 .* di(:, 4) - 2 .* J(:, 4).^2 .* di(:, 4) + 2 * J(:, 1) .* J(:, 3) .* di(:, 1) + 2 .* J(:, 1) .* J(:, 2) .* di(:, 3) + 3 * J(:, 1) .* J(:, 4) .* di(:, 3) - J(:, 3) .* J(:, 2) .* di(:, 3) - 2 * J(:, 3) .* J(:, 4) .* di(:, 3) + 2 * J(:, 1) .* J(:, 2) .* di(:, 2) - J(:, 1) .* J(:, 4) .* di(:, 2) + 3 * J(:, 3) .* J(:, 2) .* di(:, 2) - 2 * J(:, 3) .* J(:, 4) .* di(:, 2) + 2 * J(:, 2) .* J(:, 4) .* di(:, 4))) / (24 * sqrt(2) * 1) ./ (g.areaT(:) * 2);
stiff(:, 1, 3) = -(2 * 2^(1 / 2) * J(:, 3).^2 .* arE(:, 1) .* arE(:, 3) .* di(:, 1) .* s(:, 1) .* s(:, 3) - 2 * 2^(1 / 2) * J(:, 1).^2 .* arE(:, 1) .* arE(:, 3) .* di(:, 1) .* s(:, 1) .* s(:, 3) - 2 * 2^(1 / 2) * J(:, 2).^2 .* arE(:, 1) .* arE(:, 3) .* di(:, 4) .* s(:, 1) .* s(:, 3) + 2 * 2^(1 / 2) * J(:, 4).^2 .* arE(:, 1) .* arE(:, 3) .* di(:, 4) .* s(:, 1) .* s(:, 3) + 2 * 2^(1 / 2) * J(:, 1) .* J(:, 3) .* arE(:, 1) .* arE(:, 3) .* di(:, 1) .* s(:, 1) .* s(:, 3) - 2 * 2^(1 / 2) * J(:, 1) .* J(:, 2) .* arE(:, 1) .* arE(:, 3) .* di(:, 3) .* s(:, 1) .* s(:, 3) - 2^(1 / 2) * J(:, 1) .* J(:, 4) .* arE(:, 1) .* arE(:, 3) .* di(:, 3) .* s(:, 1) .* s(:, 3) + 3 * 2^(1 / 2) * J(:, 3) .* J(:, 2) .* arE(:, 1) .* arE(:, 3) .* di(:, 3) .* s(:, 1) .* s(:, 3) + 2 * 2^(1 / 2) * J(:, 3) .* J(:, 4) .* arE(:, 1) .* arE(:, 3) .* di(:, 3) .* s(:, 1) .* s(:, 3) - 2 * 2^(1 / 2) * J(:, 1) .* J(:, 2) .* arE(:, 1) .* arE(:, 3) .* di(:, 2) .* s(:, 1) .* s(:, 3) + 3 * 2^(1 / 2) * J(:, 1) .* J(:, 4) .* arE(:, 1) .* arE(:, 3) .* di(:, 2) .* s(:, 1) .* s(:, 3) - 2^(1 / 2) * J(:, 3) .* J(:, 2) .* arE(:, 1) .* arE(:, 3) .* di(:, 2) .* s(:, 1) .* s(:, 3) + 2 * 2^(1 / 2) * J(:, 3) .* J(:, 4) .* arE(:, 1) .* arE(:, 3) .* di(:, 2) .* s(:, 1) .* s(:, 3) + 2 * 2^(1 / 2) * J(:, 2) .* J(:, 4) .* arE(:, 1) .* arE(:, 3) .* di(:, 4) .* s(:, 1) .* s(:, 3)) / (24 * sqrt(2) * 1) ./ (g.areaT(:) * 2);

stiff(:, 2, 1) = -(2^(1 / 2) * arE(:, 1) .* arE(:, 2) .* s(:, 1) .* s(:, 2) .* (2 * J(:, 1).^2 .* di(:, 1) - 2 * J(:, 3).^2 .* di(:, 1) + 2 * J(:, 2).^2 .* di(:, 4) - 2 * J(:, 4).^2 .* di(:, 4) + 2 * J(:, 1) .* J(:, 3) .* di(:, 1) + 2 * J(:, 1) .* J(:, 2) .* di(:, 3) - J(:, 1) .* J(:, 4) .* di(:, 3) + 3 * J(:, 3) .* J(:, 2) .* di(:, 3) - 2 * J(:, 3) .* J(:, 4) .* di(:, 3) + 2 * J(:, 1) .* J(:, 2) .* di(:, 2) + 3 * J(:, 1) .* J(:, 4) .* di(:, 2) - J(:, 3) .* J(:, 2) .* di(:, 2) - 2 * J(:, 3) .* J(:, 4) .* di(:, 2) + 2 * J(:, 2) .* J(:, 4) .* di(:, 4))) / (24 * sqrt(2) * 1) ./ (g.areaT(:) * 2);
stiff(:, 2, 2) = (arE(:, 2).^2 .* s(:, 2).^2 .* (6 * J(:, 1).^2 .* di(:, 1) + 2 * J(:, 3).^2 .* di(:, 1) + 6 * J(:, 2).^2 .* di(:, 4) + 2 * J(:, 4).^2 .* di(:, 4) - 6 * J(:, 1) .* J(:, 3) .* di(:, 1) + 6 * J(:, 1) .* J(:, 2) .* di(:, 3) - 3 * J(:, 1) .* J(:, 4) .* di(:, 3) - 3 * J(:, 3) .* J(:, 2) .* di(:, 3) + 2 * J(:, 3) .* J(:, 4) .* di(:, 3) + 6 * J(:, 1) .* J(:, 2) .* di(:, 2) - 3 * J(:, 1) .* J(:, 4) .* di(:, 2) - 3 * J(:, 3) .* J(:, 2) .* di(:, 2) + 2 * J(:, 3) .* J(:, 4) .* di(:, 2) - 6 * J(:, 2) .* J(:, 4) .* di(:, 4))) / (24 * 1^2) ./ (g.areaT(:) * 2);
stiff(:, 2, 3) = -(arE(:, 2) .* arE(:, 3) .* s(:, 2) .* s(:, 3) .* (2 * J(:, 1).^2 .* di(:, 1) + 2 * J(:, 3).^2 .* di(:, 1) + 2 * J(:, 2).^2 .* di(:, 4) + 2 * J(:, 4).^2 .* di(:, 4) - 6 * J(:, 1) .* J(:, 3) .* di(:, 1) + 2 * J(:, 1) .* J(:, 2) .* di(:, 3) - J(:, 1) .* J(:, 4) .* di(:, 3) - 5 * J(:, 3) .* J(:, 2) .* di(:, 3) + 2 * J(:, 3) .* J(:, 4) .* di(:, 3) + 2 * J(:, 1) .* J(:, 2) .* di(:, 2) - 5 * J(:, 1) .* J(:, 4) .* di(:, 2) - J(:, 3) .* J(:, 2) .* di(:, 2) + 2 * J(:, 3) .* J(:, 4) .* di(:, 2) - 6 * J(:, 2) .* J(:, 4) .* di(:, 4))) / (24 * 1 * 1) ./ (g.areaT(:) * 2);

stiff(:, 3, 1) = -(2 .* 2.^(1 / 2) .* J(:, 3).^2 .* arE(:, 1) .* arE(:, 3) .* di(:, 1) .* s(:, 1) .* s(:, 3) - 2 .* 2.^(1 / 2) .* J(:, 1).^2 .* arE(:, 1) .* arE(:, 3) .* di(:, 1) .* s(:, 1) .* s(:, 3) - 2 .* 2.^(1 / 2) .* J(:, 2).^2 .* arE(:, 1) .* arE(:, 3) .* di(:, 4) .* s(:, 1) .* s(:, 3) + 2 .* 2.^(1 / 2) .* J(:, 4).^2 .* arE(:, 1) .* arE(:, 3) .* di(:, 4) .* s(:, 1) .* s(:, 3) + 2 .* 2.^(1 / 2) .* J(:, 1) .* J(:, 3) .* arE(:, 1) .* arE(:, 3) .* di(:, 1) .* s(:, 1) .* s(:, 3) - 2 .* 2.^(1 / 2) .* J(:, 1) .* J(:, 2) .* arE(:, 1) .* arE(:, 3) .* di(:, 3) .* s(:, 1) .* s(:, 3) + 3 .* 2.^(1 / 2) .* J(:, 1) .* J(:, 4) .* arE(:, 1) .* arE(:, 3) .* di(:, 3) .* s(:, 1) .* s(:, 3) - 2.^(1 / 2) .* J(:, 3) .* J(:, 2) .* arE(:, 1) .* arE(:, 3) .* di(:, 3) .* s(:, 1) .* s(:, 3) + 2 .* 2.^(1 / 2) .* J(:, 3) .* J(:, 4) .* arE(:, 1) .* arE(:, 3) .* di(:, 3) .* s(:, 1) .* s(:, 3) - 2 .* 2.^(1 / 2) .* J(:, 1) .* J(:, 2) .* arE(:, 1) .* arE(:, 3) .* di(:, 2) .* s(:, 1) .* s(:, 3) - 2.^(1 / 2) .* J(:, 1) .* J(:, 4) .* arE(:, 1) .* arE(:, 3) .* di(:, 2) .* s(:, 1) .* s(:, 3) + 3 .* 2.^(1 / 2) .* J(:, 3) .* J(:, 2) .* arE(:, 1) .* arE(:, 3) .* di(:, 2) .* s(:, 1) .* s(:, 3) + 2 .* 2.^(1 / 2) .* J(:, 3) .* J(:, 4) .* arE(:, 1) .* arE(:, 3) .* di(:, 2) .* s(:, 1) .* s(:, 3) + 2 .* 2.^(1 / 2) .* J(:, 2) .* J(:, 4) .* arE(:, 1) .* arE(:, 3) .* di(:, 4) .* s(:, 1) .* s(:, 3)) / (24 .* sqrt(2) .* 1) ./ (g.areaT(:) * 2);
stiff(:, 3, 2) = -(arE(:, 2) .* arE(:, 3) .* s(:, 2) .* s(:, 3) .* (2 * J(:, 1).^2 .* di(:, 1) + 2 * J(:, 3).^2 .* di(:, 1) + 2 * J(:, 2).^2 .* di(:, 4) + 2 * J(:, 4).^2 .* di(:, 4) - 6 * J(:, 1) .* J(:, 3) .* di(:, 1) + 2 * J(:, 1) .* J(:, 2) .* di(:, 3) - 5 * J(:, 1) .* J(:, 4) .* di(:, 3) - J(:, 3) .* J(:, 2) .* di(:, 3) + 2 * J(:, 3) .* J(:, 4) .* di(:, 3) + 2 * J(:, 1) .* J(:, 2) .* di(:, 2) - J(:, 1) .* J(:, 4) .* di(:, 2) - 5 * J(:, 3) .* J(:, 2) .* di(:, 2) + 2 * J(:, 3) .* J(:, 4) .* di(:, 2) - 6 * J(:, 2) .* J(:, 4) .* di(:, 4))) / (24 * 1 * 1) ./ (g.areaT(:) * 2);
stiff(:, 3, 3) = (arE(:, 3).^2 .* s(:, 3).^2 .* (2 .* J(:, 1).^2 .* di(:, 1) + 6 .* J(:, 3).^2 .* di(:, 1) + 2 .* J(:, 2).^2 .* di(:, 4) + 6 .* J(:, 4).^2 .* di(:, 4) - 6 .* J(:, 1) .* J(:, 3) .* di(:, 1) + 2 .* J(:, 1) .* J(:, 2) .* di(:, 3) - 3 .* J(:, 1) .* J(:, 4) .* di(:, 3) - 3 .* J(:, 3) .* J(:, 2) .* di(:, 3) + 6 .* J(:, 3) .* J(:, 4) .* di(:, 3) + 2 .* J(:, 1) .* J(:, 2) .* di(:, 2) - 3 .* J(:, 1) .* J(:, 4) .* di(:, 2) - 3 .* J(:, 3) .* J(:, 2) .* di(:, 2) + 6 .* J(:, 3) .* J(:, 4) .* di(:, 2) - 6 .* J(:, 2) .* J(:, 4) .* di(:, 4))) / (24 .* 1.^2) ./ (g.areaT(:) * 2);

end