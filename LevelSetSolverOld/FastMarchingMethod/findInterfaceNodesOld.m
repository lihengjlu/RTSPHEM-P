function [isInitial] = findInterfaceNodesOld(grid, lsf)
% Finds the nodes of a cartesian grid that are adjacent to the interface LSF=0.

assert(isa(grid, 'CartesianGrid'), ...
    'Argument not of class ''CartesianGrid''.');

numNodes = grid.nodes;
isInitial = false(numNodes, 1);
val = NaN(3, numNodes);
val(1, :) = lsf(:)';

for d = 1:grid.dimension

    forwardIndex = grid.getForward(1:numNodes, d);
    val(2, ~isnan(forwardIndex)) = ...
        lsf(forwardIndex(~isnan(forwardIndex)));

    backwardIndex = grid.getBackward(1:numNodes, d);
    val(3, ~isnan(backwardIndex)) = ...
        lsf(backwardIndex(~isnan(backwardIndex)));

    isInitial(min(val).*max(val) <= 0) = true;
end


end
