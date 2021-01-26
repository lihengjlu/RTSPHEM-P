%> @file testDNPP_force.m Additional force term @f$\vec{f}@f$ for the Stokes problem.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%>
%> @param t evaluation time [scalar]
%> @param x sample point in domain @f$[2\times 1]@f$
%> @retval ret see function description @f$[2\times 1]@f$

function ret = testDNPP_force(t, x)
msg = 'HyPHM: RTFM.';

assert(isequal(size(x), [2, 1]), msg)
assert(isscalar(t), msg)

cx = cos(pi*x(1));
cy = cos(pi*x(2));
sx = sin(pi*x(1));
sy = sin(pi*x(2));

ret = [t * (-cx * sy) + pi / 2 * sin(2 * pi * x(1)) + t^2 / pi * (sy - cx) * sx; ...
    t * (sx * cy) + pi / 2 * sin(2 * pi * x(2)) + t^2 / pi * (sy - cx) * cy];

end
