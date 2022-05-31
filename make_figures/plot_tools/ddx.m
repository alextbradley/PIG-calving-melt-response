function dMdx = ddx(M, dx);
% Return the derivative of M with respect to the first argument using second order centred finite differences in the interior and one sides finite differences at boundary points. 
dMdx = zeros(size(M));
dMdx(2:end-1, :)    = (M(3:end, :) - M(1:end-2,:)) / 2 / dx;
dMdx(1,:) = (-3/2 * M(1,:) + 2 * M(2,:) - 1/2 * M(3,:)) / 2 /dx;
dMdx(end,:) = (1/2 * M(end-2,:) - 2 * M(end-1,:) + 3/2 * M(end,:)) / 2 /dx;
end
