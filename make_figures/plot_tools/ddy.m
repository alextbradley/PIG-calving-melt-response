function dMdy = ddy(M, dy);
% Return the derivative of M with respect to the second argument using second order centred finite differences in the interior and one sides finite differences at boundary points. 
dMdy = zeros(size(M));
dMdy(:, 2:end-1) = (M(:,3:end) - M(:,1:end-2))/ 2 /dy;
dMdy(:, 1)     = (-3/2 * M(:,1) + 2 * M(:, 2) - 1/2 * M(:, 3))/dy;
dMdy(:, end)   = (1/2 * M(:, end-2) - 2 * M(:, end-1) + 3/2 * M(:, end)) / dy;
end
