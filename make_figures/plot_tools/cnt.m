function [fig,c] = cnt(M)
fig = clf; contourf(M, 20, 'linestyle', 'none'); c = colorbar; axis equal
end
