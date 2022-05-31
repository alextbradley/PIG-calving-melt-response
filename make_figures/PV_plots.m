% Produce plots relating to potential voriticity for the default case:
% (i)   u.grad[(f + zeta)/H)
% (ii)  u.grad(f/h)
% (iii) u.grad(zeta/h)

% NB: Many of the data files referred to in this script are too large to be hosted online. These files are hostedinternally as BAS.
% Please email Alex Bradley (aleey@bas.ac.uk) to obtain a copy.
%Alex Bradley (aleey@bas.ac.uk) 27/05/2021. MIT license.

addpath('plot_tools');
%
% Flags
%
gendata = 1; %specify whether to pass through the generate data loop
save_flag = 0;
free_surf = 1; %apply the free surface correction to velocities and water column thickness

% Data locations
rootdir = '/data/oceans_output/shelf/aleey/mitgcm/APIGi_'; %not in git repo
topodir = '../gendata/topo_files/';
bathy_path = '../gendata/bathy_files/bathymetry_H400.shice';


%grid details
nx=120; % number of grid cells along longitudinal direction
ny=320; % number of grid cells along latitudinal direction
nz=110; % number of vertical grid cells
dx=400;
dy=400;
dz=10;
X = 0:dx:(nx-1)*dx;
Y = 0:dx:(ny-1)*dy;
Z = 0:dz:(nz-1)*dz;
[XX,YY] = meshgrid(X,Y);
YYt = YY';
idx = (YYt < 30e3); %inner cavity definition

%parameters
secs_per_year = 365.25*24*60*60;
density_ice = 918.0;

%time details
ntout1 = 6;
ntout2 = 12; %define time period to average over

%
% Generate data loop
%
%uncalved
run_nos = "077";
extent = 84;
W = 100; %ridge gap
numplot = 1;

%ridge calved
%run_nos = "084";
%extent = 50;
%W = 100;
%numplot = 2; 

H = 400; %ridge height (always 400);
sz = length(run_nos);
if gendata

%load bathy
fid = fopen(bathy_path);
bathy = fread(fid, 'real*8', 'b');
bathy = reshape(bathy, [nx, ny]);
%bathy(bathy == 0) = nan;

%draft
topo_fname=  ['shelfice_topo_H' num2str(H) '_W' num2str(W) '_extent' num2str(extent) 'km.bin'];
topo_fid = fopen(strcat(topodir, '/',topo_fname));
topo = fread(topo_fid, 'real*8', 'b');
topo = reshape(topo, [nx, ny]);
fname =  strcat(rootdir, run_nos, '/run/state2D.nc');
etan = ncread(fname, 'ETAN', [1,1,ntout1], [Inf, Inf, 1 + ntout2 - ntout1]);
etan = mean(etan, 3);

%water column thickness
if free_surf
H = topo - bathy + etan;
else
H = topo - bathy;
end

%Velocities
UVEL_fname = strcat(rootdir, run_nos, '/run/stateUvel.nc');
UVEL = ncread(UVEL_fname, 'UVEL', [1,1,1,ntout1], [Inf, Inf, Inf,  1+ntout2 - ntout1]);
UVEL = mean(UVEL, 4);
VVEL_fname = strcat(rootdir, run_nos, '/run/stateVvel.nc');
VVEL = ncread(VVEL_fname, 'VVEL', [1,1,1,ntout1], [Inf, Inf, Inf,  1+ntout2 - ntout1]);
VVEL = mean(VVEL, 4);

%% Processed quantities 
% BSF
vvel = squeeze(sum(VVEL, 3)) * dz; %units m^2 /s
stream=zeros(size(vvel));
stream(nx,:)=vvel(nx,:)*dx;
for p=nx-1:-1:1
 stream(p,:)=stream(p+1,:) + vvel(p,:)*dx;
end
stream = stream/1e6; %convert to sv
streamsm = smooth2a(stream, 2,2);


% Barotropic velocities
vvel = squeeze(sum(VVEL, 3)) * dz;
uvel = squeeze(sum(UVEL, 3)) * dz;
if free_surf
%add the etan correction to the top cell
for p = 1:nx
for q = 1:ny
idx_top = find(squeeze(VVEL(p,q,:)) ~= 0, 1, 'First'); %find index of top cell
if ~isempty(idx_top)
vvel(p,q) = vvel(p,q) + VVEL(p,q,idx_top)*etan(p,q); %add the free surface corrections
uvel(p,q) = uvel(p,q) + UVEL(p,q,idx_top)*etan(p,q); 
end 
end %end loop over q = 1:ny
end %end loop over p = 1:nx
end
vvel = vvel ./ H;
uvel = uvel ./ H; 

%Put velocities onto cell centres
uvelC = zeros(nx,ny);
uvelC(1:end-1,:) = (uvel(1:end-1,:) + uvel(2:end,:))/2;
uvelC(end,:)     = (-uvel(end-1,:) + 3*uvel(end,:))/2;

vvelC = zeros(nx,ny);
vvelC(:,1:end-1) = (vvel(:,1:end-1) + vvel(:,2:end))/2;
vvelC(:,end)     = (-vvel(:,end-1) + 3*vvel(:,end))/ 2; 

%interpolate velocities onto cell corners from (default) cell edges 
%uvelC = zeros(nx,ny);
%uvelC(:, 2:end) = (uvel(:,1:end-1) + uvel(:,2:end))/2;
%uvelC(:, 1) = (3*uvel(:,1) - uvel(:,2))/2;

%vvelC = zeros(nx,ny);
%vvelC(2:end, :) = (vvel(1:end-1,:) + vvel(2:end,:))/2;
%vvelC(1,:) = (3*vvel(1,:) - vvel(2,:))/2;

vvel = vvelC;
uvel = uvelC;

%velocity gradients
dvdx = ddx(vvel, dx); 
dudy = ddy(uvel, dy); 

%smooth dudy
dudyS = dudy;
for i = 1:nx; dudyS(i,:) = smooth(dudyS(i,:)); end
for j = 1:ny; dudyS(:,j) = smooth(dudyS(:,j)); end
dudy = dudyS;


% PV 
zeta = dvdx - dudy;
f = 2*7.2921*1e-5*sind(-75);
BPV = (f + zeta) ./ H;

% Barotropic vorticity gradients
dBPVdx = ddx(BPV,dx);
dBPVdy = ddy(BPV,dy);
u_dot_grad_BPV = uvel .* dBPVdx + vvel .* dBPVdy;

% u . grad(f/H)
fOnH = f ./H;
dfOnHdx = ddx(fOnH, dx);
dfOnHdy = ddy(fOnH, dy);
u_dot_grad_fOnH = uvel .* dfOnHdx + vvel .* dfOnHdy;

% u. grad(zeta/H)
zetaOnH = zeta ./ H;
dzetaOnHdx = ddx(zetaOnH, dx);
dzetaOnHdy = ddy(zetaOnH, dy);
u_dot_grad_zetaOnH = uvel .* dzetaOnHdx + vvel .* dzetaOnHdy;

check = u_dot_grad_fOnH + u_dot_grad_zetaOnH - u_dot_grad_BPV; %should equal 0 

end

if 1
figure(numplot); clf;
% u velocity  
uvelS = saturate(uvel, 0.165, -0.15);
ax(1) = subplot(2,4,1); contourf(X/1e3, Y/1e3, uvelS', 50, 'linestyle', 'none'); colorbar; title('u (m/s)');xlabel('x (km)'); ylabel('y (km)');
colormap(ax(1), redblue);

% v velocity
vvelS = saturate(vvel, 0.525, -0.5);
ax(3) = subplot(2,4,3); contourf(X/1e3, Y/1e3, vvelS', 50, 'linestyle', 'none'); colorbar; title('v (m/s)');xlabel('x (km)'); ylabel('y (km)');
colormap(ax(3), redblue);

%dudy
dudyS = saturate(dudy, 3.15*1e-5, -3*1e-5); 
ax(2) = subplot(2,4,2); contourf(X/1e3, Y/1e3, dudyS', 50, 'linestyle', 'none'); colorbar; title('du/dy (1/s)');xlabel('x (km)'); ylabel('y (km)');
colormap(ax(2), redblue);

%dvdx
dvdxS = saturate(dvdx, 2.6*1e-4, -2.5*1e-4);
ax(4) = subplot(2,4,4); contourf(X/1e3, Y/1e3, dvdxS', 50, 'linestyle', 'none'); colorbar; title('dv/dx (1/s)');xlabel('x (km)'); ylabel('y (km)');
colormap(ax(4), redblue);

%u . grad(BPV)
satval = 0.5*1e-11;
u_dot_grad_BPVS = saturate(u_dot_grad_BPV, 1.05*satval, -satval);
ax(5) = subplot(2,4,5); contourf(X/1e3, Y/1e3,1e11 * u_dot_grad_BPVS', 50, 'linestyle', 'none'); colorbar; title('$10^{11} \times u. \nabla [(f + \zeta)/ H]$', 'interpreter', 'latex');xlabel('x (km)'); ylabel('y (km)');
colormap(ax(5), redblue);

% u . grad(f/ H)
u_dot_grad_fOnHS = saturate(u_dot_grad_fOnH, 1.05*satval, -satval);
ax(6) = subplot(2,4,6); contourf(X/1e3, Y/1e3, 1e11 * u_dot_grad_fOnHS', 50, 'linestyle', 'none'); colorbar; title('$10^{11} \times u. \nabla (f / H)$', 'interpreter', 'latex');xlabel('x (km)'); ylabel('y (km)');
colormap(ax(6), redblue);

% u . grad(zeta / H)
u_dot_grad_zetaOnHS = saturate(u_dot_grad_zetaOnH, 1.05*satval, -satval);
ax(7) = subplot(2,4,7); contourf(X/1e3, Y/1e3, 1e11 * u_dot_grad_zetaOnHS', 50, 'linestyle', 'none'); colorbar; title('$10^{11} \times u. \nabla (\zeta / H)$', 'interpreter', 'latex');xlabel('x (km)'); ylabel('y (km)');
colormap(ax(7), redblue);

zetaS = saturate(zeta, 2*1.08*1e-4, -2*1e-4);
ax(8) = subplot(2,4,8); contourf(X/1e3, Y/1e3, zetaS', 20, 'linestyle', 'none'); colorbar;xlabel('x (km)'); ylabel('y (km)'); colorbar; title('\zeta')
colormap(ax(8), redblue);
clear ax
end
%% repeat final four figures in plot of own
fig = figure(numplot + 1); clf; 
fig.Position(3:4) = [1200, 420];
ax(1) = subplot(1,5,1); contourf(X/1e3, Y/1e3,1e11 * u_dot_grad_BPVS', 50, 'linestyle', 'none'); colorbar; title('$10^{11} \times u. \nabla [(f + \zeta)/ H]$', 'interpreter', 'latex');xlabel('x (km)'); ylabel('y (km)');
colormap(ax(1), redblue);
axnew(1) = axes(); axnew(1).Position = ax(1).Position;
[c,h] = contour(axnew(1), X/1e3, Y/1e3,streamsm', [-0.7, -0.5, -0.3, -0.1,0], 'k');
clabel(c,h); axnew(1).Visible = 'off';

% u . grad(f/ H)
ax(2) = subplot(1,5,2); contourf(X/1e3, Y/1e3, 1e11 * u_dot_grad_fOnHS', 50, 'linestyle', 'none'); colorbar; title('$10^{11} \times u. \nabla (f / H)$', 'interpreter', 'latex');xlabel('x (km)'); ylabel('y (km)');
colormap(ax(2), redblue);
axnew(2) = axes(); axnew(2).Position = ax(2).Position;
[c,h] = contour(axnew(2), X/1e3, Y/1e3,streamsm', [-0.7, -0.5, -0.3, -0.1,0], 'k');
clabel(c,h); 
axnew(2).Visible = 'off';

% u . grad(zeta / H)
ax(3) = subplot(1,5,3); contourf(X/1e3, Y/1e3, 1e11 * u_dot_grad_zetaOnHS', 50, 'linestyle', 'none'); colorbar; title('$10^{11} \times u. \nabla (\zeta / H)$', 'interpreter', 'latex');xlabel('x (km)'); ylabel('y (km)');
colormap(ax(3), redblue);
axnew(3) = axes(); axnew(3).Position = ax(3).Position;
[c,h] = contour(axnew(3), X/1e3, Y/1e3,streamsm', [-0.7, -0.5, -0.3, -0.1,0], 'k');
clabel(c,h); 
axnew(3).Visible = 'off';

ax(4) = subplot(1,5,4); contourf(X/1e3, Y/1e3, zetaS', 20, 'linestyle', 'none'); colorbar;xlabel('x (km)'); ylabel('y (km)'); colorbar; title('\zeta')
colormap(ax(4), redblue);
axnew(4) = axes(); axnew(4).Position = ax(4).Position;
[c,h] = contour(axnew(4), X/1e3, Y/1e3,streamsm', [-0.7, -0.5, -0.3, -0.1,0], 'k');
clabel(c,h); 
axnew(4).Visible = 'off';

ax(5) =subplot(1,5,5); contourf(X/1e3, Y/1e3, BPV', 20, 'linestyle', 'none'); colorbar;xlabel('x (km)'); ylabel('y (km)'); colorbar; title('$(f + \zeta)/H$', 'interpreter', 'latex')
colormap(ax(5), parula);
axnew(5) = axes(); axnew(5).Position = ax(5).Position;
[c,h] = contour(axnew(5), X/1e3, Y/1e3,streamsm', [-0.7, -0.5, -0.3, -0.1,0], 'k');
clabel(c,h); 
axnew(5).Visible = 'off';

for i = 1:5
axnew(i).Position = ax(i).Position;
%axis equal
end
fig = gcf; fig.Position(3:4) = [1200, 311]; %ensures axis equal
