% Make figure 6 in the updated ms: pv and shear for the ridge crest calving case
% (a) Plot of  water column thickness with barotropic stream function overlain
% (b) Plot of barotropic potential vorticity with barotropic stream function overlain
% (c) Zonal cross section of meridional velocity at the ridge crest
% 
% NB: The data files referred to in this script are too large to be hosted online. These files are hosted internally as BAS.
% Please email Alex Bradley (aleey@bas.ac.uk) to obtain a copy.
%Alex Bradley (aleey@bas.ac.uk) 06/07/2022. MIT license.

%
% Flags
%
gendata = 1; %specify whether to pass through the generate data loop
save_flag = 0;

%
% Preliminaries
%
addpath("plot_tools");
plot_defaults
label_size = 11;
ax_fontsize = 10;
ridge_idx = 125; %index of line for zonal cross section

% Data locations
rootdir = '/data/oceans_output/shelf/aleey/mitgcm/APIGi/APIGi_'; %not in git repo
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
lambda = 7.61*1e-4;%constants in liquidus
gamma = 5.73*1e-2;
T0 = 8.32*1e-2;

%time details
ntout1 = 6;
ntout2 = 12; %define time period to average over


%
% Generate data loop
%
run_nos = "084";
sz = length(run_nos);
extent = 50;
H = 400; %ridge height (always 400);
W = 100; %ridge gap

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


%Velocities
UVEL_fname = strcat(rootdir, run_nos, '/run/stateUvel.nc');
UVEL = ncread(UVEL_fname, 'UVEL', [1,1,1,ntout1], [Inf, Inf, Inf,  1+ntout2 - ntout1]);
UVEL = mean(UVEL, 4);
VVEL_fname = strcat(rootdir, run_nos, '/run/stateVvel.nc');
VVEL = ncread(VVEL_fname, 'VVEL', [1,1,1,ntout1], [Inf, Inf, Inf,  1+ntout2 - ntout1]);
VVEL = mean(VVEL, 4);


%compute BSF
vvel = squeeze(sum(VVEL, 3)) * dz; %units m^2 /s
stream=zeros(size(vvel));
stream(nx,:)=vvel(nx,:)*dx;
for p=nx-1:-1:1
 stream(p,:)=stream(p+1,:) + vvel(p,:)*dx;
end
stream = stream/1e6; %convert to sv
streamsm = smooth2a(stream, 2,2);

% Barotropic velocities
H = topo-bathy;
vvel = squeeze(sum(VVEL, 3)) * dz;
uvel = squeeze(sum(UVEL, 3)) * dz;
vvel = vvel ./ H;
uvel = uvel ./ H;

%Put velocities onto cell centres
uvelC = zeros(nx,ny);
uvelC(1:end-1,:) = (uvel(1:end-1,:) + uvel(2:end,:))/2;
uvelC(end,:)     = (-uvel(end-1,:) + 3*uvel(end,:))/2;
vvelC = zeros(nx,ny);
vvelC(:,1:end-1) = (vvel(:,1:end-1) + vvel(:,2:end))/2;
vvelC(:,end)     = (-vvel(:,end-1) + 3*vvel(:,end))/ 2;
vvel = vvelC;
uvel = uvelC;

%velocity gradients and smooth
dvdx = ddx(vvel, dx);
dudy = ddy(uvel, dy);
dudyS = dudy;
for i = 1:nx; dudyS(i,:) = smooth(dudyS(i,:)); end
for j = 1:ny; dudyS(:,j) = smooth(dudyS(:,j)); end
dudy = dudyS;

% PV quantities
zeta = dvdx - dudy;
f = 2*7.2921*1e-5*sind(-75);
BPV = (f + zeta) ./ H;
end %end gendata loop

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plots %%%%%%%%%%%%%%%%%%%%%%%
% generate plot positions, starting with top four
positions = [0.08, 0.17, 0.55, 0.38;
	     0.08, 0.6, 0.55, 0.38;
	     0.73, 0.2, 0.18, 0.7];
cmap = lighter_blue_parula(100,0.2);
colormap(cmap);

fig = figure(1);clf; fig.Position(3:4) = [1000,250];

% plot 1:  1/h and BSF contours
ax(1) = subplot('Position', positions(2,:)); hold on; box on;
ax(1).FontSize =10;
column_thickness = topo - bathy;
invcolthick = 1e3*1./column_thickness;
invcolthick = saturate(invcolthick, 10,1);
contourf((max(Y) - Y)/1e3,X/1e3,invcolthick, 50, 'linestyle', 'none');
cbar(1) = colorbar;
cbar(1).Label.Interpreter = 'latex';
cbar(1).Label.String = '$1/h~(10^{-3}~\mathrm{m}^{-1})$';
cbar(1).Label.FontSize = 10;
ax(1).YTick = [0,24,48];
ax(1).XLim = [0,128.1];
ax(1).YLim = [-0.5, 48.1];
ax(1).XTick = [];
ylabel('$x$~(km)', 'Interpreter', 'latex', 'FontSize' ,12);
txt(1) =text(ax(1),-20,46, '(a)', 'Interpreter', 'latex', 'FontSize',  12);
%add bsf contours
axnew = axes;
axnew.Position = ax(1).Position;
[C,h] =contour((max(Y) - Y)/1e3,X/1e3, streamsm, [-0.7, -0.5, -0.3, -0.1], 'k');
clabel(C,h);
xticks([]);
yticks([]);
set(axnew, 'color', 'none')
axnew.XLim = ax(1).XLim;
axnew.YLim = ax(1).YLim;


% plot 2:  f + zeta / H and BSF contours
ax(2) = subplot('Position', positions(1,:)); hold on; box on;
ax(2).FontSize =10;
bpvsat = 1e6* BPV; bpvsat = saturate(bpvsat,.5,-2);
contourf((max(Y) - Y)/1e3, X/1e3,bpvsat, 50, 'linestyle', 'none');
cbar(2) = colorbar;
%cbar(1).Location = 'northoutside';
%cbar(1).Position(end) = 0.01;
%cbar(1).Position(2) = 0.92;
cbar(2).Label.String = '$10^{6}\times[(f + \zeta )/ h]$';
cbar(2).Label.Interpreter = 'latex';
cbar(2).Label.FontSize = 10;
ax(2).YTick = [0,24,48];
ax(2).XLim = [0,128.1]; %make sure box appears
ax(2).YLim = [-0.5, 48.1];
%ax(2).XTick = [];
ylabel('$x$~(km)', 'Interpreter', 'latex', 'FontSize' ,12);
xlabel('$y$~(km)', 'Interpreter', 'latex', 'FontSize' ,12);
txt(2) =text(ax(2),-20,46, '(b)', 'Interpreter', 'latex', 'FontSize',  12);

%add bsf contours
axnew = axes;
axnew.Position = ax(2).Position;
[C,h] =contour((max(Y) - Y)/1e3,X/1e3, streamsm, [-0.7, -0.5, -0.3, -0.1], 'k');
clabel(C,h);
xticks([]);
yticks([]);
set(axnew, 'color', 'none')
axnew.XLim = ax(1).XLim;
axnew.YLim = ax(1).YLim;

% plot 3:  zonal cross section of meridional velocity at the ridge crest
ax(3) = subplot('Position', positions(3,:)); hold on; box on;
ax(3).FontSize =10;
UZS = squeeze(UVEL(:,ridge_idx,:));
for p = 1:nx
for q = 1:nz
        if bathy(p,ridge_idx) > -Z(q)
        UZS(p,q) = nan;
        end
        if topo(p,ridge_idx)< -Z(q)
        UZS(p,q) = nan;
        end
end
end
UZS = saturate(UZS,0.22, -0.05);
contourf(X/1e3,-Z, UZS',60, 'linestyle', 'none');
hold on
ylabel('depth (m)', 'Interpreter', 'latex', 'FontSize', 12);
xlabel('$x$~(km)', 'Interpreter', 'latex', 'FontSize' ,12);

cbar(3) = colorbar;
cbar(3).Label.String = '$u~(\mathrm{m}~\mathrm{s}^{-1}$ )';
cbar(3).Label.Interpreter = 'latex';
cbar(3).Label.FontSize = 10;
cbar(3).Position(1) = 0.915;
cbar(3).Position(2) = 0.206;
ax(3).YLim = [-701,-600];
ax(3).YTick = -700:50:-600;
ax(3).XLim = [0, 48];
ax(3).XTick = [0,24,48];
txt(3) =text(ax(3),-18,-591, '(c)', 'Interpreter', 'latex', 'FontSize',  12);
ax(3).Position = positions(3,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if save_flag
saveas(gcf, '../../PIG-calving-melt-response-tex/figures/figureF_ridge_crest_calved_characteristics.png');
end
