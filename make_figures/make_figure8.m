%Make figure 8 in the PIG calving melt rate manuscript: results from the W = 200, P600 runs.
%Column 1: Melt rate and BSF contours 
%Column 2: Zonal Sections velocity
%Column 3: Zonal Sections temperature
%Column 4: Zonal Section of meridional velocity
%
% NB: Many of the data files referred to in this script are too large to be hosted online. These files are hosted internally as BAS.
% Please email Alex Bradley (aleey@bas.ac.uk) to obtain a copy.
%Alex Bradley (aleey@bas.ac.uk) 27/05/2021. MIT license.

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
figure(1); clf;
fig = gcf; fig.Position(3:4) = [1000, 800];

m% Data locations
rootdir = '/data/oceans_output/shelf/aleey/mitgcm/APIGi/APIGi_'; %nb: not in git repo
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
ntout1 = 11;
ntout2 = 12; %define time period to average over

%cross section info
zonal_idx = 50; %20km downstream
%zonal_idx = 76; %30 km downstream of gl (entrance to "inner cavity")
%zonal_idx = 126; %index of top of bump
merid_idx = floor(nx/2); %index of midpt in x

%
% Generate data loop
%
run_nos = ["125", "126", "127", "128", "129", "130", "131", "132", "133", "134"];
sz = length(run_nos);
extent = [84,80,75,70,65,60,55,50,45,40];
H = 400; %ridge height (always 400);
W = 100; %ridge gap

%generate data loop
if gendata

%setup storage
melt_scenarios = cell(1,sz);
topo_scenarios = cell(1,sz);
zonal_theta_scenarios = cell(1,sz);
zonal_salt_scenarios = cell(1,sz);
zonal_vvel_scenarios = cell(1,sz);
zonal_uvel_scenarios = cell(1,sz);
merid_uvel_scenarios = cell(1,sz);
merid_salt_scenarios = cell(1,sz);
merid_theta_scenarios = cell(1,sz);
bsf_scenarios = cell(1,sz);

%load unchanging bathy
fid = fopen(bathy_path);
bathy = fread(fid, 'real*8', 'b');
bathy = reshape(bathy, [nx, ny]);
bathy(bathy == 0) = nan;

%loop over runs
for i = 1:sz
%draft
topo_fname=  ['shelfice_topo_H' num2str(H) '_W' num2str(W) '_extent' num2str(extent(i)) 'km.bin'];
topo_fid = fopen(strcat(topodir, '/',topo_fname));
topo = fread(topo_fid, 'real*8', 'b');
topo = reshape(topo, [nx,ny]);
topo_scenarios{i} = topo;


%melt rates
state2D_fname = strcat(rootdir, run_nos(i), '/run/state2D.nc');
melt = ncread(state2D_fname, 'SHIfwFlx', [1, 1, ntout1], [Inf, Inf, 1+ntout2- ntout1]);
melt = mean(melt, 3); %average over months ntout1 to ntout2
melt = -melt * secs_per_year / density_ice;
melt_scenarios{i} = melt;

%Theta
Theta_fname = strcat(rootdir, run_nos(i), '/run/stateTheta.nc');
Theta = ncread(Theta_fname, 'THETA', [1,1,1,ntout1], [Inf, Inf, Inf, 1+ntout2 - ntout1]);
Theta = mean(Theta, 4);

%Salinity
Salt_fname = strcat(rootdir, run_nos(i), '/run/stateSalt.nc');
Salt = ncread(Salt_fname, 'SALT', [1,1,1,ntout1], [Inf, Inf, Inf, 1+ntout2 - ntout1]);
Salt = mean(Salt, 4);

%Velocities
UVEL_fname = strcat(rootdir, run_nos(i), '/run/stateUvel.nc');
UVEL = ncread(UVEL_fname, 'UVEL', [1,1,1,ntout1], [Inf, Inf, Inf,  1+ntout2 - ntout1]);
UVEL = mean(UVEL, 4);
VVEL_fname = strcat(rootdir, run_nos(i), '/run/stateVvel.nc');
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
bsf_scenarios{i} = stream;

%zonal xsection quantities
TZS = squeeze(Theta(:,zonal_idx, :));
SZS = squeeze(Salt(:, zonal_idx, :));
VZS = squeeze(VVEL(:, zonal_idx, :));
UZS = squeeze(UVEL(:, zonal_idx, :));

%remove bump and topo
for p = 1:nx
for q = 1:nz
        if bathy(p,zonal_idx) > -Z(q)
        TZS(p,q) = nan;
        SZS(p,q) = nan;
        VZS(p,q) = nan;
        UZS(p,q) = nan;
        end
        if topo(p, zonal_idx)< -Z(q)
        TZS(p,q) = nan;
        SZS(p,q) = nan;
        VZS(p,q) = nan;
        UZS(p,q) = nan;
        end
end
end
zonal_theta_scenarios{i} = TZS;
zonal_salt_scenarios{i} = SZS;
zonal_vvel_scenarios{i} = VZS;
zonal_uvel_scenarios{i} = UZS;

%merid section quantities
UMS = squeeze(UVEL(merid_idx, :,:));
VMS = squeeze(VVEL(merid_idx, :,:));
SMS = squeeze(Salt(merid_idx, :,:));
TMS = squeeze(Theta(merid_idx, :, :));
for p = 1:ny
for q = 1:nz
        if bathy(merid_idx,p) > -Z(q)
        TMS(p,q) = nan;
        SMS(p,q) = nan;
	UMS(p,q) = nan;
	VMS(p,q) = nan;
        end
        if topo(merid_idx,p)< -Z(q)
        TMS(p,q) = nan;
        SMS(p,q) = nan;
	UMS(p,q) = nan;
	VMS(p,q) = nan;
        end
end
end
merid_uvel_scenarios{i} = UMS;
merid_salt_scenarios{i} = SMS;
merid_theta_scenarios{i} = TMS;
merid_vvel_scenarios{i} = VMS;

end %end loop over runs
end %end gendata loop

%%%%%%%%%%%%%%%%% Plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

width = 0.45;
ncols = 2;
colgap = 0.07;
startx = (1 -width*ncols - (ncols-1)*colgap)/2;
starty = 0.01;
height = 1/(sz+1);
positions = zeros(4, ncols, sz);
for p = 1:sz
for q = 1:ncols
positions(:,q,p) = [startx + (q-1)*colgap + (q-1)*width, starty + (p-1)*height, width, height];
end
end
colmap = lighter_blue_parula(100,0.2);
colbar_ypos = 0.93;
colbar_xpos = [startx, startx + colgap + width,startx + 2*colgap + 2*width,startx + 3*colgap + 3*width];
colbar_width = width; cbar_fontsize = 12;


% Column 1
%
for p = 1:sz
ax(1,p) = subplot('Position', squeeze(positions(:,1,sz+1-p)));
melt = cell2mat(melt_scenarios(p));
melt = saturate(melt, 120, 0);
topo = squeeze(cell2mat(topo_scenarios(p)));
melt(topo == 0) = nan;
contourf(-Y/1e3,X/1e3,melt, 50, 'linestyle', 'none');
box on
hold on
%add ridge crest and inner cavity
plot(-[50, 50], [0, max(X/1e3)],'--', 'color', [0.5, 0.5, 0.5])
plot(-[30, 30], [0, max(X/1e3)],'w--')
colormap(ax(1,p), colmap)
xticks([]);
yticks([]);
if p == 1
        a = colorbar;
        a.Location = 'northoutside';
        a.Position(1)=colbar_xpos(1);
        a.Position(2)=colbar_ypos;
        a.Label.String = 'melt rate (m/yr)';
        a.FontSize = 10;
        a.Label.FontSize = cbar_fontsize;
        a.Label.Interpreter = 'latex';
end

%add bsf
stream = cell2mat(bsf_scenarios(p));
streamsm = smooth2a(stream, 2,2);
axnew = axes;
axnew.Position = ax(1,p).Position;
[C,h] =contour(-Y/1e3,X/1e3, streamsm, [-0.7, -0.5, -0.3, -0.1], 'k');
%clabel(C,h);
hold on
streamsm(1:4,:) = nan; streamsm(end-3:end,:) = nan; streamsm(:,1:60) = nan; streamsm(:,end-4:end) = nan; %remove borders and near Gl where stream is messy
[C,h] =contour(-Y/1e3,X/1e3, streamsm, [0,0], 'r');
[C,h] =contour(-Y/1e3,X/1e3, streamsm, [0.05,0.05], 'g');


xticks([]);
yticks([]);
set(axnew, 'color', 'none')

end


%
% Column 2
%

for p = 1:sz
ax(2,p) = subplot('Position', squeeze(positions(:,2,sz+1-p)));
TMS = cell2mat(merid_theta_scenarios(p));
TMS = saturate(TMS, 1.3, -1.2);
contourf(max(Y) - Y,-Z, TMS',30, 'linestyle', 'none');
hold on
topo = cell2mat(topo_scenarios(p));
plot(max(Y) - Y, topo(merid_idx, :), 'k', 'linewidth', 1)
%fill the topo
%fx = [max(Y) - Y,flip(max(Y) - Y)];
%fy = [topo(merid_idx,:), zeros(size(topo(merid_idx,:)))];
%fill(fx,fy, 'w')

plot(max(Y) - Y, bathy(merid_idx, :), 'k', 'linewidth', 1)
xlim([4*1e4, Y(end)])
ylim([-1100,-300])
colormap(ax(2,p), colmap)
xticks([]);
yticks([]);
if p == 1
        b = colorbar;
        b.Location = 'northoutside';
        b.Position(1)=colbar_xpos(2);
        b.Position(2)=colbar_ypos;
        b.FontSize = 10;
        b.Label.String = '$\Theta$ (${}^\circ$C)';
        b.Label.FontSize = cbar_fontsize;
        b.Label.Interpreter = 'latex';

end

%add ridge crest and inner cavity
plot(Y(end)-1e3*[50, 50], [-1100,-300],'--', 'color', [0.5, 0.5, 0.5])
plot(Y(end)-1e3*[30, 30], [-1100,-300],'w--')

%add salinity contours
axnew = axes;
axnew.Position = ax(2,p).Position;
SMS = cell2mat(merid_salt_scenarios(p));
[C,h] =contour(max(Y) - Y,-Z , SMS', 34.2:0.2:34.6, 'k');
xticks([]);
yticks([]);
set(axnew, 'color', 'none')
xlim([4*1e4, Y(end)])
ylim([-1100,-300])

end %end loop over runs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if save_flag
saveas(gcf, '../../PIG-calving-melt-response-tex/figures/figureH_melt_response_sections_P600_W200.png')
end

