%Make figure 3 in the Idealized PIG calving manuscript: results from the default run.
%(a) Plot of melt rate with boundary current overlain
%(b) Plot the bottom temperature with the bottom current overlain
%(c) Plot of the water column thickness with barotropic stream function overlain
%(d) Plot of the barotropic potential vorticity with barotropic stream function overlain
%(e) Zonal cross section of temperature near the GL with 34.2, 34.4, 34.6 salinity contours 
%(f) Meridional cross section of temperature along centre of the domain with 34.2, 34.4, 34.6 salinity contours.
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
fig = gcf; fig.Position(3:4) = [900,900];
NS_idx = 125; %index of line for zonal cross section
EW_idx = 60; %index of line for meridional cross section

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
run_nos = "077"; 
sz = length(run_nos);
extent = 84;
H = 400; %ridge height (always 400);
W = 100; %ridge gap

%generate data loop
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

%melt rates
state2D_fname = strcat(rootdir, run_nos, '/run/state2D.nc');
melt = ncread(state2D_fname, 'SHIfwFlx', [1, 1, ntout1], [Inf, Inf, 1+ntout2- ntout1]);
melt = mean(melt, 3); %average over months ntout1 to ntout2
melt = -melt * secs_per_year / density_ice;
melt(topo == 0) = nan;

%Theta
Theta_fname = strcat(rootdir, run_nos, '/run/stateTheta.nc');
Theta = ncread(Theta_fname, 'THETA', [1,1,1,ntout1], [Inf, Inf, Inf, 1+ntout2 - ntout1]);
Theta = mean(Theta, 4);

%Salinity
Salt_fname = strcat(rootdir, run_nos, '/run/stateSalt.nc');
Salt = ncread(Salt_fname, 'SALT', [1,1,1,ntout1], [Inf, Inf, Inf, 1+ntout2 - ntout1]);
Salt = mean(Salt, 4);

%Velocities
UVEL_fname = strcat(rootdir, run_nos, '/run/stateUvel.nc');
UVEL = ncread(UVEL_fname, 'UVEL', [1,1,1,ntout1], [Inf, Inf, Inf,  1+ntout2 - ntout1]);
UVEL = mean(UVEL, 4);
VVEL_fname = strcat(rootdir, run_nos, '/run/stateVvel.nc');
VVEL = ncread(VVEL_fname, 'VVEL', [1,1,1,ntout1], [Inf, Inf, Inf,  1+ntout2 - ntout1]);
VVEL = mean(VVEL, 4);


%boundary layer quantities
Nb = 3; %number of grid pts to take mean over
Sbl = nan(nx,ny); Tbl = nan(nx,ny); Ubl = nan(nx, ny); Vbl = nan(nx,ny);
Sbot = nan(nx,ny); Tbot = nan(nx,ny); Ubot = nan(nx,ny); Vbot = nan(nx,ny);
for p = 1:nx
for q = 1:ny
	%work out the bottom and velocity
	idx = find((bathy(p,q) - (-Z) < 0), 1, 'last'); %gives you the index of first Z grid point above the bathymetry
         Tbot(p,q) = double(mean(Theta(p,q,idx-Nb+1:idx)));
         Ubot(p,q) = double(mean(UVEL(p,q,idx-Nb+1:idx)));
         Vbot(p,q) = double(mean(VVEL(p,q,idx-Nb+1:idx)));


        if topo(p, q) < 0 %if we're in the cavity
                idxtop = find((topo(p,q) - (-Z)) > 0, 1, 'first'); %gives you the index of first Z grid point above the bathymetry
                idxtop = find(Theta(p,q,:) ~= 0);
                idxtop = idxtop(1);
                Sbl(p,q) = double(mean(Salt(p,q,idxtop:idxtop+Nb-1)));
                Tbl(p,q) = double(mean(Theta(p,q,idxtop:idxtop+Nb-1)));
                Ubl(p,q) = double(mean(UVEL(p,q,idxtop:idxtop+Nb-1)));
                Vbl(p,q) = double(mean(VVEL(p,q,idxtop:idxtop+Nb-1)));

                if 1 %account for partial cell in the mean calculation
                draft = topo(p,q);
                partial_cell_frac = abs(rem(draft, dz)) / dz;
                draft_rounded = draft + abs(rem(draft, dz));
                [~,idx_top] = min(abs(-Z - draft_rounded));
                vec = [partial_cell_frac,1,1]';
                Sbl(p,q) = sum(vec.*squeeze(Salt(p,q,idxtop:idxtop+Nb-1)))/sum(vec);
                Tbl(p,q) = sum(vec.*squeeze(Theta(p,q,idxtop:idxtop+Nb-1)))/sum(vec);

                Ubl(p,q) = sum(vec.*squeeze(UVEL(p,q,idxtop:idxtop+Nb-1)))/sum(vec);
                Vbl(p,q) = sum(vec.*squeeze(VVEL(p,q,idxtop:idxtop+Nb-1)))/sum(vec);
                end

        end
end %end loop over y grid
end %end loop over x grid

%compute BSF
vvel = squeeze(sum(VVEL, 3)) * dz; %units m^2 /s
stream=zeros(size(vvel));
stream(nx,:)=vvel(nx,:)*dx;
for p=nx-1:-1:1
 stream(p,:)=stream(p+1,:) + vvel(p,:)*dx;
end
stream = stream/1e6; %convert to sv

end %end gendata loop

% Compute the transport across the ridge
vvel_ridge = VVEL(:,125,:);  %125 is index of ridge crest
vvel_ridge_neg = vvel_ridge; vvel_ridge_neg(vvel_ridge_neg > 0) = 0; %set positive entries (out of cavity) to zero
neg_flux = sum(sum(vvel_ridge_neg))*dx*dy/1e9; %flux into the cavity in Sv
fprintf("total flux into the cavity is %.3f Sv \n", neg_flux)

%%%% potential vorticity
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plots %%%%%%%%%%%%%%%%%%%%%%%
% generate plot positions, starting with top four
width     = 0.20;
gap       = 0.04;
startx    = 0.07;
starty    = 0.4;
height    = 0.5;
positions = [startx, starty, width, height;
	     startx + gap + width, starty, width, height;
	     startx + 2*gap + 2*width, starty, width, height;
	     startx + 3*gap + 3*width, starty, width, height];
	    
width2 = 0.4;
height2 = 0.25; 
start2 = 0.08;
gap2  = 0.1;
positions = [positions;
	    start2, start2, width2*0.9, height2;
	    start2 + width2 + gap2, start2, width2, height2]; 

%
% Plot 1: Melt rate and BL velocities
%
ax(1) = subplot('Position', positions(1,:)); hold on; box on
ax(1).FontSize = 10;
contourf(X/1e3,Y/1e3,melt', 50, 'linestyle', 'none');
cmap = lighter_blue_parula(100,0.2); 
colormap(ax(1), cmap);

%add the velocity arrows
idxX = 1:8:120;
idxY = 1:8:320;
[XX,YY] = meshgrid(X,Y);
XX = XX/1e3; YY = YY/1e3;
velscale =15;
quiver(XX(idxY, idxX),YY(idxY, idxX),velscale *Ubl(idxX, idxY)', velscale*Vbl(idxX, idxY)', 'autoscale', 'off', 'color', 'k')
plot(X/1e3, 50*ones(length(X), 1), 'w--', 'linewidth', 1.5)
%plot(X/1e3, Y(NS_idx)*ones(length(X))/1e3, 'w', 'linewidth', 1.5)
plot(X(EW_idx)*ones(1,230)/1e3, Y(1:230)/1e3, 'm--', 'linewidth', 1.5);

%colorbar and arrow
cbar(1) = colorbar;
cbar(1).Location = 'northoutside';
cbar(1).Label.String = 'melt rate (m/yr)';
cbar(1).Label.Interpreter = 'latex';
cbar(1).Position(end) = 0.01;
cbar(1).Position(2) = 0.92;
cbar(1).Label.FontSize = 11;
%plot([20, 30], [120,120], 'k', 'linewidth', 1);
%text(31, 120, '0.6 m/s')
xlabel('$x$~(km)', 'Interpreter', 'latex', 'FontSize' ,12);
ylabel('$y$~(km)', 'Interpreter', 'latex', 'FontSize' ,12)
txt(1) = text(-8,130, '(a)', 'Interpreter', 'latex', 'FontSize',  12);

%
% Plot 2: Bottom temp
%
ax(2) = subplot('Position', positions(2,:)); hold on; box on
ax(2).FontSize = 10;
Tbotsm = smooth2a(Tbot, 2,2);
contourf(X/1e3,Y/1e3,Tbotsm', 50, 'linestyle', 'none');
cmap = lighter_blue_parula(100,0.2); 
colormap(ax(2), cmap);

plot(X/1e3, 50*ones(length(X), 1), 'w--', 'linewidth', 1.5)
cbar(2) = colorbar;
cbar(2).Location = 'northoutside';
cbar(2).Position(end) = 0.01;
cbar(2).Position(2) = 0.92;
cbar(2).Label.String = 'Bottom temp~(${}^\circ$C)';
cbar(2).Label.Interpreter = 'latex';
cbar(2).Label.FontSize = 11;
yticks([])
quiver(XX(idxY, idxX),YY(idxY, idxX),velscale *Ubot(idxX, idxY)', velscale*Vbot(idxX, idxY)', 'autoscale', 'off', 'color', 'k')
xlabel('$x$~(km)', 'Interpreter', 'latex', 'FontSize' ,12);

txt(2) = text(-8,130, '(b)', 'Interpreter', 'latex', 'FontSize',  12);

%
% Plot 3: 1/h and BSF contours
%
ax(3) = subplot('Position', positions(3,:)); hold on; box on;
ax(3).FontSize = 10;
column_thickness = topo - bathy;
plot(X/1e3, 50*ones(length(X), 1), 'w--', 'linewidth', 1.5)
contourf(X/1e3,Y/1e3,1e3* (1./column_thickness)', 20, 'linestyle', 'none');
colormap(ax(3), cmap);
cbar(3) = colorbar;
cbar(3).Location = 'northoutside';
cbar(3).Position(end) = 0.01;
cbar(3).Position(2) = 0.92;
cbar(3).Label.Interpreter = 'latex';
cbar(3).Label.String = '$1/h~(10^{-3}~\mathrm{m}^{-1})$';
cbar(3).Label.FontSize = 11;
yticks([])
xlabel('$x$~(km)', 'Interpreter', 'latex', 'FontSize' ,12);
txt(3) =text(-8,130, '(c)', 'Interpreter', 'latex', 'FontSize',  12);

%add bsf
streamsm = smooth2a(stream, 2,2);
axnew = axes;
axnew.Position = ax(3).Position;
[C,h] =contour(X/1e3,Y/1e3, streamsm', [-0.7, -0.5, -0.3, -0.1], 'k');
clabel(C,h);
hold on
streamsm(1:4,:) = nan; streamsm(end-3:end,:) = nan; streamsm(:,1:20) = nan; streamsm(:,end-4:end) = nan; %remove borders and near Gl where stream is messy
[C,h] =contour(X/1e3,Y/1e3, streamsm', [0,0], 'r');
clabel(C,h);
[C,h] =contour(X/1e3,Y/1e3, streamsm', [0.05,0.05], 'g');
clabel(C,h);

xticks([]);
yticks([]);
axnew.YTick = ax(1).YTick;
set(axnew, 'color', 'none')


%
% Plot 4: f + zeta / H
%
ax(4) = subplot('Position', positions(4,:)); hold on; box on
ax(4).FontSize = 10;
contourf(X/1e3, Y/1e3,1e7* BPV', 20, 'linestyle', 'none'); 
colormap(ax(4), cmap);
cbar(4) = colorbar;
cbar(4).Location = 'northoutside';
cbar(4).Position(end) = 0.01;
cbar(4).Position(2) = 0.92;
cbar(4).Label.String = '$10^{7}\times[(f + \zeta )/ h]$';
cbar(4).Label.Interpreter = 'latex';
cbar(4).Label.FontSize = 11;
yticks([])
ax(4).XTick = [0,20,40];
xlabel('$x$~(km)', 'Interpreter', 'latex', 'FontSize' ,12);
txt(4) = text(-8,130, '(d)', 'Interpreter', 'latex', 'FontSize',  12);

%add bsf contours
axnew = axes;
axnew.Position = ax(4).Position;
[C,h] =contour(X/1e3,Y/1e3, streamsm', [-0.7, -0.5, -0.3, -0.1], 'k');
clabel(C,h);
hold on
[C,h] =contour(X/1e3,Y/1e3, streamsm', [0,0], 'r');
clabel(C,h);
[C,h] =contour(X/1e3,Y/1e3, streamsm', [0.05,0.05], 'g');
clabel(C,h);
xticks([]);
yticks([]);
set(axnew, 'color', 'none')



%
% Plot 5: meridional cross section
%
ax(5) = subplot('Position', positions(6,:)); hold on; box on;
ax(5).FontSize = 10;
SMS = squeeze(Salt(EW_idx, :,:));
TMS = squeeze(Theta(EW_idx, :, :));
for p = 1:ny
for q = 1:nz
        if bathy(EW_idx,p) > -Z(q)
        TMS(p,q) = nan;
        SMS(p,q) = nan;
        end
        if topo(EW_idx,p)< -Z(q)
        TMS(p,q) = nan;
        SMS(p,q) = nan;
        UMS(p,q) = nan;
        VMS(p,q) = nan;
        end
end
end
contourf(max(Y)/1e3 - Y/1e3,-Z, TMS',30, 'linestyle', 'none');
colormap(ax(5), cmap);
hold on
plot(max(Y)/1e3 - Y/1e3, topo(EW_idx, :), 'k', 'linewidth', 1)
plot(max(Y)/1e3 - Y/1e3, bathy(EW_idx, :), 'k', 'linewidth', 1)
xlim([4*1e4, Y(end)]/1e3)
ylim([-1100,-300])
ylabel('depth (m)', 'Interpreter', 'latex', 'FontSize', 12);
xlabel('$y$~(km)', 'Interpreter', 'latex', 'FontSize' ,12);
cbar(5) = colorbar;
cbar(5).Location = 'north';
cbar(5).Position = [0.65, 0.26, 0.3,0.01];
cbar(5).Label.String = '$\Theta$~(${}^\circ$C)';
cbar(5).Label.Interpreter = 'latex';
cbar(5).Label.FontSize = 11;
%plot((max(Y)/1e3 - 20)*[1,1], [bathy(3,NS_idx),topo(3,NS_idx_idx)], 'm--', 'linewidth', 1.5)


axnew = axes;
axnew.Position = ax(5).Position;
[C,h] =contour(max(Y)/1e3 - Y/1e3,-Z , SMS', 34.2:0.2:34.6, 'k');
%clabel(C,h);

xticks([]);
yticks([]);
set(axnew, 'color', 'none')
xlim([4*1e4, Y(end)]/1e3)
ylim([-1100,-300])
ax(5).YTick = -1000:200:-400;
%shift labels because we plotted in reverse
ax(5).XTick = max(Y)/1e3 - (80:-20:0);
ax(5).XTickLabels = {"80", "60","40", "20", "0"};

txt(5) = text(30,-280, '(f)', 'Interpreter', 'latex', 'FontSize',  12);

%
% Plot (e): zonal cross section of temperature and salinity
% 
ax(6) = subplot('Position', positions(5,:)); hold on; box on;
ax(6).FontSize = 10;
SZS = squeeze(Salt(:, NS_idx,:));
TZS = squeeze(Theta(:, NS_idx, :));
UZS = squeeze(UVEL(:,NS_idx,:));
VZS = squeeze(VVEL(:,NS_idx,:));
for p = 1:nx
for q = 1:nz
        if bathy(p,NS_idx) > -Z(q)
        TZS(p,q) = nan;
        SZS(p,q) = nan;
        UZS(p,q) = nan;
        VZS(p,q) = nan;
        end
        if topo(p,NS_idx)< -Z(q)
        TZS(p,q) = nan;
        SZS(p,q) = nan;
        UZS(p,q) = nan;
        VZS(p,q) = nan;
        end
end
end
TZS  = saturate(TZS, max(max(TMS)), min(min(TMS)));
UZS = saturate(UZS,0.22, -0.05);
contourf(X/1e3,-Z, UZS',60, 'linestyle', 'none');
colormap(ax(6), cmap);
hold on
ylabel('depth (m)', 'Interpreter', 'latex', 'FontSize', 12);
xlabel('$x$~(km)', 'Interpreter', 'latex', 'FontSize' ,12);

cbar(6) = colorbar;
%colormap(ax(6), redblue)
%cbar(6).Location = 'northoutside';
%cbar(6).Position(1) = positions(5,1) + 0.02;
%cbar(6).Position(2) = 0.89;
cbar(6).Label.String = '$u~(\mathrm{m}~\mathrm{s}^{-1}$ )';
cbar(6).Label.Interpreter = 'latex';
cbar(6).Label.FontSize = 11;
%cbar(6).Position(3) = widthsect - 0.04;
%cbar(6).Position(4) = 0.02;

ax(6).YLim = [-700,-600];
ax(6).YTick = -700:50:-600;
%ylim([bathy(3,NS_idx),topo(3,NS_idx)])
%axnew = axes;
%axnew.Position = ax(6).Position;
%[C,h] =contour(X/1e3,-Z , SZS', 34.2:0.2:34.6, 'k');
%clabel(C,h);

%set(axnew, 'color', 'none')
%ylim([bathy(3,NS_idx),topo(3,NS_idx)])

txt(6) = text(-12,-596, '(e)', 'Interpreter', 'latex', 'FontSize',  12);

%
% Save
%
if save_flag
%saveas(gcf, "plots/figure3", 'epsc')
saveas(gcf, "plots/figure3.png")
end


