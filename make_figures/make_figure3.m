%Make figure 3 in the Idealized PIG calving manuscript: results from the default run.
%(a) Plot of melt rate with boundary current overlain
%(b) Plot of the water column thickness with barotropic stream function overlain
%(c) Plot of the bottom temperature with bottom current overlain
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
fig = gcf; fig.Position(3:4) = [1085, 540];
NS_idx = 60; %index of line running North south
EW_idx = 50; %index of line running east west

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
width = 0.11;
widthsect = 0.25;
gap = 0.025;
startx = (1 - 5*width - 4*gap)/2 - 0.05;
startx = 0.05 ; %
starty = 0.085;
height = 0.8;
heightsmall = 0.35; 
positions = [startx, starty, width, height;
	     startx + gap + width, starty, width, height;
	     startx + 2*gap + 2*width, starty, width, height;
	     startx + 3*gap + 3*width, starty, width, height;
	     startx + 4*gap + 4*width + 0.05, starty, widthsect, heightsmall;
	     startx + 4*gap + 4*width + 0.05, starty + height/2+0.045, widthsect, heightsmall];
	    
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
plot(X/1e3, Y(EW_idx)*ones(length(X))/1e3, 'w', 'linewidth', 1.5)
plot(X(NS_idx)*ones(1,230)/1e3, Y(1:230)/1e3, 'm--', 'linewidth', 1.5);

%colorbar and arrow
A = colorbar;
A.Location = 'northoutside';
A.Label.String = 'melt rate (m/yr)';
A.Label.Interpreter = 'latex';
A.Position(end) = 0.02;%A.Position(end) - 0.02;
A.Position(2) = 0.89;
A.Label.FontSize = 11;
%plot([20, 30], [120,120], 'k', 'linewidth', 1);
%text(31, 120, '0.6 m/s')
xlabel('$x$~(km)', 'Interpreter', 'latex', 'FontSize' ,12);
ylabel('$y$~(km)', 'Interpreter', 'latex', 'FontSize' ,12)
text(-10,141, '(a)', 'Interpreter', 'latex', 'FontSize',  12)

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
b = colorbar;
b.Location = 'northoutside';
b.Position(end) = 0.02; %b.Position(end) - 0.02;
b.Position(2) = 0.89;%b.Position(2) + 0.03;
b.Label.String = 'Bottom temp~(${}^\circ$C)';
b.Label.Interpreter = 'latex';
b.Label.FontSize = 11;
yticks([])
quiver(XX(idxY, idxX),YY(idxY, idxX),velscale *Ubot(idxX, idxY)', velscale*Vbot(idxX, idxY)', 'autoscale', 'off', 'color', 'k')
xlabel('$x$~(km)', 'Interpreter', 'latex', 'FontSize' ,12);

text(-10,141, '(b)', 'Interpreter', 'latex', 'FontSize',  12)

%
% Plot 3: 1/h and BSF contours
%
ax(3) = subplot('Position', positions(3,:)); hold on; box on;
ax(3).FontSize = 10;
column_thickness = topo - bathy;
plot(X/1e3, 50*ones(length(X), 1), 'w--', 'linewidth', 1.5)
contourf(X/1e3,Y/1e3,1e3* (1./column_thickness)', 20, 'linestyle', 'none');
colormap(ax(3), cmap);
a = colorbar;
a.Location = 'northoutside';
a.Position(end) = 0.02;%a.Position(end) - 0.02;
a.Position(2) = 0.89;%a.Position(2) + 0.03;
a.Label.Interpreter = 'latex';
a.Label.String = '$1/h~(10^{-3}~\mathrm{m}^{-1})$';
a.Label.FontSize = 11;
yticks([])
xlabel('$x$~(km)', 'Interpreter', 'latex', 'FontSize' ,12);
text(-10,141, '(c)', 'Interpreter', 'latex', 'FontSize',  12)

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
set(axnew, 'color', 'none')


%
% Plot 4: f + zeta / H
%
ax(4) = subplot('Position', positions(4,:)); hold on; box on
ax(4).FontSize = 10;
contourf(X/1e3, Y/1e3,1e7* BPV', 20, 'linestyle', 'none'); 
colormap(ax(4), cmap);
cc = colorbar;
cc.Location = 'northoutside';
cc.Position(end) = 0.02; %b.Position(end) - 0.02;
cc.Position(2) = 0.89;%b.Position(2) + 0.03;
cc.Label.String = '$10^{7}\times[(f + \zeta )/ h]$';
cc.Label.Interpreter = 'latex';
cc.Label.FontSize = 11;
yticks([])
ax(4).XTick = [0,20,40];
xlabel('$x$~(km)', 'Interpreter', 'latex', 'FontSize' ,12);
text(-10,141, '(d)', 'Interpreter', 'latex', 'FontSize',  12)

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
ax(5) = subplot('Position', positions(5,:)); hold on; box on;
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
c = colorbar;
c.Location = 'north';
c.Position(1) = positions(5,1) + 0.02;
c.Position(3) = widthsect - 0.04;
c.Position(4) = 0.02;
c.Position(2) = 0.4;
c.Label.String = '$\Theta$~(${}^\circ$C)';
c.Label.Interpreter = 'latex';
c.Label.FontSize = 11;
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
%shift labels because we plotted in reverse
ax(5).XTick = max(Y)/1e3 - (80:-20:0);
ax(5).XTickLabels = {"80", "60","40", "20", "0"};

text(18,-300, '(f)', 'Interpreter', 'latex', 'FontSize',  12)

%
% Plot 6
% 
ax(6) = subplot('Position', positions(6,:)); hold on; box on;
ax(6).FontSize = 10;
SZS = squeeze(Salt(:, NS_idx,:));
TZS = squeeze(Theta(:, NS_idx, :));
for p = 1:nx
for q = 1:nz
        if bathy(p,NS_idx) > -Z(q)
        TZS(p,q) = nan;
        SZS(p,q) = nan;
        end
        if topo(p,NS_idx)< -Z(q)
        TZS(p,q) = nan;
        SZS(p,q) = nan;
        end
end
end
TZS  = saturate(TZS, max(max(TMS)), min(min(TMS)));
contourf(X/1e3,-Z, TZS',30, 'linestyle', 'none');
colormap(ax(6), cmap);
hold on
ylabel('depth (m)', 'Interpreter', 'latex', 'FontSize', 12);
xlabel('$x$~(km)', 'Interpreter', 'latex', 'FontSize' ,12);
d = colorbar;
d.Location = 'northoutside';
d.Position(1) = positions(5,1) + 0.02;
d.Position(2) = 0.89;
d.Label.String = '$\Theta$~(${}^\circ$C)';
d.Label.Interpreter = 'latex';
d.Label.FontSize = 11;
d.Position(3) = widthsect - 0.04;
d.Position(4) = 0.02;
ylim([bathy(3,NS_idx),topo(3,NS_idx)])
axnew = axes;
axnew.Position = ax(6).Position;
[C,h] =contour(X/1e3,-Z , SZS', 34.2:0.2:34.6, 'k');
clabel(C,h);

xticks([]);
yticks([]);
set(axnew, 'color', 'none')
ylim([bathy(3,NS_idx),topo(3,NS_idx)])

txte = text(-12,-600, '(e)', 'Interpreter', 'latex', 'FontSize',  12);

%
% Save
%
if save_flag
%saveas(gcf, "plots/figure3", 'epsc')
saveas(gcf, "plots/figure3.png")
end


