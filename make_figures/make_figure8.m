%Make figure 8 in the PIG calving manuscripts:
%(a) Mean melt rate change with calving for P700 (W = 100, 150, 200)
%(b) Mean melt rate change with calving for P800 (W = 100, 150, 200)
%(c) Velocity-thermal driving decompostion for P800, W = 100
%(d) Velocity-thermal driving decompostion for P800, W = 150
%(c) Velocity-thermal driving decompostion for P800, W = 100
%(e) Velocity-thermal driving decompostion for P800, W = 200


%
% NB: Many of the data files referred to in this script are too large to be hosted online. These files are hosted internally as BAS.
% Please email Alex Bradley (aleey@bas.ac.uk) to obtain a copy.

%Alex Bradley (aleey@bas.ac.uk) 27/05/2021. MIT license.

%
% Flags
%
gendata = 1; %specify whether to pass through the generate data loop
saveflag = 0;

%
% Preliminaries
%
addpath("plot_tools");
plot_defaults
label_size = 11;
ax_fontsize = 10;
figure(1); clf;
fig = gcf; fig.Position(3:4) = [1000, 520];


%
% Data locations
%
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
ntout2 = 8; %define time period to average over

run_nos =  ["185", "186", "187", "188", "189", "190", "191", "192", "193", "194";  %W = 200, P = 700
	    "197", "198", "199", "200", "201", "202", "203", "204", "205", "206";  %W = 150, P = 700
	    "207", "208", "209", "210", "211", "212", "213", "214", "215", "216";  %W = 100, P = 700
	    "217", "218", "219", "220", "221", "222", "223", "224", "225", "226";  %W = 200, P = 800
	    "227", "228", "229", "230", "231", "232", "233", "234", "235", "236";  %W = 150, P = 800
	    "237", "238", "239", "240", "241", "242", "243", "244", "245", "246"]; %W = 100, P = 800
sz = size(run_nos);
extent = [84,80,75,70,65,60,55,50,45,40];
W = [200,150, 100, 200, 150, 100]; %ridge gap
P = [700, 700, 700, 800, 800, 800];
H = 400;

%generate data loop
if gendata

%load unchanging bathy
fid = fopen(bathy_path);
bathy = fread(fid, 'real*8', 'b');
bathy = reshape(bathy, [nx, ny]);
bathy(bathy == 0) = nan;

%setup storage
melt_scenarios = cell(sz);
Ubl_scenarios  = cell(sz);
Vbl_scenarios  = cell(sz);
Tbl_scenarios  = cell(sz);
Sbl_scenarios  = cell(sz);
topo_scenarios = cell(sz);

%loop over runs
for i = 1:sz(1)
for j = 1:sz(2)
%draft
topo_fname=  ['shelfice_topo_H' num2str(H) '_W' num2str(W(i)) '_extent' num2str(extent(j)) 'km.bin'];
topo_fid = fopen(strcat(topodir, '/',topo_fname));
topo = fread(topo_fid, 'real*8', 'b');
topo = reshape(topo, [nx, ny]);
topo_scenarios{i,j} = topo;

%melt rates
state2D_fname = strcat(rootdir, run_nos(i,j), '/run/state2D.nc');
melt = ncread(state2D_fname, 'SHIfwFlx', [1, 1, ntout1], [Inf, Inf, 1+ntout2- ntout1]);
melt = mean(melt, 3); %average over months ntout1 to ntout2
melt = -melt * secs_per_year / density_ice;
melt_scenarios{i,j} = melt;


%Theta
Theta_fname = strcat(rootdir, run_nos(i,j), '/run/stateTheta.nc');
Theta = ncread(Theta_fname, 'THETA', [1,1,1,ntout1], [Inf, Inf, Inf, 1+ntout2 - ntout1]);
Theta = mean(Theta, 4);


%Salinity
Salt_fname = strcat(rootdir, run_nos(i,j), '/run/stateSalt.nc');
Salt = ncread(Salt_fname, 'SALT', [1,1,1,ntout1], [Inf, Inf, Inf, 1+ntout2 - ntout1]);
Salt = mean(Salt, 4);

%Velocities
UVEL_fname = strcat(rootdir, run_nos(i,j), '/run/stateUvel.nc');
UVEL = ncread(UVEL_fname, 'UVEL', [1,1,1,ntout1], [Inf, Inf, Inf,  1+ntout2 - ntout1]);
UVEL = mean(UVEL, 4);
VVEL_fname = strcat(rootdir, run_nos(i,j), '/run/stateVvel.nc');
VVEL = ncread(VVEL_fname, 'VVEL', [1,1,1,ntout1], [Inf, Inf, Inf,  1+ntout2 - ntout1]);
VVEL = mean(VVEL, 4);
%boundary layer quantities
Nb = 2; %number of grid pts to take mean over
Sbl = nan(nx,ny); Tbl = nan(nx,ny); Ubl = nan(nx, ny); Vbl = nan(nx,ny);
for p = 1:nx
for q = 1:ny
   if topo(p, q) < 0 %if we're in the cavity
        draft = topo(p,q);
        partial_cell_frac = abs(rem(draft, dz)) / dz;
        draft_rounded = draft + abs(rem(draft, dz));
        [~,idxtop] = min(abs(-Z - draft_rounded));
        vec = [partial_cell_frac,1-partial_cell_frac]'; %issue weights so that it's computed as the average over dz
        Sbl(p,q) = sum(vec.*squeeze(Salt(p,q,idxtop:idxtop+Nb-1)))/sum(vec);
        Tbl(p,q) = sum(vec.*squeeze(Theta(p,q,idxtop:idxtop+Nb-1)))/sum(vec);
        Ubl(p,q) = sum(vec.*squeeze(UVEL(p,q,idxtop:idxtop+Nb-1)))/sum(vec);
        Vbl(p,q) = sum(vec.*squeeze(VVEL(p,q,idxtop:idxtop+Nb-1)))/sum(vec);
    end
end %end loop over y grid
end %end loop over x grid

Sbl_scenarios{i,j} = Sbl;
Tbl_scenarios{i,j} = Tbl;
Ubl_scenarios{i,j} = Ubl;
Vbl_scenarios{i,j} = Vbl;

end %end loop over i
end %end loop over j
end %end generate data loop

%
% Plots
%
positions = [0.06, 0.55, 0.4, 0.43; 
	     0.56, 0.55, 0.4, 0.43;
	     0.06, 0.1, 0.27, 0.33;
	     0.38, 0.1, 0.27, 0.33;
	     0.70, 0.1, 0.27, 0.33];
plotcols = [plotcolor3; plotcolor1; plotcolor2]; %permute for consitency with fig 6
plotcols = parula(4);
%
% (a) decomposition for P 700
%
ax(1) = subplot('Position', positions(1,:));hold on
ylabel('Melt rate (m/yr)', 'FontSize', 12, 'Interpreter', 'latex');
xlabel('$l_c$ (km)', 'FontSize', 12, 'Interpreter', 'latex');
box on

ave_melt = zeros(3,sz(2));
for i = 1:3
for j = 1:sz(2)
melt = cell2mat(melt_scenarios(i,j));
ave_melt(i,j) = mean(melt(idx));
end %end loop over front positions
end %end loop over i = 1:3
m100_700 = ave_melt(3,:);

i = 6;
for j = 1:sz(2)
melt = cell2mat(melt_scenarios(i,j));
m100_800(j) = mean(melt(idx));
end

plot([34,34],  [0,100], 'k--', 'linewidth', 1.5, 'HandleVisibility', 'off'); %plot the location of top of ridge
for i = 1:3
plot(84 - extent, ave_melt(i,:), 'o-', 'color', plotcols(i,:), 'markerfacecolor', plotcols(i,:))
end
ax(1).YLim = [25, 65];
ax(1).XLim = [0,45];
grid on
ax(1).XTick = [0:10:40];
ax(1).YTick = 25:10:65;
P700txt = text(ax(1), 0.3, 27, '$P=700$~m', 'interpreter', 'latex', 'FontSize', 12);

%add sensitivity as an inset
axin = axes;axin.Position = [0.27, 0.62, 0.17, 0.17];
plot(axin, 84 - extent, m100_700 ./ m100_800, 'ko-', 'markerfacecolor', 'k')
axin.XLabel.String = '$l_c$~(km)';
axin.XLabel.Interpreter = 'latex';
axin.YLabel.String= 'sensitivity';
axin.YLabel.Interpreter = 'latex';
axin.YLabel.FontSize  = 10;

%
% (b)
%

ax(2) = subplot('Position', positions(2,:));hold on
ylabel('Melt rate (m/yr)', 'FontSize', 12, 'Interpreter', 'latex');
xlabel('$l_c$ (km)', 'FontSize', 12, 'Interpreter', 'latex');

ave_melt = zeros(3,sz(2));
for i = 1:3
for j = 1:sz(2)
melt = cell2mat(melt_scenarios(i+3,j));
ave_melt(i,j) = mean(melt(idx));
end %end loop over front positions
end %end loop over i = 1:3

plot([34,34],  [0,100], 'k--', 'linewidth', 1.5, 'HandleVisibility', 'off'); %plot the location of top of ridge
for i = 1:3
plot(84 - extent, ave_melt(i,:), 'o-', 'color', plotcols(i,:), 'markerfacecolor', plotcols(i,:))
end
ax(2).YLim = ax(1).YLim;
ax(2).XLim = ax(1).XLim;
grid on
ax(2).XTick = [0:10:40];
ax(2).YTick = ax(1).YTick;
box on
P800txt = text(ax(2), 0.3, 27, '$P=800$~m', 'interpreter', 'latex', 'FontSize', 12);
legend({"$W=200$~m", "$W=150$~m", "$W=100$~m"}, 'location', 'northeast', 'FontSize', 12, 'Interpreter', 'latex');



%
% decompositions
%
ax(3) = subplot('Position', positions(5,:));
xlabel('$l_c$ (km)', 'FontSize', 12, 'Interpreter', 'latex');
hold on
box on
ax(4) = subplot('Position', positions(4,:));
xlabel('$l_c$ (km)', 'FontSize', 12, 'Interpreter', 'latex');
box on
hold on
ax(5) = subplot('Position', positions(3,:)); %switch to W = 100 on left
xlabel('$l_c$ (km)', 'FontSize', 12, 'Interpreter', 'latex');
ylabel('Relative change', 'FontSize', 12, 'Interpreter', 'latex');
box on
hold on

for i = 1:3
%set up storage
relmelt        = zeros(1,sz(2));
relmelt_noTemp = zeros(1,sz(2));
relmelt_noVel  = zeros(1,sz(2));


%get the baseline velocity, theta, salt
Ubl_baseline = cell2mat(Ubl_scenarios(i+3,1));
Vbl_baseline = cell2mat(Vbl_scenarios(i+3,1));
Tbl_baseline = cell2mat(Tbl_scenarios(i+3,1));
Sbl_baseline = cell2mat(Sbl_scenarios(i+3,1));
Tl_baseline  = T0 + lambda*cell2mat(topo_scenarios(i+3,1)) - gamma*Sbl_baseline;
UdT_baseline = sqrt(Ubl_baseline.^2 + Vbl_baseline.^2) .* (Tbl_baseline - Tl_baseline);

for j = 1:sz(2)

%get the current velocity, theta, salt
Ubl = cell2mat(Ubl_scenarios(i+3,j));
Vbl = cell2mat(Vbl_scenarios(i+3,i));
Tbl = cell2mat(Tbl_scenarios(i+3,j));
Sbl = cell2mat(Sbl_scenarios(i+3,j));
Tl  = T0 + lambda*topo - gamma*Sbl; %liquidus temperature in BL

%compute relative melt contributions
UdT = sqrt(Ubl.^2 + Vbl.^2) .* (Tbl - Tl);
UdT_noVel =  sqrt(Ubl_baseline.^2 + Vbl_baseline.^2) .* (Tbl - Tl);
UdT_noTemp = sqrt(Ubl.^2 + Vbl.^2) .* (Tbl_baseline - Tl_baseline);

relmelt(j) = nanmean(UdT(idx)) / nanmean(UdT_baseline(idx));
relmelt_noTemp(j) = nanmean(UdT_noTemp(idx)) / nanmean(UdT_baseline(idx));
relmelt_noVel(j) = nanmean(UdT_noVel(idx)) / nanmean(UdT_baseline(idx));
end %end loop over runs

plot(ax(i+2), [34,34],  [0.6, 1.6], 'k--', 'linewidth', 1.5, 'handlevisibility', 'off'); %plot the location of top of ridge
plot(ax(i+2), 84 - extent, relmelt, 'o-', 'color', plotcolor1, 'markerfacecolor', plotcolor1);
plot(ax(i+2),84 - extent, relmelt_noTemp, 'o-', 'color', plotcolor2, 'markerfacecolor', plotcolor2);
plot(ax(i+2), 84 - extent, relmelt_noVel, 'o-', 'color', plotcolor3, 'markerfacecolor', plotcolor3);
ax(i+2).XLim = [0, 45];
grid(ax(i+2), 'on')
end

%tidy
ax(3).YLim = [0.8, 1.2];
ax(4).YLim = [0.8, 1.2];
ax(5).YLim = [0.8, 1.6];
legend(ax(3),{"$\mathcal{M}$", "$U_e$", "$\Delta T_e$"}, 'location', 'southwest','interpreter', 'latex', 'FontSize', 12)
W100txt = text(ax(5), 0.5, 1.55, '$W=100$~m', 'Interpreter', 'latex', 'FontSize', 11.5);
W150txt = text(ax(4), 0.5, 1.175, '$W=150$~m', 'Interpreter', 'latex', 'FontSize', 11.5);
W200txt = text(ax(3), 0.5, 1.175, '$W=200$~m', 'Interpreter', 'latex', 'FontSize', 11.5);

%plot labels
ta = text(ax(1), -6.5, 65, '(a)','Interpreter', 'latex', 'FontSize', 12);
tb = text(ax(2), -6.5, 65, '(b)','Interpreter', 'latex', 'FontSize', 12);
tc = text(ax(5), -10, 1.69, '(c)', 'Interpreter', 'latex', 'FontSize', 12);
td = text(ax(4), -8, 1.24, '(d)', 'Interpreter', 'latex', 'FontSize', 12);
te = text(ax(3), -8, 1.24, '(e)', 'Interpreter', 'latex', 'FontSize', 12);


%v
% save?
%
set(gcf, 'color', 'w')
if saveflag
%saveas(gcf, 'plots/figure8.eps', 'epsc');
saveas(gcf, "plots/figure8.png");
end

%%

