%Make figure 13 of the PIG calving manuscript: mean melt rate, millgate decomposition for the two regions. Neeed to run make_figure10 first to get the 'snap distance'


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
plotcolor2 = [0,1,1]; %overwrite defaults
plotcolor3 = [1,0,1];
label_size = 11;
ax_fontsize = 11;
figure(3); clf;
fig = gcf; fig.Position(3:4) = [1085, 540];


%
% Data locations
%
rootdir = '/data/oceans_output/shelf/aleey/mitgcm/rPIG/rPIG_'; %output data NOT in github repo (contact for copy)
topodir = '../gendata_realistic/topo_files/';
bathypath = '../gendata_realistic/bathy_files/bathymetry.shice';

%grid details
nx=360; % number of grid cells along longitudinal direction
ny=320; % number of grid cells along latitudinal direction
nz=120; % number of vertical grid cells
dx=400;
dy=400;
dz=10;
X = ncread(strcat(rootdir, "078", '/run/state2D.nc'), 'LONGITUDE');
Y = ncread(strcat(rootdir, "078", '/run/state2D.nc'), 'LATITUDE');
[XX,YY] = meshgrid(X,Y);
YYt = YY';
Z = 0:dz:(nz-1)*dz;


%parameters
secs_per_year = 365.25*24*60*60;
density_ice = 918.0;
lambda = 7.61*1e-4;%constants in liquidus
gamma = 5.73*1e-2;
T0 = 8.32*1e-2;

%time details
ntout1 = 9;
ntout2 = 12; %define time period to average over



%
% Generate data loop
%
if gendata
run_nos = ["078", "082", "083", "084", "085", "086"];
%run_nos = ["079", "088", "089", "090", "091", "092"]; %2012 boundary conditions (use this to compute sensitivity to P in realistic runs)
%run_nos = ["141", "142", "143", "144", "145", "146"];
sz = length(run_nos);

%setup storage
melt_scenarios = cell(1,sz);
topo_scenarios = cell(1,sz);
Tbl_scenarios = cell(1,sz);
Sbl_scenarios = cell(1,sz);
ustar_scenarios = cell(1,sz);

%load bathy
bathyfid = fopen(bathypath);
bathy = fread(bathyfid, 'real*8', 'b');
bathy = reshape(bathy, [nx,ny]);
bathy = double(bathy);

%loop over runs
for i = 1:sz
%draft
topo_fname=  ['shelfice_topo_scn', num2str(i), '.shice'];
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

%u star
ubl = ncread(state2D_fname, 'SHIuLoc', [1, 1, ntout1], [Inf, Inf, 1+ntout2- ntout1]);
vbl = ncread(state2D_fname, 'SHIvLoc', [1, 1, ntout1], [Inf, Inf, 1+ntout2- ntout1]);
ubl = mean(ubl,3);
vbl = mean(vbl,3);
ustar = sqrt(ubl.^2 + vbl.^2);
ustar_scenarios{i} = ustar;

%Theta
Theta_fname = strcat(rootdir, run_nos(i), '/run/stateTheta.nc');
Theta = ncread(Theta_fname, 'THETA', [1,1,1,ntout1], [Inf, Inf, Inf, 1+ntout2 - ntout1]);
Theta = mean(Theta, 4);
Salt_fname = strcat(rootdir, run_nos(i), '/run/stateSalt.nc');
Salt = ncread(Salt_fname, 'SALT', [1,1,1,ntout1], [Inf, Inf, Inf, 1+ntout2 - ntout1]);
Salt = mean(Salt, 4);

Sbl = zeros(nx,ny);
Tbl = zeros(nx,ny);
Nb = 3; %number of grid cells to average over
for p = 1:nx
for q = 1:ny
        if topo(p, q) < 0 %if we're in the cavity
                draft = topo(p,q);
                partial_cell_frac = abs(rem(draft, dz)) / dz;
                draft_rounded = draft + abs(rem(draft, dz));
                [~,idxtop] = min(abs(-Z - draft_rounded));
                vec = [partial_cell_frac,1,1]';
                Sbl(p,q) = sum(vec.*squeeze(Salt(p,q,idxtop:idxtop+Nb-1)))/sum(vec);
                Tbl(p,q) = sum(vec.*squeeze(Theta(p,q,idxtop:idxtop+Nb-1)))/sum(vec);
	end
end %end loop over y
end %end loop over x
Tbl_scenarios{i} = Tbl;
Sbl_scenarios{i} = Sbl;

end %end loop over runs
end

%
% Plots
%
positions = [0.1,0.1,0.4,0.85;
             0.6, 0.55, 0.38, 0.4;
             0.6, 0.1, 0.38, 0.4];


realistic_inner_cavity_definition; %bring inner cavity definition into scope
in1 = inpolygon(XX',YY', a1,b1);
in2 = inpolygon(XX',YY', a2,b2);

topo = cell2mat(topo_scenarios(1));
idx1 = (topo < 0) & in1;
idx2 = (topo < 0) & in2;
for i = 1:sz
melt = cell2mat(melt_scenarios(i));
ave_melt1(i) = mean(melt(idx1));
ave_melt2(i) = mean(melt(idx2));
end
ax(1) = subplot(2,2,1); hold on; box on; grid on
plot(snap_distance(4)*ones(1,2)/1e3, [0, 100], 'k--', 'HandleVisibility', 'off')
plot(snap_distance/1e3, ave_melt1, 'o-', 'color', plotcolor2, 'markerfacecolor', plotcolor2);

ax(3) = subplot(2,2,3); hold on; box on; grid on
plot(snap_distance(4)*ones(1,2)/1e3, [0, 100], 'k--', 'HandleVisibility', 'off')
plot(snap_distance/1e3, ave_melt2, 'o-', 'color', plotcolor3, 'markerfacecolor', plotcolor3);

%tidy these
ax(1).YLim = [39, 44];
ax(1).YTick = 39:44;
ax(3).YLim = [72, 77];
ax(3).YTick = 72:77;
ax(1).XLabel.String = '$l_c$ (km)';
ax(3).XLabel.String = '$l_c$ (km)';
ax(1).YLabel.Interpreter = 'latex';
ax(1).XLabel.Interpreter = 'latex';
ax(3).YLabel.Interpreter = 'latex';
ax(3).XLabel.Interpreter = 'latex';
ax(1).YLabel.String = 'mean melt rate (m/yr)';
ax(3).YLabel.String = 'mean melt rate (m/yr)';
ax(1).XLabel.FontSize = 12;
ax(1).YLabel.FontSize = 12;
ax(3).XLabel.FontSize = 12;
ax(3).YLabel.FontSize = 12;
tN = text(ax(1), 37.5, 39.4, "north box", "Interpreter", "latex", "FontSize", 14, "Color", 'k');
tS = text(ax(3), 1, 76.4, "south box", "Interpreter", "latex", "FontSize", 14, "Color", 'k');
txt2009 = text(ax(1), 1.5,39.3, '2009', 'Interpreter', 'latex', 'FontSize', 10, 'rotation', 90);
txt2020 = text(ax(1), 18, 39.7, '2020', 'Interpreter', 'latex', 'FontSize', 10, 'rotation', 90);

txt2009a = text(ax(3), 1.5, 72.6, '2009', 'Interpreter', 'latex', 'FontSize', 10, 'rotation', 90);
txt2020a = text(ax(3), 18, 74.4, '2020', 'Interpreter', 'latex', 'FontSize', 10, 'rotation', 90);


%
% decompositions
%

%set up storage
relmelt        = zeros(2,sz);
relmelt_noTemp = zeros(2,sz);
relmelt_noVel  = zeros(2,sz);
%get the baseline velocity, theta, salt
topo = cell2mat(topo_scenarios(1));
ustar_baseline = cell2mat(ustar_scenarios(1));
Tbl_baseline = cell2mat(Tbl_scenarios(1));
Sbl_baseline = cell2mat(Sbl_scenarios(1));
Tl_baseline  = T0 + lambda*topo - gamma*Sbl_baseline;
UdT_baseline = ustar_baseline .* (Tbl_baseline - Tl_baseline);

for j = 1:sz

%get the current velocity, theta, salt
ustar = cell2mat(ustar_scenarios(j));
Tbl = cell2mat(Tbl_scenarios(j));
Sbl = cell2mat(Sbl_scenarios(j));
topo = cell2mat(topo_scenarios(j));
Tl  = T0 + lambda*topo - gamma*Sbl; %liquidus temperature in BL

%compute relative melt contributions
UdT = ustar .* (Tbl - Tl);
UdT_noVel =  ustar_baseline .* (Tbl - Tl);
UdT_noTemp = ustar.* (Tbl_baseline - Tl_baseline);


relmelt(1,j) = nanmean(UdT(idx1)) / nanmean(UdT_baseline(idx1));
relmelt_noTemp(1,j) = nanmean(UdT_noTemp(idx1)) / nanmean(UdT_baseline(idx1));
relmelt_noVel(1,j) = nanmean(UdT_noVel(idx1)) / nanmean(UdT_baseline(idx1));

relmelt(2,j) = nanmean(UdT(idx2)) / nanmean(UdT_baseline(idx2));
relmelt_noTemp(2,j) = nanmean(UdT_noTemp(idx2)) / nanmean(UdT_baseline(idx2));
relmelt_noVel(2,j) = nanmean(UdT_noVel(idx2)) / nanmean(UdT_baseline(idx2));

end

plot_defaults %rerun to return to default colours
for p = 1:2
ax(2*p) = subplot(2,2,2*p); hold on;box on; grid on;
plot(snap_distance(4)*ones(1,2)/1e3, [0.95, 1.15], 'k--', 'HandleVisibility', 'off')
plot(snap_distance/1e3, relmelt(p,:), 'o-', 'color', plotcolor1, 'markerfacecolor', plotcolor1);
plot(snap_distance/1e3, relmelt_noTemp(p,:), 'o-', 'color', plotcolor2, 'markerfacecolor', plotcolor2);
plot(snap_distance/1e3, relmelt_noVel(p,:), 'o-', 'color', plotcolor3, 'markerfacecolor', plotcolor3);
end
l = legend(ax(2), {"$\mathcal{M}$", "$U_e$", "$\Delta T_e$"}, 'location', 'northwest','interpreter', 'latex', 'FontSize', 12);

ax(2).YLim = [0.98, 1.15];
ax(4).YLim = [0.98, 1.08];
ax(2).XLabel.String = '$l_c$ (km)';
ax(4).XLabel.String = '$l_c$ (km)';
ax(2).YLabel.Interpreter = 'latex';
ax(2).XLabel.Interpreter = 'latex';
ax(4).YLabel.Interpreter = 'latex';
ax(4).XLabel.Interpreter = 'latex';
ax(2).YLabel.String = 'relative effect';
ax(4).YLabel.String = 'relative effect';
ax(2).XLabel.FontSize = 12;
ax(2).YLabel.FontSize = 12;
ax(4).XLabel.FontSize = 12;
ax(4).YLabel.FontSize = 12;


%plot labels
ta = text(ax(1), -10, 44.4, '(a)', 'Interpreter', 'latex', 'FontSize', 12);
tb = text(ax(2), -10, 1.166, '(b)', 'Interpreter', 'latex', 'FontSize', 12);
tc = text(ax(3), -10, 77.4, '(c)', 'Interpreter', 'latex', 'FontSize', 12);
td = text(ax(4), -10, 1.09, '(d)', 'Interpreter', 'latex', 'FontSize', 12);

% add the inset with the % change in melt rate
axnew = axes;
hold on; box on;
axnew.Position= [0.17, 0.75, 0.12, 0.15];
s1 = sum(sum(idx1)); s2 = sum(sum(idx2)); st = s1 + s2;
pc1 = (ave_melt1 - ave_melt1(1)) / ave_melt1(1) * 100; %percetnage change in box north
pc2 = (ave_melt2 - ave_melt2(1)) / ave_melt2(1) * 100; %south box
pct = pc1 * (s1 / (s1 + s2)) + pc2 * (s2 / (s1 + s2)); %weighted percentage change (mean in inner cavity)
ms = 2;
plot(snap_distance/1e3, pc1, 'co-', 'MarkerSize', ms, 'MarkerFaceColor', 'c');
plot(snap_distance/1e3, pc2, 'mo-', 'MarkerSize', ms, 'MarkerFaceColor', 'm');
plot(snap_distance/1e3, pct, 'ko-', 'MarkerSize', ms, 'MarkerFaceColor', 'k');
xlim([0, 50])
xlabel('$\ell_c$ (km)', 'Interpreter', 'latex');
ylabel('perc. change', 'Interpreter', 'latex');




%
% Save flag
%
set(gcf, 'color', 'w')
if saveflag
%saveas(gcf, 'plots/figure13.eps', 'epsc')
saveas(gcf, 'plots/figure13.png')
end
