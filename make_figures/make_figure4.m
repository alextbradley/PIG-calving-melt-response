%Make figure 4 in the idealized PIG calving manuscript: plots of (a) mean melt rate with calving and (b) decomposition of changes in melting
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
label_size = 13;
ax_fontsize = 13;
figure(1); clf; 
fig = gcf; fig.Position(3:4) = [920, 390];

%
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
run_nos = ["077", "078", "079", "080", "081", "082", "083", "084", "085", "086"];
sz = length(run_nos);
extent = [84,80,75,70,65,60,55,50,45,40];
H = 400; %ridge height (always 400);
W = 100; %ridge gap

%generate data loop
if gendata

%load unchanging bathy
fid = fopen(bathy_path);
bathy = fread(fid, 'real*8', 'b');
bathy = reshape(bathy, [nx, ny]);
bathy(bathy == 0) = nan;

%setup storage
melt_scenarios = cell(1,sz);
Ubl_scenarios  = cell(1,sz);
Vbl_scenarios  = cell(1,sz);
Tbl_scenarios  = cell(1,sz);
Sbl_scenarios  = cell(1,sz);
topo_scenarios = cell(1,sz); 
uvel_baro_scenarios = cell(1,sz);
vvel_baro_scenarios = cell(1,sz);


%loop over runs
for i = 1:sz
%draft
topo_fname=  ['shelfice_topo_H' num2str(H) '_W' num2str(W) '_extent' num2str(extent(i)) 'km.bin'];
topo_fid = fopen(strcat(topodir, '/',topo_fname));
topo = fread(topo_fid, 'real*8', 'b');
topo = reshape(topo, [nx, ny]);
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

%boundary layer quantities
Nb = 2; %number of grid pts to take mean over (hard coded)
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

Sbl_scenarios{i} = Sbl;
Tbl_scenarios{i} = Tbl;
Ubl_scenarios{i} = Ubl;
Vbl_scenarios{i} = Vbl;

%compute barotropic velocities
vvel_barotropic = zeros(nx,ny);
uvel_barotropic = zeros(nx,ny);
for p = 1:nx
for q = 1:ny
vvel = squeeze(VVEL(p,q,:));
vvel_barotropic(p,q) = mean(vvel(vvel ~= 0));
uvel = squeeze(UVEL(p,q,:));
uvel_barotropic(p,q) = mean(uvel(uvel ~= 0));

end
end
vvel_baro_scenarios{i} = vvel_barotropic;
uvel_baro_scenarios{i} = uvel_barotropic;
end




end %end generate data loop

%
% Plots
%


%
% Plot 1: Mean inner cavity melt rate with calving
%
Positions = [0.07, 0.12, 0.41, 0.76;
	    0.57, 0.12, 0.41, 0.76];

axa= subplot('Position',Positions(1,:)); grid on; hold on; ax = gca; box on
ave_melt = zeros(1,sz);
for i = 1:sz
melt = cell2mat(melt_scenarios(i));
ave_melt(i) = mean(melt(idx));
end
plot([34,34],  [40,80], 'k--', 'linewidth', 1.5); %plot the location of top of ridge
plot(84 - extent, ave_melt, 'o-', 'color', plotcolor1, 'markerfacecolor', plotcolor1);
xlabel('$l_c$ (km)', 'Interpreter', 'Latex','FontSize', ax_fontsize);
ylabel('Inner cavity melt rate (m/yr)', 'Interpreter', 'latex','FontSize', ax_fontsize);
xlim([0, 84 - 40]);
txa = text(-7, 84, '(a)', 'FontSize', ax_fontsize, 'Interpreter', 'latex');

%add axis with yf on
axa.XTick = [0,10,20,30,40];
axa2 = axes; axa2.Position = axa.Position; axa2.Color = 'none'; axa2.XAxisLocation = 'top';
axa2.XTick = [0,10,20,30,40]/44;
axa2.YTick = [];
axa2.XTickLabel = {"84", "74", "64", "54", "44"};
axa2.XLabel.String = '$y_f$ (km)';
axa2.XLabel.Interpreter= 'latex';
axa2.XLabel.FontSize = ax_fontsize;

%
% Plot 2: Decomposition
%
axb = subplot('Position', Positions(2,:)); grid on; hold on; ax = gca; box on

%set up storage
relmelt        = zeros(1,sz);
relmelt_noTemp = zeros(1,sz);
relmelt_noVel  = zeros(1,sz);

%get the baseline velocity, theta, salt
Ubl_baseline = cell2mat(Ubl_scenarios(1));
Vbl_baseline = cell2mat(Vbl_scenarios(1));
Tbl_baseline = cell2mat(Tbl_scenarios(1));
Sbl_baseline = cell2mat(Sbl_scenarios(1));
Ubaro_baseline = cell2mat(uvel_baro_scenarios(1));
Vbaro_baseline = cell2mat(vvel_baro_scenarios(1));


Tl_baseline  = T0 + lambda*topo - gamma*Sbl_baseline;
UdT_baseline = sqrt(Ubl_baseline.^2 + Vbl_baseline.^2) .* (Tbl_baseline - Tl_baseline);
UdT_baro_baseline = sqrt(Ubaro_baseline.^2 + Vbaro_baseline.^2) .* (Tbl_baseline - Tl_baseline);
for i = 1:sz

%get the current velocity, theta, salt
Ubl = cell2mat(Ubl_scenarios(i));
Vbl = cell2mat(Vbl_scenarios(i));
Tbl = cell2mat(Tbl_scenarios(i));
Sbl = cell2mat(Sbl_scenarios(i));
Tl  = T0 + lambda*topo - gamma*Sbl; %liquidus temperature in BL
Ubaro = cell2mat(uvel_baro_scenarios(i));
Vbaro = cell2mat(vvel_baro_scenarios(i));


%compute relative melt contributions
UdT = sqrt(Ubl.^2 + Vbl.^2) .* (Tbl - Tl);
UdT_noVel =  sqrt(Ubl_baseline.^2 + Vbl_baseline.^2) .* (Tbl - Tl);
UdT_noTemp = sqrt(Ubl.^2 + Vbl.^2) .* (Tbl_baseline - Tl_baseline);
UdT_barotropic = sqrt(Ubaro.^2 + Vbaro.^2) .* (Tbl_baseline - Tl_baseline);

relmelt(i) = nanmean(UdT(idx)) / nanmean(UdT_baseline(idx)); 
relmelt_noTemp(i) = nanmean(UdT_noTemp(idx)) / nanmean(UdT_baseline(idx));
relmelt_noVel(i) = nanmean(UdT_noVel(idx)) / nanmean(UdT_baseline(idx));
relmelt_baro(i)  = nanmean(UdT_barotropic(idx)) / nanmean(UdT_baro_baseline(idx));

end %end loop over runs
plot([34,34],  [0.4, 1.8], 'k--', 'linewidth', 1.5, 'handlevisibility', 'off'); %plot the location of top of ridge
plot(84 - extent, relmelt, 'o-', 'color', plotcolor1, 'markerfacecolor', plotcolor1);
plot(84 - extent, relmelt_noTemp, 'o-', 'color', plotcolor2, 'markerfacecolor', plotcolor2);
plot(84 - extent, relmelt_baro, 'o--', 'color', plotcolor2, 'markerfacecolor', plotcolor2);
plot(84 - extent, relmelt_noVel, 'o-', 'color', plotcolor3, 'markerfacecolor', plotcolor3);

%tidy
ylim([0.3, 1.8])
xlabel('$l_c$ (km)', 'Interpreter', 'Latex', 'FontSize', ax_fontsize);
ylabel('Relative change', 'Interpreter', 'latex', 'FontSize', ax_fontsize)
xlim([0, 84 - 40]);
legend({"$\mathcal{M}$", "$U_e$", "$\bar{U}_{e}$", "$\Delta T_e$"}, 'location', 'southwest','interpreter', 'latex', 'FontSize', ax_fontsize)

txb = text(-7,1.95, '(b)', 'FontSize', ax_fontsize, 'Interpreter', 'latex');

%add second axis with yf on
axb2 = axes; 
axb2.Position = axb.Position; 
axb2.Color = 'none'; 
axb2.XAxisLocation = 'top';
axb2.XTick = [0,10,20,30,40]/44;
axb2.YTick = [];
axb2.XTickLabel = {"84", "74", "64", "54", "44"};
axb2.XLabel.String = '$y_f$ (km)';
axb2.XLabel.Interpreter= 'latex';
axb2.XLabel.FontSize = ax_fontsize;
%
% Save figure
%
if save_flag
saveas(gcf, '../../PIG-calving-melt-response-tex/figures/figureD_baseline_decomposition.png')
end

% % if you'd like to plot f/h for each of the situations
%figure(1); clf; hold on; box on;
%A = parula(sz+1);
%for i =1:sz
%topo = cell2mat(topo_scenarios(i));
%h = topo-bathy;
%h_lat = h(60,:);
%f_on_h = 2*7.292*1e-5 * sind(-75)./h_lat;
%plot(128 - Y/1e3, f_on_h, 'color', A(i,:));
%legendinfo{i} = ['$\ell_c = ' num2str(84 - extent(i)) '$ km'];
%xlabel('y');
%end
%legend(legendinfo, 'interpreter', 'latex', 'location', 'southwest');
%txt = text(128, -1e-6, 'grounding line', 'rotation', 90);
