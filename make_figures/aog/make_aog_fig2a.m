% 29/08/22
%
% Make Figure 2a in the AoG paper: plot of the % change in inner cavity melt rate as a function of calved length (figure 14c of the JGRO paper)
% NB: must run make_aog_fig1a first to get the calved length.

%
% Data locations
%
rootdir = '/data/oceans_output/shelf/aleey/mitgcm/rPIG/rPIG_'; %output data NOT in github repo (contact for copy)
topodir = '../../gendata_realistic/topo_files/';
bathypath = '../../gendata_realistic/bathy_files/bathymetry.shice';

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

% misc
ntout1 = 10;
ntout2 = 12; %define time period to average over
secs_per_year = 365*24*60^2;
density_ice = 918.0;

% get data
run_nos = ["078", "082", "083", "084", "085", "086"];
sz = length(run_nos);
melt_scenarios = cell(1,sz);
topo_scenarios = cell(1,sz);

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

end

% inner cavity
cav_def; %bring inner cavity definition into scope (a1,b1,a2,b2)
in1 = inpolygon(XX',YY', a1,b1);
in2 = inpolygon(XX',YY', a2,b2);
topo = cell2mat(topo_scenarios(1));
idx1 = (topo < 0) & in1;
idx2 = (topo < 0) & in2;
idx = ( idx1 | idx2);

for i = 1:sz
melt = cell2mat(melt_scenarios(i));
ave_melt1(i) = mean(melt(idx1)); %average melt in the north part
ave_melt2(i) = mean(melt(idx2)); %average melt in the southern part
end
s1 = sum(sum(idx1)); s2 = sum(sum(idx2)); st = s1 + s2;
pc1 = (ave_melt1 - ave_melt1(1)) / ave_melt1(1) * 100; %percetnage change in box north
pc2 = (ave_melt2 - ave_melt2(1)) / ave_melt2(1) * 100; %south box
pct = pc1 * (s1 / (s1 + s2)) + pc2 * (s2 / (s1 + s2)); %weighted percentage change (mean in inner cavity)

% make plot
figure(1); clf; hold on; box on
plot(snap_distance/1e3, pct, 'ko--','LineWidth', 1, 'MarkerSize', 4, 'MarkerFaceColor', 'k');
grid on
ylim([0, 10])

xlabel('calved length (km)');
ylabel(sprintf('%% change in melt rate'))

