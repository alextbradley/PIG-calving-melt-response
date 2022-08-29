% Make figure 1b in the AoG paper: melt rate in the 2009 configuration and cumulative melt rate anomalies

% NB: Many of the data files referred to in this script are too large to be hosted online. These files are hosted internally as BAS.
% Please email Alex Bradley (aleey@bas.ac.uk) to obtain a copy.
%Alex Bradley (aleey@bas.ac.uk) 27/08/2022. MIT license.

%
% Flags
%
gendata = 1; %specify whether to pass through the generate data loop

%
% Preliminaries
%
addpath("../plot_tools");
plot_defaults
label_size = 11;
ax_fontsize = 7.5;
figure(1); clf;
fig = gcf; fig.Position(3:4) = [1085, 540];
%plot colours
background_color = 0.8*[1,1,1]; %ice color
background_color = [1,1,1];
ocean_color = [226,249, 255]/255;
cols = [7,43, 184;
        184, 7, 43]/255; %colours for each different bc
bathycontourcolor = [255,0,255]/255;

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

%put onto the latlon grid
xshift = 870;
yshift = -15000;
lambdashift = -101;
phishift = -0.07;
XX = XX + xshift;
YY = YY + yshift;
[phi, lambda] = polarstereo_inv(XX,YY);
lambda = lambda + lambdashift;
phi = phi + phishift;
lambda = double(lambda);
phi = double(phi);


%parameters
secs_per_year = 365.25*24*60*60;
density_ice = 918.0;


%
% Generate data loop
%
if gendata
run_nos = ["078", "082", "083", "084", "085", "086"]; ntout1 = 10; ntout2 = 12;
%run_nos = ["141", "142", "143", "144", "145", "146"]; ntout1 = 6; ntout2 = 7;
sz = length(run_nos);

%setup storage
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
end

%
% Make the plot
%
width = 0.30;
height = 0.44;
gapx = 0.01;
gapy = 0.01;
ncols = 3;
nrows = 2;
startx = (1 -width*ncols - (ncols-1)*gapx)/2;
starty = 1-(1 -height*nrows - (nrows -1)*gapy)/2;
starty = 0.97;
positions = zeros(4, nrows* ncols);
for i = 1:nrows*ncols
q = 1 + mod(i-1,ncols); %index in x direction
p = ceil(i/ncols); %index in y directio
positions(:,i) = [startx + (q-1)*gapx + (q-1)*width, starty - p*height - (p-1)*gapy, width, height];
end
cmap = lighter_blue_parula(100,0.1);
labels = ["(a)", "(b)", "(c)", "(d)", "(e)", "(f)"];

%loop over runs
for q = 1:sz
icetopo = cell2mat(topo_scenarios(q));
melt = cell2mat(melt_scenarios(q));

%if baseline, store as such. Otherwise, get the anomaly
if q == 1
melt_diff = melt; 
else
melt_diff = melt - cell2mat(melt_scenarios(1)); %
end
melt_diff(icetopo == 0) = nan;
if q == 1
melt_diff_sat = saturate(melt_diff, 126,0);
else
melt_diff_sat = saturate(melt_diff, 22,-20);
end


%make plot
axs(q) = subplot('Position',positions(:,q));
contourf(lambda,phi,melt_diff_sat', 20, 'linestyle', 'none')
hold on
contour(lambda,phi,icetopo', [0,0], 'k', 'linewidth', 1.5)
if q == 1
        colormap(axs(q),cmap);
else
        colormap(axs(q),redblue);
end
xlim([-102.3, -99])
ylim([-75.45,-74.75])
set(axs(q),'Color',background_color)
        
%add the open ocean
bathynoice = bathy;
bathynoice(icetopo < 0) = 0; %remove cavity
bathynoice(bathy == 0) = 0; %remove grounded ice
bathynoice(bathynoice ~= 0) = 1; %make constant
cf = contourf(lambda, phi, bathynoice',[1,1]);
idx = find((cf(1,:) == 1));
c1 = cf(1,2:idx(2)-1);
c2 = cf(2,2:idx(2)-1);
fill(c1,c2,ocean_color, 'linewidth', 1.5, 'edgecolor', 'k');
        
%add 750m bathymetric contour
%contour(lambda, phi, -1e-2*bathy', -1e-2*[-750, -750], 'color',bathycontourcolor, 'linewidth' ,1, 'linestyle', '--'); %he weird "* -1e-2" is to get bring -750 into range of colourar, avoid skewing it
grid on
axs(q).YTick = -75.4:0.1:-74.8;
if mod(q,3) == 1
axs(q).YTickLabel = {"24'", "18'",  "12'", "6'",  "75S", "54'", "48'"};
else
axs(q).YTickLabel = cell(length(axs(q).YTick),1);
end

axs(q).XTick = [-102.0:0.5:-99.2];
if floor((q-1)/3) == 1
axs(q).XTickLabel = {"102W",  "30'", "101W", "30'", "100W", "30'"};
else
axs(q).XTickLabel = cell(length(axs(q).XTick),1);
end

%add the definitions of the two cavity regions
realistic_inner_cavity_definition; %bring inner cavity definition into scope (a1,b1,a2,b2)
[XXns,YYns] = meshgrid(X,Y);

in1 = inpolygon(XXns',YYns', a1,b1);
in2 = inpolygon(XXns',YYns', a2,b2);
idx1 = (icetopo < 0) & in1;
idx2 = (icetopo < 0) & in2;
%if ((q == 1) || (q == 6))
melt_copy = melt_diff_sat;
melt_copy(idx1) = 1;
melt_copy(~idx1) =0;
%contour(lambda,phi,melt_copy',[1,1], '--', 'linewidth', 1.75, 'linecolor', 'c')
melt_copy = melt_diff_sat;
melt_copy(idx2) = 1;
melt_copy(~idx2) =0;
%contour(lambda,phi,melt_copy',[1,1], '--', 'linewidth', 1.75, 'linecolor', 'm')
%end

end %end loop over runs
fig = gcf; fig.Position(3:4) = [790, 500];
set(gcf, 'color', 'w');
for i = 1:6; axs(i).FontSize = 7.5;end

