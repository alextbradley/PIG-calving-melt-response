%Make figure 14 of the PIG melt response paper: Evolution of cross sectional temperature for both north and south boxes.



% NB: Many of the data files referred to in this script are too large to be hosted online. These files are hosted internally as BAS.
% Please email Alex Bradley (aleey@bas.ac.uk) to obtain a copy.
%Alex Bradley (aleey@bas.ac.uk) 27/05/2021. MIT license.

%
% Flags
%
gendata = 1; %specify whether to pass through the generate data loop
save_flag = 0;
addpath('plot_tools');
plot_defaults

%
% Data locations
%
rootdir = '/data/oceans_output/shelf/aleey/mitgcm/rPIG_'; %output data NOT in github repo (contact for copy)
topodir = '../gendata_realistic/topo_files/';
bathypath = '../gendata_realistic/bathy_files/bathymetry.shice';

%
% Preliminaries
%
nx = 360;
ny = 320;
nz = 110;
dx = 400;
dy = 400;
dz = 10;
X = ncread(strcat(rootdir, "078", '/run/state2D.nc'), 'LONGITUDE');
Y = ncread(strcat(rootdir, "078", '/run/state2D.nc'), 'LATITUDE'); %stereographic co-ords
[XX,YY] = meshgrid(X,Y);
x = dx:dx:nx*dx;
y = dy:dy:ny*dy; %stereographic co-ords with zero origin
Z = 0:dz:(nz-1)*dz;
[xx,yy] = meshgrid(x,y);


%load bathy
bathyfid = fopen(bathypath);
bathy = fread(bathyfid, 'real*8', 'b');
bathy = reshape(bathy, [nx,ny]);
bathy = double(bathy);
bathy(bathy ==0)  = nan;
figure(1); clf;hold on
 contourf(XX, YY, bathy', 30, 'linestyle','none')

% plot cross sections
fid = fopen("../gendata_realistic/topo_files/shelfice_topo_scn1.shice");
topo = fread(fid, 'real*8', 'b');
topo = reshape(topo, [nx,ny]);
contour(XX,YY, topo', [0,0], 'k', 'linewidth', 2)
realistic_inner_cavity_definition; %bring inner cavity definition into scope (a1,b1,a2,b2)
in1 = inpolygon(XX',YY', a1,b1);
in2 = inpolygon(XX',YY', a2,b2);
idx1 = (topo < 0) & in1;
idx2 = (topo < 0) & in2;
A = zeros(nx,ny);
A(idx1) = 1;
[cin1,~] = contour(X,Y, A',[1,1], 'c');
cin1 = cin1(:, (cin1(1,:)~=1)); %remove levels
A(idx2) = 1;
A(~idx2) = 0;
[cin2,~] = contour(X,Y, A',[1,1], 'm');
cin2 = cin2(:, (cin2(1,:)~=1)); %remove levels

%centreline
np = 100;
xidx = round(linspace(102,334,np));
yidx = round(linspace(160,51, np));
plot(X(xidx), Y(yidx), 'k--');

%north line
xidxN = round(linspace(102,334,np));
yidxN = round(linspace(170,81, np));
s_transect{1} = sqrt((X(xidxN) - X(xidxN(1))).^2 + (Y(yidxN)-Y(yidxN(1))).^2);
pn = plot(X(xidxN), Y(yidxN), 'c--');

%south line
xidxS = round(linspace(102,334,np));
yidxS = round(linspace(140,31, np));
s_transect{2} = sqrt((X(xidxS) - X(xidxS(1))).^2 + (Y(yidxS)-Y(yidxS(1))).^2); %distance along transect
ps = plot(X(xidxS), Y(yidxS), 'm--');

if gendata
run_nos = ["078", "082", "083", "084", "085", "086"]; ntout1 = 11; ntout2 = 12;
%run_nos = ["141", "142", "143", "144", "145", "146"]; ntout1 = 6; ntout2 = 7;
sz = length(run_nos);

%setup storage
theta_scenarios = cell(sz,2);
topo_scenarios = cell(1,sz);
topo_transect_scenarios = cell(sz,2);
bathy_transect_scenarios = cell(sz,2);

%loop over runs
for i = 1:length(run_nos)
%get topo
topo_fname=  ['shelfice_topo_scn', num2str(i), '.shice'];
topo_fid = fopen(strcat(topodir, '/',topo_fname));
topo = fread(topo_fid, 'real*8', 'b');
topo = reshape(topo, [nx,ny]);
topo_scenarios{i} = topo;

%get theta
theta = ncread(strcat(rootdir, run_nos(i), '/run/stateTheta.nc'),'THETA',  [1,1,1,ntout1], [Inf, Inf, Inf, 1 + ntout2 - ntout1]);
theta = mean(theta, 4);

%put it on the transect 
theta_transectN = nan(np,nz);
theta_transectS = nan(np,nz);
topo_transectN  = nan(np,1);
topo_transectS  = nan(np,1);
bathy_transectN  = nan(np,1);
bathy_transectS  = nan(np,1);

for p = 1:np
for q = 1:nz
if bathy(xidxN(p), yidxN(p)) <= -Z(q) && topo(xidxN(p), yidxN(p)) >= -Z(q)
theta_transectN(p,q) = theta(xidxN(p), yidxN(p),q);
end
if bathy(xidxS(p), yidxS(p)) <= -Z(q) && topo(xidxS(p), yidxS(p)) >= -Z(q)
theta_transectS(p,q) = theta(xidxS(p), yidxS(p),q);
end

%get the topo and bathy along transect
%if topo(xidxN(p), yidxN(p)) < 0
topo_transectN(p) = topo(xidxN(p), yidxN(p));
%end
if bathy(xidxN(p), yidxN(p)) < 0
bathy_transectN(p) = bathy(xidxN(p), yidxN(p));
end
%if topo(xidxS(p), yidxS(p)) < 0
topo_transectS(p) = topo(xidxS(p), yidxS(p));
%end
if bathy(xidxS(p), yidxS(p)) < 0
bathy_transectS(p) = bathy(xidxS(p), yidxS(p));
end

end %end loop over z levels
end %end loop over transect pts

theta_scenarios{i,1} = theta_transectN;
theta_scenarios{i,2} = theta_transectS;
topo_transect_scenarios{i,1} = topo_transectN;
topo_transect_scenarios{i,2} = topo_transectS;
bathy_transect_scenarios{i,1} = bathy_transectN;
bathy_transect_scenarios{i,2} = bathy_transectS;
end %end loop over scenarios (i)
end %end gendata loop

%
% Plot setup
%
ncols = 2;
colgap = 0.1;
width = 0.35;
startx = (1 -width*ncols - (ncols-1)*colgap)/2;
starty = 0.07;
hgap = 0.02;
height = 1/(sz+1) - starty/sz  - hgap- 0.02;
positions = zeros(4, ncols, sz);

for p = 1:sz
for q = 1:ncols
positions(:,q,p) = [startx + (q-1)*colgap + (q-1)*width, starty + (p-1)*hgap +  (p-1)*height, width, height];
end
end
positions(2,:,end) = positions(2,:,end) + 0.1;

%
% Make plot
%
figure(2); clf;
for r = 1:2 %north and south boxes
for i = 1:sz
ax(i,r) = subplot('Position', squeeze(positions(:,r,sz+1 - i))); hold on

%if i == 1

%get the temp
theta_diff = cell2mat(theta_scenarios(i,r));
theta_diff = saturate(theta_diff, 1.2, -1.2);
colmap = parula(100);
%else
%theta_diff = cell2mat(theta_scenarios(i,r)) - cell2mat(theta_scenarios(i-1,r));
%theta_diff = saturate(theta_diff, 0.166, -0.15);
%colmap = redblue(100);
%end

%plot the temp
S = cell2mat(s_transect(r));
contourf(S/1e3,-Z,theta_diff', 20, 'linestyle', 'none')
if i == 1
c = colorbar;
c.Location = 'south';
c.Position(2) =  ax(1).Position(2) + ax(1).Position(4) + 0.04;
c.Label.String = 'Potential temp. (${}^\circ$C)';
c.Label.Interpreter = 'latex';
c.Label.FontSize = 10;
c.Label.Position(2) = 4.2;
elseif i == 2
d = colorbar;
d.Location = 'north';
d.Position(2) =  ax(2).Position(2) + ax(2).Position(4) + 0.04;
d.Label.String = 'Potential temp. anomaly (${}^\circ$C)';
d.Label.Interpreter = 'latex';
d.Label.FontSize = 10;
d.Label.Position(2) = 4.2;
end

%add the bathy and topo
bathyline = cell2mat(bathy_transect_scenarios(i,r));
topoline = cell2mat(topo_transect_scenarios(i,r));
plot(S/1e3, bathyline, 'k', 'linewidth', 1.5);
plot(S/1e3, topoline, 'k', 'linewidth', 1.5);

%tidy plot
xlim([30, max(S/1e3)])
colormap(ax(i,r), colmap)
if r == 1
ylim([-1000, -300])
else
ylim([-1100, -400])
end

if i == 1
ax(i,r).YLabel.String = 'depth (m)';
ax(i,r).YLabel.Interpreter= 'latex';
ax(i,r).YLabel.FontSize = 11;
end


if i == 6
ax(i,r).XLabel.String = 'Distance along transect (km)';
ax(i,r).XLabel.Interpreter = 'latex';
ax(i,r).XLabel.FontSize = 11;
else
ax(i,r).XTickLabel = cell(1,length(ax(i).XTickLabel));
end 
box on
end
end

%
% Macro tidy up
%
fig = gcf;
fig.Position(3:4) = [940, 530];

%
% save
%
set(gcf, 'color', 'w')
if save_flag
%saveas(gcf, "plots/figure14.eps", "epsc")
saveas(gcf, "plots/figure14.png")
end
