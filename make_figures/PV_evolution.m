% Produce plots of the evolution of (quantities related to) PV with calving


% NB: Many of the data files referred to in this script are too large to be hosted online. These files are hostedinternally as BAS.
% Please email Alex Bradley (aleey@bas.ac.uk) to obtain a copy.
%Alex Bradley (aleey@bas.ac.uk) 27/05/2021. MIT license.

addpath('plot_tools');
%
% Flags
%
gendata = 1; %specify whether to pass through the generate data loop
save_flag = 0;


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
WCT_scenarios    = cell(1,sz); %water column thickness
stream_scenarios = cell(1,sz); %barotropic streamfunction
zeta_scenarios   = cell(1,sz); %relative vorticity
BPV_scenarios    = cell(1,sz); %barotropic potential vorticies
WCTterm_scenarios= cell(1,sz); %= u . grad(f/H)
RelViscterm_scenarios = cell(1,sz); %=u . grad(zeta/H)
u_dot_grad_BPV_scenarios = cell(1,sz); %= u.grad(BPV)

for i = 1:sz
%bathy
fid = fopen(bathy_path);
bathy = fread(fid, 'real*8', 'b');
bathy = reshape(bathy, [nx, ny]);
%bathy(bathy == 0) = nan;

%draft
topo_fname=  ['shelfice_topo_H' num2str(H) '_W' num2str(W) '_extent' num2str(extent(i)) 'km.bin'];
topo_fid = fopen(strcat(topodir, '/',topo_fname));
topo = fread(topo_fid, 'real*8', 'b');
topo = reshape(topo, [nx, ny]);

%water column thickness
wct = topo - bathy;
WCT_scenarios{i} = wct;

%Velocities
UVEL_fname = strcat(rootdir, run_nos(i), '/run/stateUvel.nc');
UVEL = ncread(UVEL_fname, 'UVEL', [1,1,1,ntout1], [Inf, Inf, Inf,  1+ntout2 - ntout1]);
UVEL = mean(UVEL, 4);
VVEL_fname = strcat(rootdir, run_nos(i), '/run/stateVvel.nc');
VVEL = ncread(VVEL_fname, 'VVEL', [1,1,1,ntout1], [Inf, Inf, Inf,  1+ntout2 - ntout1]);
VVEL = mean(VVEL, 4);

%% Processed quantities
% BSF
vvel = squeeze(sum(VVEL, 3)) * dz; %units m^2 /s
stream=zeros(size(vvel));
stream(nx,:)=vvel(nx,:)*dx;
for p=nx-1:-1:1
 stream(p,:)=stream(p+1,:) + vvel(p,:)*dx;
end
stream = stream/1e6; %convert to sv
streamsm = smooth2a(stream, 2,2);
stream_scenarios{i} = streamsm;


% Barotropic velocities
vvel = squeeze(sum(VVEL, 3)) * dz;
uvel = squeeze(sum(UVEL, 3)) * dz;
vvel = vvel ./ wct;
uvel = uvel ./ wct;

%Put velocities onto cell centres
uvelC = zeros(nx,ny);
uvelC(1:end-1,:) = (uvel(1:end-1,:) + uvel(2:end,:))/2;
uvelC(end,:)     = (-uvel(end-1,:) + 3*uvel(end,:))/2;
vvelC = zeros(nx,ny);
vvelC(:,1:end-1) = (vvel(:,1:end-1) + vvel(:,2:end))/2;
vvelC(:,end)     = (-vvel(:,end-1) + 3*vvel(:,end))/ 2;
vvel = vvelC;
uvel = uvelC;

%velocity gradients
dvdx = ddx(vvel, dx);
dudy = ddy(uvel, dy);
%smooth dudy
dudyS = dudy;
for p = 1:nx; dudyS(p,:) = smooth(dudyS(p,:)); end
for q = 1:ny; dudyS(:,q) = smooth(dudyS(:,q)); end
dudy = dudyS;

% PV
zeta = dvdx - dudy;
zeta_scenarios{i} = zeta;
f = 2*7.2921*1e-5*sind(-75);
BPV = (f + zeta) ./ wct;
BPV_scenarios{i} = BPV;


% Barotropic vorticity gradients
dBPVdx = ddx(BPV,dx);
dBPVdy = ddy(BPV,dy);
u_dot_grad_BPV = uvel .* dBPVdx + vvel .* dBPVdy;
u_dot_grad_BPV_scenarios{i} = u_dot_grad_BPV;

% u . grad(f/H)
fOnH = f ./wct;
dfOnHdx = ddx(fOnH, dx);
dfOnHdy = ddy(fOnH, dy);
u_dot_grad_fOnH = uvel .* dfOnHdx + vvel .* dfOnHdy;
WCTterm_scenarios{i} = u_dot_grad_fOnH;

% u. grad(zeta/H)
zetaOnH = zeta ./ wct;
dzetaOnHdx = ddx(zetaOnH, dx);
dzetaOnHdy = ddy(zetaOnH, dy);
u_dot_grad_zetaOnH = uvel .* dzetaOnHdx + vvel .* dzetaOnHdy;
RelViscterm_scenarios{i} = u_dot_grad_zetaOnH;

end %end loop over runs
end %end gendata loop

%%%%%%%%%%%%%%%%%%%%% Plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
width = 0.19;
ncols = 5;
colgap = 0.01;
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
offset = 0.02; %colbar offset from edge
colbar_xpos = startx:(colgap+width):(startx + (ncols - 1)*(colgap + width)) + offset;
colbar_width = width - 2*offset; cbar_fontsize = 12;

%
% Column 1: u.grad(BVP)
%
satval = 0.5*1e-11;
for p = 1:sz
ax(1,p) = subplot('Position', squeeze(positions(:,1,sz+1-p)));
u_dot_grad_BPV = cell2mat(u_dot_grad_BPV_scenarios(p));
u_dot_grad_BPVS = saturate(u_dot_grad_BPV, 1.05*satval, -satval);
contourf(-Y/1e3, X/1e3,1e11 * u_dot_grad_BPVS, 50, 'linestyle', 'none');
colormap(ax(1,p), redblue);
xticks([]);
yticks([]);

if p == 1
        a = colorbar;
        a.Location = 'northoutside';
        a.Position(1)=colbar_xpos(1);
        a.Position(2)=colbar_ypos;
	a.Position(3)=colbar_width;
        a.Label.String = '$u . \nabla [(f + \zeta)/H]$';
        a.FontSize = 10;
        a.Label.FontSize = cbar_fontsize;
        a.Label.Interpreter = 'latex';
end

%add bsf
stream = cell2mat(stream_scenarios(p));
streamsm = smooth2a(stream, 2,2);
axnew = axes;
axnew.Position = ax(1,p).Position;
[C,h] =contour(-Y/1e3,X/1e3, streamsm, [-0.7, -0.5, -0.3, -0.1, 0], 'k');
xticks([]);
yticks([]);
set(axnew, 'color', 'none')

end %end loop over plots

% 
% Column 2: u . grad(f/H) [water column thickness term]
%

for p = 1:sz
ax(2,p) = subplot('Position', squeeze(positions(:,2,sz+1-p)));
WCTterm = cell2mat(WCTterm_scenarios(p));
WCTterm = saturate(WCTterm, 1.05*satval, -satval);
contourf(-Y/1e3, X/1e3,1e11 * WCTterm, 50, 'linestyle', 'none');
colormap(ax(2,p), redblue);
xticks([]);
yticks([]);

if p == 1
        b = colorbar;
        b.Location = 'northoutside';
        b.Position(1)=colbar_xpos(2);
        b.Position(2)=colbar_ypos;
	b.Position(3)=colbar_width;
        b.Label.String = '$u . \nabla [f/H]$';
        b.FontSize = 10;
        b.Label.FontSize = cbar_fontsize;
        b.Label.Interpreter = 'latex';
end

%add bsf
stream = cell2mat(stream_scenarios(p));
streamsm = smooth2a(stream, 2,2);
axnew = axes;
axnew.Position = ax(2,p).Position;
[C,h] =contour(-Y/1e3,X/1e3, streamsm, [-0.7, -0.5, -0.3, -0.1, 0], 'k');
xticks([]);
yticks([]);
set(axnew, 'color', 'none')

end %end loop over plots in column 2

% 
% Column 3: u . grad(zeta/H) [relative vorticity term]
%
for p = 1:sz
ax(3,p) = subplot('Position', squeeze(positions(:,3,sz+1-p)));
RelViscterm = cell2mat(RelViscterm_scenarios(p));
RelViscterm = saturate(RelViscterm, 1.05*satval, -satval);
contourf(-Y/1e3, X/1e3,1e11 * RelViscterm, 50, 'linestyle', 'none');
colormap(ax(3,p), redblue);
xticks([]);
yticks([]);

if p == 1
        c = colorbar;
        c.Location = 'northoutside';
        c.Position(1)=colbar_xpos(3);
        c.Position(2)=colbar_ypos;
	c.Position(3)=colbar_width;
        c.Label.String = '$u . \nabla [\zeta/H]$';
        c.FontSize = 10;
        c.Label.FontSize = cbar_fontsize;
        c.Label.Interpreter = 'latex';
end

%add bsf
stream = cell2mat(stream_scenarios(p));
streamsm = smooth2a(stream, 2,2);
axnew = axes;
axnew.Position = ax(3,p).Position;
[C,h] =contour(-Y/1e3,X/1e3, streamsm, [-0.7, -0.5, -0.3, -0.1, 0], 'k');
xticks([]);
yticks([]);
set(axnew, 'color', 'none')

end %end loop over plots in column 3

% 
% Column 4: u . zeta [relative vorticity]
%

for p = 1:sz
ax(4,p) = subplot('Position', squeeze(positions(:,4,sz+1-p)));
zeta = cell2mat(zeta_scenarios(p));
zeta = saturate(zeta, 1.05*2*1e-4, -2*1e-4);
contourf(-Y/1e3, X/1e3,zeta, 50, 'linestyle', 'none');
colormap(ax(4,p), redblue);
xticks([]);
yticks([]);

if p == 1
        d = colorbar;
        d.Location = 'northoutside';
        d.Position(1)=colbar_xpos(4);
        d.Position(2)=colbar_ypos;
	d.Position(3)=colbar_width;
        d.Label.String = '$\zeta$';
        d.FontSize = 10;
        d.Label.FontSize = cbar_fontsize;
        d.Label.Interpreter = 'latex';
end

%add bsf
stream = cell2mat(stream_scenarios(p));
streamsm = smooth2a(stream, 2,2);
axnew = axes;
axnew.Position = ax(4,p).Position;
[C,h] =contour(-Y/1e3,X/1e3, streamsm, [-0.7, -0.5, -0.3, -0.1, 0], 'k');
xticks([]);
yticks([]);
set(axnew, 'color', 'none')

end %end loop over plots in column 4


% 
% Column 5: BPV [barotropic potential vorticity]
%

for p = 1:sz
ax(5,p) = subplot('Position', squeeze(positions(:,5,sz+1-p)));
BPV = cell2mat(BPV_scenarios(p));
BPV = saturate(BPV, -2*1e-7, -16*1e-7);
contourf(-Y/1e3, X/1e3,1e11 * BPV, 50, 'linestyle', 'none');
colormap(ax(5,p), parula);
xticks([]);
yticks([]);

if p == 1
        e = colorbar;
        e.Location = 'northoutside';
        e.Position(1)=colbar_xpos(5);
        e.Position(2)=colbar_ypos;
	e.Position(3)=colbar_width;
        e.Label.String = '$(f + \zeta)/H$';
        e.FontSize = 10;
        e.Label.FontSize = cbar_fontsize;
        e.Label.Interpreter = 'latex';
end

%add bsf
stream = cell2mat(stream_scenarios(p));
streamsm = smooth2a(stream, 2,2);
axnew = axes;
axnew.Position = ax(5,p).Position;
[C,h] =contour(-Y/1e3,X/1e3, streamsm, [-0.7, -0.5, -0.3, -0.1, 0], 'k');
xticks([]);
yticks([]);
set(axnew, 'color', 'none')

end %end loop over plots in column 5

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot max circulation as a function of d/dy(1/H)
max_circ = zeros(1,sz);
topo_strength = zeros(1,sz);
idx = (YY < 50*1e3);
A = jet(sz);
figure(2); clf; hold on; box on
for i = 1:sz
circ = cell2mat(stream_scenarios(i));
max_circ(1,i) = min(min(circ(idx)));

invH = 1./cell2mat(WCT_scenarios(i));
ddy_invH = ddy(invH, dy);
indx = ~isnan(ddy_invH) & ~isinf(ddy_invH);
topo_strength(1,i) = max(max(abs(ddy_invH(indx))));
plot(topo_strength(1,i), max_circ(1,i), 'o', 'MarkerFaceColor', A(i,:), 'MarkerEdgeColor', A(i,:));
end 


