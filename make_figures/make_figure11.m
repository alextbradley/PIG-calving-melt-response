%Make figure 11 of the PIG calving manuscript: description of calving experiments
% (a) Satellite image of PIG with grounding line, 2009, 2012 ice fronts, synthetic ice fronts and north south cavities labelled. Inset: ice-bathy gap along dashed line in main figure.
% (b) Satellite image of PIG with cavity water column thickness in colours.

% NB: Many of the data files referred to in this script are too large to be hosted online. These files are hosted internally as BAS.
% Please email Alex Bradley (aleey@bas.ac.uk) to obtain a copy.
%Alex Bradley (aleey@bas.ac.uk) 27/05/2021. MIT license.
saveflag = 0; %save toggle
addpath("plot_tools")
plot_defaults
plotcolor2 = [0,1,1]; %overwrite defaults 
plotcolor3 = [1,0,1];
%
% Data info
%
topodir = '../gendata_realistic/topo_files/';
bathypath = '../gendata_realistic/bathy_files/bathymetry.shice';
rootdir = '/data/oceans_output/shelf/aleey/mitgcm/rPIG/rPIG_'; %output data NOT in github repo (contact for copy)
nx = 360;
ny = 320;
dx = 400;
dy = 400;
X = ncread(strcat(rootdir, "078", '/run/state2D.nc'), 'LONGITUDE');
Y = ncread(strcat(rootdir, "078", '/run/state2D.nc'), 'LATITUDE'); %stereographic co-ords
[XX,YY] = meshgrid(X,Y);
x = dx:dx:nx*dx;
y = dy:dy:ny*dy; %stereographic co-ords with zero origin
[xx,yy] = meshgrid(x,y);

%
% get grid in image space
%
load('image_to_model_points.mat', 'ximage', 'yimage', 'xmod', 'ymod');
[model2image, image2model] = genmaps_image2model(ximage, yimage, xmod, ymod); %maps from model to image and vice
yi = zeros(size(yy));
xi = zeros(size(xx));
for i = 1:320
cc = [xx(i,:); yy(i,:)];
ci = model2image(cc);
yi(i,:) = ci(2,:);
xi(i,:) = ci(1,:);
end

%load bathy
bathyfid = fopen(bathypath);
bathy = fread(bathyfid, 'real*8', 'b');
bathy = reshape(bathy, [nx,ny]);
bathy = double(bathy);

%compute distance along contour
figure(1); clf; hold on
%bathy(bathy ==0) = nan;
[cgl,~] = contour(x,y, bathy',[0,0], 'linestyle', 'none');
[c750,~] = contour(x,y, bathy',[-750,-750], 'linestyle', 'none');

A = lines(6);
topo_scenarios = cell(6,1);
for i = 1:6
    topo_fname=  ['shelfice_topo_scn', num2str(i), '.shice'];
    topo_fid = fopen(strcat(topodir, '/',topo_fname));
    topo = fread(topo_fid, 'real*8', 'b');
    topo = reshape(topo, [nx,ny]);
    %engineer topo so that grounding line doesn't show up at zero contour (i.e. only ice front)
    for p = 2:nx-1
        for q = 2:ny-1
            if any( bathy(p,q+1) == 0 || bathy(p, q-1) == 0 || bathy(p+1,q) == 0 || bathy(p-1,q) == 0)
                topo(p,q) = nan;
            end
        end
    end

topo_scenarios{i} = topo;
end

%line cross section definition
%xidx = [102,334];
%yidx = [160,41];
%xline_idx = min(xidx):max(xidx); %x indices of points on the line
%yline_idx = round(diff(yidx)/diff(xidx) * (xline_idx - xidx(end)) + yidx(end));%corresponding y indices
%xline_idx = xline_idx(50:200);
%yline_idx = yline_idx(50:200); %remove some entries
np = 50;
xline_idx = round(linspace(102,334,np));
yline_idx = round(linspace(160,51, np));

sline =  sqrt((x(xline_idx) - x(xline_idx(1))).^2 + (y(yline_idx) - y(yline_idx(1))).^2); %arclength along line
hold on
xline = x(xline_idx);
yline = y(yline_idx);
plot(xline, yline, 'k--')

%compute distance along line
snap_distance = zeros(1,6);
xsnap = zeros(1,6);
ysnap = zeros(1,6);
for i = 1:6
[c,h] = contour(x,y,cell2mat(topo_scenarios(i))', [0,0], 'color', A(i,:));

%store this info
topo_front_data{i} = c;

%remove any lvels
c1 = c(1,:);
c2 = c(2,:);
c1 = c1(c1 ~=0);
c2 = c2(c1 ~=0); %remove level spec

%loop over every co-ordinate in the contour, find the nearest pt in the line, and ge tthe min of all of these
min_idx = 1;
min_dist = 1e10; %large to start with
for j = 1:length(c1)
[val, idx] = min(abs((xline - c1(j)).^2 + (yline - c2(j)).^2));
if val < min_dist
	min_dist = val;
	min_idx = idx;
	
end
end
plot(xline(min_idx), yline(min_idx), 'ro-', 'markerfacecolor','r');
snap_distance(i) = sline(min_idx);
xsnap(i) = xline(min_idx);
ysnap(i) = yline(min_idx);
end
snap_distance = snap_distance - snap_distance(1); %take relative to 2012 topo

%cross section
ycrossidx = 5:140;
xcrossidx = 256*ones(1,length(ycrossidx));
xcross = x(xcrossidx);
ycross = y(ycrossidx);


%get inner cavity contours
realistic_inner_cavity_definition; %bring inner cavity definition into scope (a1,b1,a2,b2)
in1 = inpolygon(XX',YY', a1,b1);
in2 = inpolygon(XX',YY', a2,b2);
idx1 = (topo < 0) & in1;
idx2 = (topo < 0) & in2;
A = zeros(nx,ny);
A(idx1) = 1;
[cin1,~] = contour(x,y, A',[1,1], 'linestyle', 'none');
cin1 = cin1(:, (cin1(1,:)~=1)); %remove levels
A(idx2) = 1;
A(~idx2) = 0;
[cin2,~] = contour(x,y, A',[1,1], 'linestyle', 'none');
cin2 = cin2(:, (cin2(1,:)~=1)); %remove levels

%north and south transect lines
xlineidxN = round(linspace(102,334,np)); %north line x indices
ylineidxN = round(linspace(170,81, np)); %north line y indices
xlineidxS = round(linspace(102,334,np)); % etc
ylineidxS = round(linspace(140,31, np)); % etc
xlineN = x(xlineidxN); 		 	 % north line x co-ordinates
ylineN = y(ylineidxN);		  	 % north line y co-ordinates
xlineS = x(xlineidxS); 		 	 % south line x co-ordinates
ylineS = y(ylineidxS);		  	 % south line y co-ordinates




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% make the plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%open the image
t = Tiff('PIG-S2-NovDec2020.tif', 'r'); %!! not in git repo!!
imageData = read(t);
fig1 = figure(1); clf; hold on
imshow(imageData);
ax = gca;
%add the reference points
load('image_to_model_points.mat', 'ximage', 'yimage', 'xmod', 'ymod');

%add the 2012 grounding line 
[model2image, image2model] = genmaps_image2model(ximage, yimage, xmod, ymod); %maps from model to image and vice versa
c_image_GL = model2image(cgl); %put gl model gl position onto image
scatter(ax, c_image_GL(1,:), c_image_GL(2,:), 4,'k', 'filled')
%contour(xi, yi, bathy', [0,0], 'k', 'linewidth', 2)

%add the fronts
colmap = lines(6);
colmap(1,:) = [1,0,0]; %make first red
colmap(end,:) = [0,1,0]; %last row green
for i = 6:-1:1 %reverse order so 2009 on top
    cfront = cell2mat(topo_front_data(i));
    c_image_front = model2image(cfront);
    scatter(c_image_front(1,:), c_image_front(2,:), 5, colmap(i,:), 'filled');
end

%add the calving front measurement line
cl = [xline; yline];
cl_image = model2image(cl);
plot(cl_image(1,:), cl_image(2,:), 'b--');

%add the north south cross section lines
clN = [xlineN; ylineN];
cl_imageN = model2image(clN);
plot(cl_imageN(1,:), cl_imageN(2,:), '--', 'color', 0.6*[1,1,1]);
clS = [xlineS; ylineS];
cl_imageS = model2image(clS);
plot(cl_imageS(1,:), cl_imageS(2,:), '--','color', 0.6*[1,1,1]);

%add the calving front measurement positions
csnap = [xsnap; ysnap];
csnap_image = model2image(csnap);
plot(csnap_image(1,:), csnap_image(2,:), 'ko', 'markerfacecolor', 'k');

%add the ridge cross section
cx = [xcross; ycross];
cx_image = model2image(cx);
plot(cx_image(1,:), cx_image(2,:), 'k--')

%add the inner cavity definition
cin1_img = model2image(cin1); 
cin2_img = model2image(cin2);
plot(cin1_img(1,:), cin1_img(2,:), '--','color', plotcolor2);
plot(cin2_img(1,:), cin2_img(2,:), '--', 'color',plotcolor3);


%tidy plot
ylim([2194, 14847])
xlim([2000,12000])

%rotate
camroll(-90);

%add the A,B pts and 2009, 2020 front labels
ax1 = gca;
ptAA = text(ax1,9800,8000, 'A', 'FontSize', 12, 'FontWeight', 'bold');
ptBB = text(ax1, 4300,8800, 'B', 'FontSize', 12, 'FontWeight', 'bold');
txtN = text(ax1, 5200, 6400, 'North', 'FontSize', 12, 'Interpreter', 'latex');
txtS = text(ax1, 8700, 6000, 'South', 'FontSize', 12, 'Interpreter', 'latex');
f2009 = text(ax1, 4000,13200, '2009', 'FontSize',12, 'color', colmap(1,:));
f2020 = text(ax1, 7500,11700, '2020', 'FontSize',12, 'color', colmap(2,:));



%add the gap width as inset axes
topo = cell2mat(topo_scenarios(1));
bathyline = nan(1,length(ycrossidx));
topoline  = nan(1,length(ycrossidx));
sline =  sqrt((x(xcrossidx) - x(xcrossidx(1))).^2 + (y(ycrossidx) - y(ycrossidx(1))).^2); %arclength along line

for i = 1:length(bathyline)
if bathy(xcrossidx(i), ycrossidx(i)) ~=0
bathyline(i) = bathy(xcrossidx(i), ycrossidx(i));
end
if topo(xcrossidx(i), ycrossidx(i)) ~=0
topoline(i)  = topo(xcrossidx(i), ycrossidx(i));
end
end
idx = ~isnan(bathyline);
sline = sline(idx);
sline = sline - sline(1);
bathyline = bathyline(idx);
topoline = topoline(idx);
ax2 = axes; hold on; box on
ax2.Position = [0.55, 0.07, 0.3, 0.2];
fill(ax2,[0, 22, 22, 0], [0, 0, 400, 400], 'm','linestyle', 'none', 'FaceAlpha', 0.3);
fill(ax2,[22, 45, 45, 22], [0, 0, 400, 400], 'c','linestyle', 'none', 'FaceAlpha', 0.3);
plot(ax2,sline/1e3,topoline - bathyline, 'color', 'k', 'linewidth', 2);

ax2.XLabel.String = 'distance (km)';
ax2.YLabel.String = 'gap (m)';
ax2.XLabel.Interpreter = 'latex';
ax2.XLabel.FontSize = 12;
ax2.YLabel.Interpreter = 'latex';
ax2.YLabel.FontSize = 12;
ptA2 = text(ax2, 0.5,30, 'A', 'FontSize', 10, 'FontWeight', 'bold');
ptB2 = text(ax2, 42.2,30, 'B', 'FontSize', 10, 'FontWeight', 'bold');

%
% Make a second figure with the 1/h contours
%
fig2 = figure(2);clf; hold on;  fig2.Position = fig1.Position;
imshow(imageData);

%add the inverse water column thickness
fid = fopen("../gendata_realistic/bathy_files/bathymetry.shice");
bathy = fread(fid, 'real*8', 'b');
bathy = reshape(bathy, [nx,ny]);
bathyng = bathy; bathyng(bathy ==0) = nan;

fid = fopen("../gendata_realistic/topo_files/shelfice_topo_scn1.shice");
topo2009 = fread(fid, 'real*8', 'b');
topo2009 = reshape(topo2009, [nx,ny]);
gap = topo2009 - bathy;
gap(bathy == 0) = nan;
invgap = (1./gap);
invgap = saturate(invgap, 0.01,0);
hold on
contourf(xi, yi, invgap', 30, 'linestyle', 'none');
d = colorbar;
d.Label.String = '$1/h$ (m\textsuperscript{-1})';
d.Label.Interpreter = 'latex';
d.Label.FontSize = 12;
d.Position = [0.85, 0.28, 0.03, 0.5];


ax = gca;
%colormap(ax, viridis);
scatter(ax, c_image_GL(1,:), c_image_GL(2,:), 4,'k', 'filled')
%contour(xi, yi, bathy', [0,0], 'k', 'linewidth', 2)

%add the fronts
for i = 6:-1:1
    cfront = cell2mat(topo_front_data(i));
    c_image_front = model2image(cfront);
    scatter(c_image_front(1,:), c_image_front(2,:), 5, colmap(i,:), 'filled');
end
%add the inner cavity definition
cin1_img = model2image(cin1); 
cin2_img = model2image(cin2);
plot(ax,cin1_img(1,:), cin1_img(2,:), '--','color', plotcolor2);
plot(ax,cin2_img(1,:), cin2_img(2,:), '--', 'color',plotcolor3);


%tidy plot
ylim([2194, 14847])
xlim([2000,12000])

%rotate
camroll(-90);


%
% save flag
%
set(fig1, 'color', 'w')
set(fig2, 'color', 'w')
if saveflag
saveas(fig1, "../../PIG-calving-melt-response-tex/figures/figureK_realistic_setup_a.png");
saveas(fig2, "../../PIG-calving-melt-response-tex/figures/figureK_realistic_setup_b.png");
end
