% Make figure 1. Four panels:
% (a) PIG bathymetry with 2009, 2020 ice fronts and approx ridge crest 
% (b) PIG water column thickness with 2009, 2020 ice fronts and approx ridge crest
% (c) ice topo and bathymetry along the ridge crest
% (d) gap along the contour
%
% Note that this figure appears as (matlab) figure 2. Code also produces contour plots of bathymetry and fronts, whose data are used in the main figure.

%Alex Bradley (aleey@bas.ac.uk) 27/05/2021. MIT license.

%
% Plot info
%
saveflag = 0;
addpath('plot_tools');
plot_defaults           %set plot defaults
positions = [0.08,0.37,0.40,0.53; 
	     0.58,0.37,0.40,0.53;	
             0.08, 0.1, 0.40, 0.25;
             0.58, 0.1, 0.40, 0.25]; %subplot positions

close all

%
% Grid
%
nx = 360;
ny = 320;
dx = 400;
dy = 400;
x = dx:dx:nx*dx;
y = dy:dy:ny*dy;
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


% 
% Bathy and topos
%
fid = fopen("../gendata_realistic/bathy_files/bathymetry.shice");
bathy = fread(fid, 'real*8', 'b');
bathy = reshape(bathy, [nx,ny]);
bathyng = bathy; bathyng(bathy ==0) = nan;

fid = fopen("../gendata_realistic/topo_files/shelfice_topo_scn1.shice");
topo2009 = fread(fid, 'real*8', 'b');
topo2009 = reshape(topo2009, [nx,ny]);

fid = fopen("../gendata_realistic/topo_files/shelfice_topo_scn2.shice");
topo2020 = fread(fid, 'real*8', 'b');
topo2020 = reshape(topo2020, [nx,ny]);

% 
% engineer topos so that gl doesn't show up at zero contour
%
for p = 2:nx-1
    for q = 2:ny-1
        if any( bathy(p,q+1) == 0 || bathy(p, q-1) == 0 || bathy(p+1,q) == 0 || bathy(p-1,q) == 0)
            topo2009(p,q) = nan;
            topo2020(p,q) = nan;
        end
    end
end

%
% plot bathy and topos and get data
%
figure(1); clf; hold on
gap = topo2009 - bathy;
%contourf(x,y,bathyng', 20, 'linestyle', 'none')
contourf(x,y,gap', 50, 'linestyle', 'none')

[cgl,~] = contour(x,y,bathy', [0,0], 'k');
[c2009,~] = contour(x,y,topo2009', [0,0], 'r');
[c2020,~] = contour(x,y,topo2020', [0,0], 'b');


%
%transect
%
%yidx = 5:140;
%xidx = 256*ones(1,length(yidx));

np = 50; %nmber of pts on each half of transect
cx = 245;
cy = 73;
xidx1 = round(linspace(257, cx, np));
yidx1 = round(linspace(5, cy, length(xidx1)));
%p1 =  plot(x(xidx1), y(yidx1), 'ro');

xidx2 = round(linspace(cx, 259, np));
yidx2 = round(linspace(cy, 140, length(xidx2)));
%p2 =  plot(x(xidx2), y(yidx2), 'bo');

xidx = [xidx1, xidx2];
yidx = [yidx1, yidx2];


sline =  sqrt((x(xidx) - x(xidx(1))).^2 + (y(yidx) - y(yidx(1))).^2); %arclength along line
xline = x(xidx);
yline = y(yidx);
plot(xline, yline, 'k--', 'linewidth', 2);


%
%get bathy and topo along line
%
bathyline = nan(1,length(yidx));
topoline  = nan(1,length(yidx));
for i = 1:length(bathyline)
if bathy(xidx(i), yidx(i)) ~=0
bathyline(i) = bathy(xidx(i), yidx(i));
end
if topo2009(xidx(i), yidx(i)) ~=0
topoline(i)  = topo2009(xidx(i), yidx(i));
end
end
idx = ~isnan(bathyline);
sline = sline(idx);
sline = sline - sline(1);
bathyline = bathyline(idx);
topoline = topoline(idx);

%
% Open image
% 
t = Tiff('PIG-S2-NovDec2020.tif', 'r'); %!! not in git repo!!
imageData = read(t);
figure(2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Make Plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
% Plot 1: bathymetry
%
ax(1) = subplot('Position',positions(1,:)); imshow(imageData);

% add bathymetry colours and some contours
hold on
contourf(xi, yi, bathyng', 30, 'linestyle', 'none'); 
c = colorbar('Location', 'northoutside');
c.Position(2) = 0.91;
c.Position(4) = 0.02;
c.Label.String = 'seabed depth (m)';
c.Label.VerticalAlignment = 'bottom';
c.Label.FontSize = 12;
c.Label.Interpreter = 'latex';
contour(xi, yi, bathyng', [-1000,-900, -800,-700, -600, -500], 'linecolor', [211,211,211]/255);


%add grounding line
[cgl,~] = contour(xi,yi,bathy', [0,0], 'k', 'linewidth',2 );

%add ice fronts
[c2009,h2009] = contour(xi,yi,topo2009', [0,0], 'linecolor', plotcolor3, 'linewidth', 2);
[c2020,h2020] = contour(xi,yi,topo2020', [0,0], 'linecolor', plotcolor1, 'linewidth', 2);

% add cross section
cline = [xline; yline];
cline_img = model2image(cline);
plot(cline_img(1,:), cline_img(2,:), 'k--');
ptA = text(ax(1), 9700,8200, 'A', 'FontSize',12, 'FontWeight', 'bold');
ptB = text(ax(1), 5300,8900, 'B', 'FontSize',12, 'FontWeight', 'bold');

% Label fronts
f2009 = text(ax(1), 4500,13500, '2009', 'FontSize',12, 'color', plotcolor3);
f2020 = text(ax(1), 8700,10200, '2020', 'FontSize',12, 'color', plotcolor1);

%tidy
ylim([3000, 14847])
xlim([2000,12000])
plot(ax(1), [11000,11000], [5000,4000], 'k', 'linewidth', 3)
scalebar = text(ax(1),10950, 6600, '10 km', 'Interpreter', 'latex', 'FontSize', 12);
camroll(-90);
%txa = text(ax(1), 250, 16000, '(a)', 'FontSize', 12, 'interpreter', 'latex');

%
% Plot 2: Water column thickness
%
ax(2) = subplot('Position',positions(2,:)); imshow(imageData);

% add bathymetry and fronta
hold on
contourf(xi, yi, gap', 30, 'linestyle', 'none'); 
d = colorbar('Location', 'northoutside');
d.Position(2) = c.Position(2);
d.Position(4) = c.Position(4);
d.Label.String = 'water column thickness (m)';
d.Label.VerticalAlignment = 'bottom';
d.Label.FontSize = 12;
d.Label.Interpreter = 'latex';

[cgl,~] = contour(xi,yi,bathy', [0,0], 'k', 'linewidth',2 );
[c2009,h2009] = contour(xi,yi,topo2009', [0,0], 'linecolor', plotcolor3, 'linewidth', 2);
[c2020,h2020] = contour(xi,yi,topo2020', [0,0], 'linecolor', plotcolor1, 'linewidth', 2);


contour(xi, yi, gap', [125,125], 'm'); 

% add cross section
cline = [xline; yline];
cline_img = model2image(cline);
plot(cline_img(1,:), cline_img(2,:), 'k--');
ptA = text(ax(2), 9700,8200, 'A', 'FontSize',12, 'FontWeight', 'bold');
ptB = text(ax(2), 5300,8900, 'B', 'FontSize',12, 'FontWeight', 'bold');

% Label fronts
%f2009 = text(ax(2), 4500,13500, '2009', 'FontSize',12, 'color', plotcolor3);
%f2020 = text(ax(2), 8700,10200, '2020', 'FontSize',12, 'color', plotcolor1);

%tidy
ylim([3000, 14847])
xlim([2000,12000])
%plot(ax(2), [11000,11000], [5000,4000], 'k', 'linewidth', 3)
%scalebar = text(ax(1),10950, 6600, '10 km', 'Interpreter', 'latex', 'FontSize', 12);
camroll(-90);

%
% Plot along the cross section
%
ax(3) = subplot('Position', positions(3,:)); grid on 
plot(sline/1e3, bathyline, 'color', 'k', 'linewidth', 2);
hold on
plot(sline/1e3, topoline,'--', 'color', 'k', 'linewidth', 2);
xlim([0, 45])
ptAA = text(ax(3), 1,-750, 'A', 'FontSize', 12, 'FontWeight', 'bold');
ptBB = text(ax(3), max(ax(3).XLim)-3 ,-750, 'B', 'FontSize', 12, 'FontWeight', 'bold');

%ax(3).XTickLabels = cell(length(ax(3).XTick), 1);

ax(4) = subplot('Position', positions(4,:)); grid on
plot(sline/1e3,topoline - bathyline, 'color', 'k', 'linewidth', 2);
xlim([0, 45])
ptAA = text(ax(4), 1,300, 'A', 'FontSize', 12, 'FontWeight', 'bold');
ptBB = text(ax(4), max(ax(4).XLim)-3 ,300, 'B', 'FontSize', 12, 'FontWeight', 'bold');


%add subplot labels
laba = text(ax(1), 200, 17000, '(a)', 'FontSize', 12,  'Interpreter', 'latex');
labb = text(ax(2), 200, 17000, '(b)', 'FontSize', 12,  'Interpreter', 'latex');
labc = text(ax(3), -8.2, -200, '(c)', 'FontSize', 12, 'Interpreter', 'latex');
labd = text(ax(4), -9, 400, '(d)', 'FontSize', 12, 'Interpreter', 'latex');

%
% tidy everything
%
ax(4).XLabel.String = 'Distance along transect (km)';
ax(3).XLabel.String = 'Distance along transect (km)';
ax(4).XLabel.Interpreter = 'latex';
ax(3).XLabel.Interpreter = 'latex';
ax(4).XLabel.FontSize = 12;
ax(3).XLabel.FontSize = 12;
ax(3).YLabel.String = 'depth (m)';
ax(3).YLabel.FontSize = 12;
ax(3).YLabel.Interpreter = 'latex';
ax(3).YLim = [-850, -200];
ax(3).YTick = [-800, -600, -400,-200];
ax(4).YTick = [0,100,200,300, 400];
ax(4).YLabel.String = 'gap thickness (m)';
ax(4).YLabel.FontSize = 12;
ax(4).YLabel.Interpreter = 'latex';
ax(3).YLabel.FontSize = 12;
ax(3).YLabel.Interpreter = 'latex';

%txb = text(ax(3), -8, 80, '(b)', 'FontSize', 12, 'interpreter', 'latex')
%txc = text(ax(4), -8, 420, '(c)', 'FontSize', 12, 'interpreter', 'latex')

fig = gcf; fig.Position(3:4) = [900, 600];
grid(ax(3), 'on')
grid(ax(4), 'on')

%
ax(3).YLabel.FontSize = 12;
ax(3).YLabel.Interpreter = 'latex';



% Save?
%
if saveflag
saveas(gcf, '../../PIG-calving-melt-response/figures/figureA_maps_and_ridge_gap.png');
end
