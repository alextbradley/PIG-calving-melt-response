% 25/08/22
%
% Make figure 1a in the annals of glaciology paper. 
% Experimental setup (as figure 11 in the JGR paper)
%


%
% Data info
%
addpath('..')
addpath('../plot_tools')
topodir = '../../gendata_realistic/topo_files/';
bathypath = '../../gendata_realistic/bathy_files/bathymetry.shice';
rootdir = '/data/oceans_output/shelf/aleey/mitgcm/rPIG/rPIG_'; %output data NOT in github repo (contact for copy)
nx = 360;
ny = 320;
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
load('../plot_tools/image_to_model_points.mat', 'ximage', 'yimage', 'xmod', 'ymod');
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

np = 50;
xline_idx = round(linspace(102,334,np));
yline_idx = round(linspace(160,51, np));

sline =  sqrt((x(xline_idx) - x(xline_idx(1))).^2 + (y(yline_idx) - y(yline_idx(1))).^2); %arclength along line
hold on
xline = x(xline_idx);
yline = y(yline_idx);
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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% make the plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%open the image
t = Tiff('../PIG-S2-NovDec2020.tif', 'r'); %!! not in git repo!!
imageData = read(t);
fig1 = figure(1); clf; hold on
imshow(imageData);
ax = gca;
%add the reference points
load('../plot_tools/image_to_model_points.mat', 'ximage', 'yimage', 'xmod', 'ymod');

%add the 2012 grounding line
[model2image, image2model] = genmaps_image2model(ximage, yimage, xmod, ymod); %maps from model to image and vice versa
c_image_GL = model2image(cgl); %put gl model gl position onto image
scatter(ax, c_image_GL(1,:), c_image_GL(2,:), 4,'k', 'filled')
%contour(xi, yi, bathy', [0,0], 'k', 'linewidth', 2)

%add the fronts
colmap = parula(7);
%colmap(1,:) = [1,0,0]; %make first red
%colmap(end,:) = [0,1,0]; %last row green
for i = 6:-1:1 %reverse order so 2009 on top
    cfront = cell2mat(topo_front_data(i));
    c_image_front = model2image(cfront);
    scatter(c_image_front(1,:), c_image_front(2,:), 3, colmap(i,:), 'filled');
end


%add the calving front measurement line
cl = [xline; yline];
cl_image = model2image(cl);
plot(cl_image(1,:), cl_image(2,:), 'k--', 'linewidth', 1.5);

%add the inner cavity definition
cin1_img = model2image(cin1);
cin2_img = model2image(cin2);
%plot(cin1_img(1,:), cin1_img(2,:), '--','color', plotcolor2);
%plot(cin2_img(1,:), cin2_img(2,:), '--', 'color',plotcolor3);
f(1) = fill(cin1_img(1,:), cin1_img(2,:), 'm', 'FaceAlpha', 0.5, 'LineStyle', 'none');
f(2) = fill(cin2_img(1,:), cin2_img(2,:), 'm' , 'FaceAlpha', 0.5, 'LineStyle', 'none');


%tidy plot
ylim([2194, 14847])
xlim([2000,12000])

%rotate
camroll(-90);

%add the A,B pts and 2009, 2020 front labels
%ax1 = gca;
%txtN = text(ax1, 5200, 6400, 'North', 'FontSize', 12, 'Interpreter', 'latex');
%txtS = text(ax1, 8700, 6000, 'South', 'FontSize', 12, 'Interpreter', 'latex');
%f2009 = text(ax1, 4000,13200, '2009', 'FontSize',12, 'color', colmap(1,:));
%f2020 = text(ax1, 7500,11700, '2020', 'FontSize',12, 'color', colmap(2,:));



