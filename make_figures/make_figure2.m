%Make figure 2 in the manuscript: schematic of the shelf geometry and hydrographic conditions
%Alex Bradley (aleey@bas.ac.uk), 21/05/2021, MIT Licence


%
% Flags
% 
save_flag = 0; 

%
% Preliminaries
%
addpath("plot_tools");
plot_defaults
fig = figure(1); clf; 
fig.Position(3:4) = [900, 420];
label_size = 12;
ax_fontsize = 12;


%
% (a) Schematic of ice shelf geometry
%

%grid
nx  = 320;
dx  = 400;
yy  = 0:dx:((nx-1)*dx); 
yyf = max(yy) - yy;
[~, idx] = min(abs(yy - 84*1e3)); %ice front index

%generate ice topo profiles
H = [100, 150, 200];
linestyles = ["-","--", "--"];
h_profiles = zeros(3, length(yy));
for i = 1:3
h_profiles(i,:)=(310 + H(i))/2.64*atan(0.17*yy/1000 - 3) + 0.47*(H(i)+400) - 1051.3;
h_profiles(i, idx+1:end) = nan;
end

%make plot
pos = [0.1, 0.12, 0.5, 0.82];
p0 = subplot('Position', pos); box on; hold on
for i = 1:3 %loop over H values
if i == 1
%add inner cavity definition
X = 98:128;
fillX = [X, flip(X)];
fillY = [-1100*ones(size(X)), -400*ones(size(X))];
fill(fillX, fillY, [167, 0, 47]/255, 'linestyle', 'none', 'FaceAlpha', 0.3)

%add shaded region of ice
fillX = [yyf(1:idx), flip(yyf(1:idx))]/1e3;
fillY = [h_profiles(i,1:idx), zeros(1,idx)];
fill(fillX, fillY, [173, 216, 230]/255, 'linewidth', 1.5)

%add calving lines before other H values
lc = [103, 124, 200]/255; 
extent = 40:5:80;
for i = 1:9
	%find the value of H here
	[~,idx] = min(abs(yy - 1e3*extent(i)));
	hmin = h_profiles(length(yy) - idx);
	plot((128 - extent(i))*ones(1,2), [h_profiles(1,idx), 0], 'color', lc)
end

else
plot((max(yy) - yy)/1e3, h_profiles(i,:), 'k', 'linestyle', linestyles(i))
end
end %end loop over h values

%add the ridge
fillX = [yy, flip(yy)]/1e3;
latg = [1.62e6:400:1.748e6-400];
bump = 400*exp(-(latg-1.67e6).^2/(2*12000^2)) - 1095;
fillY = [-1100*ones(1,length(yy)), bump];
fill(fillX, fillY, [203, 150, 80]/255, 'Linewidth', 1.5)
%plot(yyf/1e3, bump, 'k', 'linewidth', 1.5)	


xlabel('y (km)','FontSize',label_size, 'Interpreter', 'latex' )
ylabel('depth(m)', 'FontSize', label_size, 'Interpreter', 'latex')
xlim([min(yy), max(yy)]/1e3)
ylim([-1100, 0])
p0.YTick = [-1100,-1000:200:0];
p0.YLabel.String = 'depth (m)';
grid on

%add north south text
north = text(6,-210, "North", 'FontSize', 16,  'Interpreter', 'latex');
set(north, 'Rotation', 90)

south = text(123,-210, "South", 'FontSize', 16, 'Interpreter', 'latex');
set(south, 'Rotation', 90)

%add the W values
text(25,-480, "W = 200", 'FontSize', 11, 'Interpreter' , 'latex')
text(25,-580, "W = 100", 'FontSize', 11, 'Interpreter', 'latex')
text(25,-530, "W = 150", 'FontSize', 11, 'Interpreter', 'latex')

%add w arrow
%pt1 = [128-50, -650]; pt2 = [128-50, -620];
%dp = pt2 - pt1;
%quiver( pt1(1), pt1(2), dp(1), dp(2), 0, 'k')
%pt1 = [128-50, -650]; pt2 = [128-50, -680];
%dp = pt2 - pt1;
%quiver( pt1(1), pt1(2), dp(1), dp(2), 0, 'k')
%tw = text(128 - 56, -655, "W", 'FontSize', 11, 'Interpreter', 'latex');


%sort out the ticks, which are the wrong way round
p0.XTick = flip(128 - (20:20:120));
p0.XTickLabel = {'120', '100', '80', '60', '40', '20'};


%
% (b) and (c) salinity and temperature profiles
%
pos_t = [0.62, 0.12, 0.17, 0.82]; %subplot positions
pos_s = [0.80, 0.12, 0.17, 0.82];
p1 = subplot('Position',pos_t); box on; hold on; grid on
p2 = subplot('Position',pos_s); box on; hold on; grid on
depth = 0:10:1110;
P = [600, 700. 800];
linestyles = ["-", "--", "-."];

%2009 and 2012 profiles in background
load('./data/TS_2009.mat', "Z", "T2009","S2009");
plot(p1, T2009, Z, 'color', ones(3,1)*0.7)
plot(p2, S2009, Z, 'color', ones(3,1)*0.7)

load('./data/TS_2012.mat', "Z", "T2012","S2012");
plot(p1, T2012, Z, 'color', ones(3,1)*0.5)
plot(p2, S2012, Z, 'color', ones(3,1)*0.5)

for i = 1:3 %for each P value
[t_prof, s_prof] = TS_profile(depth,-1100,P(i) - 600, P(i)-600); %send and third arguments are offset from 600

plot(p1, t_prof,-depth(1:end-1), 'r', 'linestyle', linestyles(i), 'linewidth', 1.5);
axes(p1); txt =  text(0, -P(i) + 190, strcat("P = " ,num2str(P(i))), 'Color', 'r', 'Rotation', -45, 'Interpreter', 'latex');

plot(p2,s_prof, -depth(1:end-1), 'b', 'linestyle', linestyles(i), 'linewidth', 1.5);
axes(p2); txt =  text(34.3, -P(i) + 200, strcat("P = " ,num2str(P(i))), 'Color', 'b', 'Rotation', -45, 'interpreter', 'latex');
end %end loop over P values

%tidy up
p1.XLim = [-1.2, 1.4]; 
p1.YLim = [-1100,0];
p1.YLabel.FontSize = label_size;
p1.XLabel.String = 'Pot. temp. (${}^\circ$C)';
p1.XLabel.Interpreter = 'latex';
p1.XLabel.FontSize = label_size;
p1.YTick = [-1100,-1000:200:0];
p1.FontSize = 10;
p1.YTickLabel = cell(length(p2.YTickLabel),1);
p1.YLabel.Interpreter = 'latex';

p2.YTick = p1.YTick;
p2.YLim = [-1100,0];
p2.XLim = [33.9, 34.8];
p2.XTick = [34, 34.3, 34.6];
p2.YTickLabel = cell(length(p2.YTickLabel),1);
p2.XLabel.String = 'Salinity (psu)';
p2.XLabel.FontSize = label_size;
p2.FontSize = 10;
p2.XLabel.Interpreter = 'latex';
p2.YLabel.Interpreter = 'latex';

%figure labels 
text(p0, -20, -40, "(a)", 'FontSize', 12, 'Interpreter', 'latex')
text(p1,  0.93, -40, "(b)", 'FontSize', 12, 'Interpreter', 'latex')
text(p2, 34.653, -40, "(c)", 'FontSize', 12, 'Interpreter', 'latex')


%
% save
%
if save_flag 
%saveas(gcf, "plots/figure2", 'epsc')
%saveas(gcf, "plots/figure2.png")
end
