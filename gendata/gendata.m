% Create data files for the idealized PIG snapping case (interp_APIGi)
% Code written by Alex Bradley (aleey) Feb 2021  based on interp_PIG_ideal.m (Jan De Rydt c. 2013)

%close all
%add interp tools dir
addpath('interp_tools');

% Flags
RBCSvelFlag = 0;
OBCSFlag = 0;

% Input parameters
a = 400; %height of ridge
W = 100; %thickness of gap
dzT = 0; %shift of thermocline from 200m depth (positive = lowering)
dzS = 0; %shift of the halocline from 200m depth (positive = lowering)
%for extent = [84,80,75,70,65,60,55,50,45,40,35,30]*1e3; %define extent of ice shelf (from 0) [default = 84*1e3 = 84km]
for extent = [84*1e3]
%Plot flags
plot_slice_from_gl = 0; %plot ice topo slice as a function of distance from gl
plot_slice_in_stereographic = 0; %plot ice topo slice in stereographic coords 
plot_ST_profiles = 0; %plot the temperatrue and salinity profile
plot_pload_plan = 0;
plot_pload_slice = 0;
numplot = 1;

%%
%% Grid
%%
nx=120; % number of grid cells along longitudinal direction
ny=320; % number of grid cells along latitudinal direction
nz=110; % number of vertical grid cells

dlon = 400; dx=dlon; % horizontal grid resolution (meters)
dlat = 400; dy=dlat;
deltaZ = 10; % vertical grid resolution (meters)

long = [0:dlon:48e3-dlon];
lonc = long+dlon/2;
latg = [1.62e6:dlat:1.748e6-dlat];
latc = latg+dlat/2;
dz = deltaZ*ones(1,nz);
zgp1 = [0,cumsum(dz)];
zc = .5*(zgp1(1:end-1)+zgp1(2:end));
zg = zgp1(1:end-1);

hFacMin=0.05;

%%
%% Model parameters
%%

% Equation of state
%eos = 'linear';
%eos = 'jmd95z';
eos = 'mdjwf';

% Precision for writing binary files
acc = 'real*8';

% Model geometry
H = -nz*deltaZ; % nominal depth of the model

% Reference temperature and salinity
Tref = -1;
Sref = 34;

% Physical constants
gravity=9.81;
rhoConst = 1030;
si2dbar = 1e-4;

%%
%% Create ice topography
%%

LAT=[0:dlat:(ny-1)*dlat];
[~,idx_extent] = min(abs(LAT - extent)); %find index where ice passes through ice front 
% ATAN profile
%iceprofile=160*atan(0.17*LAT/1000-3)-720;
iceprofile=(-90 + W + a)/2.64*atan(0.17*LAT/1000 - 3) + 0.47*(W+a) - 1051.3;
% thicken i2ce2
%iceprofile = iceprofile;
% no negative values
%ind=icetopo<0
%icetopo(ind)=0;
%icetopo = ones(nx,1)*iceprofile-100;
icetopo = ones(nx,1)*iceprofile;

% Van de Veen profile
%H0=1150;
%U0=6743/(365*24*60*60); %375 (400m gap), 1203 (300m gap), 3041 (200m gap), 3745 (175m gap), 4581 (150m gap), 5571 (125m gap), 6743 (100m gap) 
%Ri=917;
%Rw=1030;
%A=3.5e-25;
%n=3;
%Cs=A*(1/4*Ri*gravity*(1-Ri/rhoConst))^n;
%U=(U0^(n+1)+(n+1)*Cs*(U0*H0)^n*LAT).^(1/(n+1));
%IceThick=U0*H0./U;
%h=(1-Ri/rhoConst)*IceThick;
%iceprofile=h-IceThick;
%icetopo=ones(nx,1)*iceprofile;

%figure; plot(latg,iceprofile);

% channels
nchannels = 0;
hmax = 150;
hmin = 0;
width = 20000;
   
channels=zeros(nx,ny);
if nchannels~=0;
    spacing=nx/(2*nchannels);
    cl=round(spacing:2*spacing:nx);
    h=zeros(ny,1);
    w=-floor(width/(2*dlon)):floor(width/(2*dlon));
    w=w(end-floor(width/dlon)+1:end);
    dIy1=5000/dlat;
    dIy2=50000/dlat;
    for ii=1:ny
        if LAT(ii)<5000
            h(ii)=(hmax-hmin)/(dIy1)*(ii-5000/dlat)+hmax;
        elseif LAT(ii)<40000 && LAT(ii)>=5000
            h(ii)=hmax;
        elseif LAT(ii)<90000 && LAT(ii)>=40000
            h(ii)=(hmax-hmin)/(-dIy2)*(ii-90000/dlat)+hmin;
        elseif LAT(ii)>=90000
            h(ii)=hmin;
        end
    end
    %figure; plot(latg,h);
    xsection=zeros(nx,1);
    for ii=1:ny
        for jj=1:nchannels
        chpr=h(ii)*sin(pi*([0:length(w)+1])/(length(w)+1));
        xsection(w+cl(jj))=chpr(2:end-1);
        end
        channels(:,ii)=xsection;
    end
end

% combine flat profile and channels
icetopo=icetopo+channels;
icetopo(:,idx_extent:end)=0; % create ice front

%figure; imagesc(icetopo);

[ix,iy]=find(icetopo~=0);
is=find(icetopo~=0);

%%
%% Create bathymetry (box with open boundary in the north)
%%

bathy = ones(nx,ny)*H; % ocean floor

% gaussian ridge
%a=400; %comment out if specified above
b=1.67e6;
c=12000;

bump=a*exp(-(latg-b).^2/(2*c^2));
%bump=0;

for ii=1:nx
    bathy(ii,:)=bathy(ii,:)+bump;
end

bathy(1,:) = 0; % western sidewall
bathy(:,1) = 0; % southern sidewall
bathy(end,:)=0; % eastern sidewall


%figure; imagesc(bathy);

%%
%% Check consistency of icetopo and bathymetry
%%

for ii=1:nx
   for jj=1:ny
       if icetopo(ii,jj) < bathy(ii,jj)
           icetopo(ii,jj) = 0 ;
           bathy(ii,jj)= 0 ;
       end
   end
end

%icetopo(:,1:123)=0;
%bathy(:,1:123)=0;

%%
%% check depth of water column thickness (min = 100m)
%%

for ii=1:nx
    for jj=1:ny
        if icetopo(ii,jj) - bathy(ii,jj) < 80 && bathy(ii,jj) ~= 0
            icetopo(ii,jj) = bathy(ii,jj) + 80;
        end
    end
end

%icetopo(:,1:126)=0;
%bathy(:,1:126)=0;

% plot

if plot_slice_from_gl; %plot ice topo slice as a function of distance from gl
figure(numplot);clf;  plot(LAT, icetopo(2,:), LAT, bathy(2,:)); grid on
numplot = numplot + 1;
end
if plot_slice_in_stereographic
figure(numplot); clf; hold on; plot(latg,icetopo(2,:),latg,bathy(2,:)); grid on;
numplot = numplot + 1;
end

% save
saveFolder = 'bathy_files';
if ~exist(saveFolder, 'dir');
mkdir(saveFolder);
end
fname = ['bathymetry_H' num2str(a) '.shice'];
fid=fopen(strcat(saveFolder, '/', fname),'w','b'); fwrite(fid,bathy,acc); fclose(fid);


saveFolder = 'topo_files';
if ~exist(saveFolder, 'dir');
mkdir(saveFolder);
end
fname = ['shelfice_topo_H' num2str(a) '_W' num2str(W) '_extent' num2str(extent/1e3) 'km.bin'];
fid=fopen(strcat(saveFolder, '/', fname),'w','b'); fwrite(fid,icetopo,acc); fclose(fid);

%%
%% Set Theta and Salinity boundary conditions (here: northern RBCs)
%%

[t_profile_or,s_profile_or]=TSprofile_APIGi(zgp1,H,nz,dzT,dzS);
%t_profile=smoothn(t_profile_or,10);
%s_profile=smoothn(s_profile_or,10);
t_profile=t_profile_or;
s_profile=s_profile_or;

if plot_ST_profiles
figure(numplot);clf; subplot(1,2,1); plot(t_profile,-zc,t_profile_or,-zc); grid on; hold on;
subplot(1,2,2);plot(s_profile,-zc,s_profile_or,-zc); grid on; hold on;
end
if OBCSFlag ~= 1 
    % Set masks
    mask_T=zeros(nx,ny,nz);
    mask_S=zeros(nx,ny,nz);
    if RBCSvelFlag == 1
        mask_U=zeros(nx,ny,nz);
        mask_V=zeros(nx,ny,nz);
    end

    % relax at northern boundary over 5 grid cells
    for ii=0:4
        mask_T(:,ny-ii,:)=1-ii*0.20;
        mask_S(:,ny-ii,:)=1-ii*0.20;
    end
    if RBCSvelFlag == 1
        for ii=0:4
            mask_U(:,ny-ii,:)=1-ii*0.20;
            mask_V(:,ny-ii,:)=1-ii*0.20;
        end
    end

    % Set boundary condtions
    relaxTFile=zeros(nx,ny,nz);
    relaxSFile=zeros(nx,ny,nz);
    for jj=0:4
        for ii=1:nx
            relaxTFile(ii,ny-jj,:)=t_profile;
            relaxSFile(ii,ny-jj,:)=s_profile;
        end
    end
    if RBCSvelFlag == 1
        relaxUFile=zeros(nx,ny,nz);
        relaxVFile=zeros(nx,ny,nz);
    end

    saveFolder = 'RBC_files';
    if ~exist(saveFolder, 'dir');
	mkdir(saveFolder);
    end
    fnamet = ['RBCt_dzT_' num2str(dzT) '_dzS_' num2str(dzS) '_nz' num2str(nz) '.bin'];
    fnames = ['RBCs_dzT_' num2str(dzT) '_dzS_' num2str(dzS) '_nz' num2str(nz) '.bin'];
    fid=fopen(strcat(saveFolder, '/', fnamet),'w','b'); fwrite(fid,relaxTFile,acc); fclose(fid);
    fid=fopen(strcat(saveFolder, '/', fnames),'w','b'); fwrite(fid,relaxSFile,acc); fclose(fid);

    saveFolder = 'RBC_mask_files';
    if ~exist(saveFolder, 'dir');
	mkdir(saveFolder);
    end
    fid=fopen(strcat(saveFolder,strcat('/RBCt_mask_nz', num2str(nz), '.bin')),'w','b'); fwrite(fid,mask_T,acc); fclose(fid);
    fid=fopen(strcat(saveFolder,strcat('/RBCs_mask_nz', num2str(nz), '.bin')),'w','b'); fwrite(fid,mask_S,acc); fclose(fid);
   
    if RBCSvelFlag == 1
        fid=fopen('RBCu_mask.bin','w','b'); fwrite(fid,mask_U,acc); fclose(fid);
        fid=fopen('RBCv_mask.bin','w','b'); fwrite(fid,mask_V,acc); fclose(fid);
        fid=fopen('RBCu.bin','w','b'); fwrite(fid,relaxUFile,acc); fclose(fid);
        fid=fopen('RBCv.bin','w','b'); fwrite(fid,relaxVFile,acc); fclose(fid);
    end

else
    Nb_t=ones(nx,1)*t_profile';
    Nb_s=ones(nx,1)*s_profile';

    fid=fopen('OBNt.bin','w','ieee-be'); fwrite(fid,Nb_t,'real*8'); fclose(fid);
    fid=fopen('OBNs.bin','w','ieee-be'); fwrite(fid,Nb_s,'real*8'); fclose(fid);
end 
%%
%% Initial conditions 
%%

% create hydrographic fields
levt=zeros(nx,ny,nz);
levs=zeros(nx,ny,nz);

% Northern boundary conditions applied for entire domain
for ii=1:nz
    levt(:,:,ii)=t_profile(ii);
    levs(:,:,ii)=s_profile(ii);
end

% save
saveFolder = 'hydrog_files';
if ~exist(saveFolder, 'dir');
mkdir(saveFolder);
end
fnamet = ['lev_t_dzT_' num2str(dzT) '_dzS_' num2str(dzS) '_nz' num2str(nz) '.shice'];
fnames = ['lev_s_dzT_' num2str(dzT) '_dzS_' num2str(dzS) '_nz' num2str(nz) '.shice'];
fid=fopen(strcat(saveFolder,'/', fnamet),'w','b'); fwrite(fid,levt,acc); fclose(fid);
fid=fopen(strcat(saveFolder,'/', fnames),'w','b'); fwrite(fid,levs,acc); fclose(fid);

%%
%% Compute potential field underneath ice shelf
%%

dzm = abs([zg(1)-zc(1) .5*diff(zc)]);
dzp = abs([.5*diff(zc) zc(end)-zg(end)]);
%hFacMin=0.05; % in principle needs to be the same as parameter in data file for MITgcm

for ks=1:length(ix)
t0 = squeeze(levt(ix(ks),iy(ks),:)); 
s0 = squeeze(levs(ix(ks),iy(ks),:));
p = abs(zc(:))*gravity*rhoConst*si2dbar;
dp = p;
kp = 0;

while rms(dp) > 1e-13
  phiHydF(1) = 0;
  p0 = p;
  kp = kp+1;

    switch eos
     case 'mdjwf'
      drho = densmdjwf(s0,t0,p(:,end))-rhoConst;
     otherwise
      error(sprintf('unknown eostype: %s',eos))
    end
    
    for k = 1:nz 
        phiHydC(k)   = phiHydF(k) + dzm(k)*gravity*drho(k)/rhoConst;
        phiHydF(k+1) = phiHydC(k) + dzp(k)*gravity*drho(k)/rhoConst;
    end
    
  switch eos
    case 'mdjwf'
        p = (gravity*rhoConst*abs(zc(:)) + phiHydC(:)*rhoConst)*si2dbar;%/gravity/rhoConst;
      otherwise
          error(sprintf('unknown eostype: %s',eos))
    end
  dp = p-p0;
end
  
% find appropriate level
   zloc = icetopo(ix(ks),iy(ks));
   kl = max(find(abs(zg)<=abs(zloc)));
   % kl = max(find(abs(zg-hFacMin*zg)<=abs(zloc)));
  % hfloc= squeeze(hf(ix(ks),iy(ks),:));
  % kl = min(find(hfloc>0));
   if isempty(kl);
     kl = 0;
     ph(ks) = 0;
   end
   if kl > 0
       klp1 = min(kl+1,nz);
       %% compute hFacC
        dHfrac=(-zg(kl)-zloc)/(zg(klp1)-zg(kl));
        if dHfrac < hFacMin
            drloc=0;
        elseif dHfrac > 1-hFacMin
            drloc=1;
        else
            drloc=dHfrac;
        end
       %%
       dph=phiHydF(klp1)-phiHydF(kl);
       ph(ks) = phiHydF(kl)+drloc*dph;
   end
  % disp(sprintf('kl0 = %u, kl = %u',kl0,kl));
end

pload = zeros(nx,ny);
for ks=1:length(ix)
      pload(ix(ks),iy(ks)) = ph(ks)*rhoConst;
end

if plot_pload_plan
figure(numplot); clf; imagesc(lonc,latc,pload'); set(gca,'YDir','normal'); axis equal;
numplot = numplot + 1;
end
if plot_pload_slice
figure(numplot); clf; plot(latc,pload(2,:),latc,icetopo(2,:));
end

saveFolder = 'pload_files';
if ~exist(saveFolder, 'dir');
mkdir(saveFolder);
end

fname = ['pload_H' num2str(a) '_W' num2str(W) '_extent' num2str(extent/1e3) 'km.'];
fname = strcat(saveFolder, '/', fname, eos);
fid=fopen(fname,'w','b'); fwrite(fid,pload,acc); fclose(fid);
end
