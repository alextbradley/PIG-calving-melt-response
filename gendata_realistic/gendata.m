%
% generate input files for rPIG case
% Based on the interp_PISOMIP code by Paul Holland and interp_PIG code of Jan de Rydt. Designed for c65t.
% Alex Bradley 05/12/20

close all
clear all
addpath('MITgcm_matlab_utils'); %specify path for matlab utilities
datapath = 'data/'; %path to locate data
savepath = 'gendata/'; %path to save output to
%% Flags

DiffFlag=0;
OBCSFlag=0;
WindFlag=0;
SnappedFlag = 0;
Benchmarking = 0; %flag to specify comparing output files to the benchmark JDR 2009 files
    wind='real';
    %wind='synth';
LatLon=0;
STFlag=3; % options:    '0: artificial S and T profiles' -> SPECIFY DEPTH OF PYCNOCLINE AND BOTTOM WATER PROPERTIES !!!,
%                       '1: realistic S and T profiles from 2009',
%                       '2: realistic S and T profiles from 2009, reconstructed from model output'
%                       '3: realistic S and T profiles from 2012, reconstructed from model output'

Geometry ='JDR'; %options: JDR: geometry provided by Jan de Rydt
%
%% grid settings
%

nx=360;
ny=320;
nz=120;

delX = 400;
delY = 400;
delZ = 120;

% Make grids
%[x,y,deltaX,deltaY,deltaZ]=PrepareGrid(nx,ny,LatLon);
deltaZ = delZ; %hack before we've figured out PrepareGrid
dz = deltaZ*ones(1,nz);
zgp1 = [0,cumsum(dz)];
zc = .5*(zgp1(1:end-1)+zgp1(2:end));
zg = zgp1(1:end-1);
dz = diff(zgp1);
%% Model Parameters
%
%eos = 'linear';    % equation of state 
%eos = 'jmd95z';
eos = 'mdjwf';

acc = 'real*8';    % precision to write out

% Model geometry
%H = -nz*deltaZ; % nominal depth of the model

% Reference temperature and salinity
Tref = 1;
Sref = 34.8;

% Physical constants
gravity=9.81;
rhoConst = 1030;
si2dbar = 1e-4;

%
%% Ice shelf and bathymetry data
if strcmp(Geometry, 'JDR')
    if (nx ~= 360 ) || (ny ~= 320) || (nz ~= 120)
        error('grid size incompatible with JDR dataset')
    else
        if SnappedFlag
            error('need to add snapped data');
        else
            load(strcat(datapath,'JDR_icetopo'), 'icetopo');
            load(strcat(datapath,'JDR_bathy'), 'bathy');
        end
    end
end
[ix,iy]=find(icetopo~=0);
is=find(icetopo~=0);

fid=fopen(strcat(savepath, 'shelfice_topo.bin'),'w','b'); fwrite(fid,icetopo,acc); fclose(fid);
fid=fopen(strcat(savepath,'bathymetry.shice'),'w','b'); fwrite(fid,bathy,acc); fclose(fid);


%% Theta and Salinity boundary conditions
%note JDR uses code "TSProfile" with input STFlag generate t and s
%profiles, rather than if statements
if STFlag == 1
    load(strcat(datapath, '2009_refTS_JDR.mat'));
    t_profile = squeeze(ref_temp);
    s_profile = squeeze(ref_salt);
elseif STFlag == 2    
    load(strcat(datapath, '2009forcing_refTS_constructed_from_JDR_output.mat'));
    t_profile = squeeze(ref_temp);
    s_profile = squeeze(ref_salt);
elseif STFlag == 3
    load(strcat(datapath, '2012forcing_refTS_constructed_from_JDR_output.mat'));
    t_profile = squeeze(ref_temp);
    s_profile = squeeze(ref_salt);
    
end

% Set masks
mask_T=zeros(nx,ny,nz);
mask_S=zeros(nx,ny,nz);
mask_U=zeros(nx,ny,nz);
mask_V=zeros(nx,ny,nz);

% relax at western/northern boundary over 5 grid cells
for ii=0:4
    mask_T(1+ii,:,:)=1-ii*0.125;
    mask_S(1+ii,:,:)=1-ii*0.125;
    mask_T(:,ny-ii,:)=1-ii*0.125;
    mask_S(:,ny-ii,:)=1-ii*0.125;
end


% Set boundary condtions
relaxTFile=zeros(nx,ny,nz);
relaxSFile=zeros(nx,ny,nz);
relaxUFile=zeros(nx,ny,nz);
relaxVFile=zeros(nx,ny,nz);

for ii=1:5
    for jj=1:ny
        relaxTFile(ii,jj,:)=t_profile;
        relaxSFile(ii,jj,:)=s_profile;
    end
end
for jj=ny-4:ny
    for ii=1:nx
        relaxTFile(ii,jj,:)=t_profile;
        relaxSFile(ii,jj,:)=s_profile;
    end
end


fid=fopen(strcat(savepath,'RBCt_mask.bin'),'w','b'); fwrite(fid,mask_T,acc); fclose(fid);
fid=fopen(strcat(savepath,'RBCs_mask.bin'),'w','b'); fwrite(fid,mask_S,acc); fclose(fid);
fnameRBCt = strcat(savepath,'RBCt_STFlag', num2str(STFlag) ,'.bin');
fnameRBCs = strcat(savepath,'RBCs_STFlag', num2str(STFlag) ,'.bin');
fid=fopen(fnameRBCt,'w','b'); fwrite(fid,relaxTFile,acc); fclose(fid);
fid=fopen(fnameRBCs,'w','b'); fwrite(fid,relaxSFile,acc); fclose(fid);
%
%% Initial conditions
%%

% create hydrographic fields
levt=zeros(nx,ny,nz);
levs=zeros(nx,ny,nz);

% Western boundary conditions applied for entire domain
for ii=1:nz
    levt(:,:,ii)=t_profile(ii);
    levs(:,:,ii)=s_profile(ii);
end

% save
fnamelev_t = strcat(savepath,'lev_t_STFlag', num2str(STFlag) ,'.shice');
fnamelev_s = strcat(savepath,'lev_s_STFlag', num2str(STFlag) ,'.shice');
fid=fopen(fnamelev_t,'w','b'); fwrite(fid,levt,'real*8'); fclose(fid);
fid=fopen(fnamelev_s,'w','b'); fwrite(fid,levs,'real*8'); fclose(fid);

%% Compute potential field underneath ice shelf
%%

dzm = abs([zg(1)-zc(1) .5*diff(zc)]);
dzp = abs([.5*diff(zc) zc(end)-zg(end)]);

for ks=1:length(ix)

t0 = squeeze(levt(ix(ks),iy(ks),:));
s0 = squeeze(levs(ix(ks),iy(ks),:));
p = abs(zc(:))*gravity*rhoConst*si2dbar;
dp = p;
kp = 0;
tol1=1;
tol2=2;

while tol1 > 0
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
        p = (gravity*rhoConst*abs(zc(:)) + phiHydC(:)*rhoConst)*si2dbar;
      otherwise
          error(sprintf('unknown eostype: %s',eos))
    end
  dp = p-p0;
  tol2=tol1;
  %tol1=rms(dp);
  tol1 = sqrt(mean(dp.^2)); %need signal processing toolbox for above calc
  if tol1==tol2; break; end
end

% find appropriate level
   zloc = icetopo(ix(ks),iy(ks));
   kl = max(find(abs(zg)<=abs(zloc)));
   if isempty(kl);
     kl = 0;
     ph(ks) = 0;
   else
     ph(ks) = phiHydF(kl);
   end
end



pload = zeros(nx,ny);
for ks=1:length(ix)
      pload(ix(ks),iy(ks)) = ph(ks)*rhoConst;
end

fid=fopen(strcat(savepath,['pload.' eos]),'w','b'); fwrite(fid,pload,acc);fclose(fid);

%% Benchmarking
%compare output files with those provided by JDR (useful test of the code)
if Benchmarking
    disp('comparing output files with the benchmark')
    filenames = [strcat(savepath,"bathymetry.shice"),
                strcat(savepath, "lev_s.shice"),
                strcat(savepath, "lev_t.shice"),
                strcat(savepath, "pload.mdjwf"),
                fnameRBCs, 
                fnameRBCt, 
                strcat(savepath,  "RBCs_mask.bin"),
                strcat(savepath,  "RBCt_mask.bin"),
                strcat(savepath,   "shelfice_topo.bin")];
            
    filenames_benchmark = ["benchmark/bathymetry.shice",
                "benchmark/lev_s.shice",
                "benchmark/lev_t.shice",
                "benchmark/pload.mdjwf",
                "benchmark/RBCs.bin",
                "benchmark/RBCt.bin",
                "benchmark/RBCs_mask.bin",
                "benchmark/RBCt_mask.bin",
                "benchmark/shelfice_topo.bin"];
            
    for i = 1:length(filenames)
        %load our file
        fid = fopen(filenames(i));
        var = fread(fid,'real*8', 'b');
        
        %load benchmark file
        fid_benchmark = fopen(filenames_benchmark(i));
        var_benchmark = fread(fid_benchmark,'real*8', 'b');
        
        if any((var_benchmark - var) > 1e-5)
            disp(strcat("failed benchmark test: ", filenames(i)));
            warning('failed test')
        else
            disp(strcat("passed benchmark test: ", filenames(i)));
        end
    end
end
        
   


