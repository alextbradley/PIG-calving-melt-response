function [t_profile,s_profile] = TS_profile(zgp1,H, dzT, dzS)
% Create the piecewise linear TS profiles used in APIGi
% Inputs:
% zgp1 	: 	(nz x 1) array - Mesh
% H	:	Scalar - Depth of sea bed
% nz 	:	Integer - Number of grid points in z
% dzT 	:	Shift of thermocline away from 600m default (positive = lowering)
% dzS 	: 	Shift of halocline away from 600m defualt (positive = lowering)

nz = length(zgp1);
dT=0;%[-0.1 0];%[-0.4 -0.2 0 0.2 0.4];
dS=0;%[-0.0125 0];%[-0.05 -0.025 0 0.025 0.05];

% T profile
ii=2;
while ii<=length(zgp1)
   if zgp1(ii)<200+dzT
       t_profile(ii-1)=-1;
   elseif zgp1(ii)<=600+dzT
       t_profile(ii-1)=(2.2+dT)/400*(zgp1(ii)-200-dzT)-1;
   else
       t_profile(ii-1)=1.2+dT;
   end
   ii=ii+1;
end
ind=t_profile==0; t_profile(ind)=eps;

% S profile 
ii=2;
while ii<=length(zgp1)
   if zgp1(ii)<=200+dzS
       s_profile(ii-1)=34;
   elseif zgp1(ii)<600+dzS
       s_profile(ii-1)=(0.7+dS)/400*(zgp1(ii)-200-dzS)+34;
   else
       s_profile(ii-1)=(34.7+dS);
   end
   ii=ii+1;
end
end
