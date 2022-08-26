clc
clear

%**********************************
%****** PHI SURFACE - GIMP 30m ****
%**********************************
[gimp,R]=geotiffread('/Volumes/Extreme Pro/GIMP_DEM/gimpdem_90m_v01.1.tif');
X1=R.XLimWorld(1):R.DeltaX:R.XLimWorld(2);
Y1=R.YLimWorld(2):R.DeltaY:R.YLimWorld(1);

X=0.5*(X1(1:end-1)+X1(2:end));
Y=0.5*(Y1(1:end-1)+Y1(2:end));

%X=repmat(X ,length(Y),1);
%Y=repmat(Y',1,size(X,2));

%PLOT TO TEST
figure(1),clf
mypcolor(X,Y',gimp,0,2000);colorbar;axis equal;axis tight

%???
%gimp(find(gimp<35))=nan;

save GIMP_90m X1 Y1 X Y gimp 