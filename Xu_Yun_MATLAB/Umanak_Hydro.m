clear

%**********************************
%****** PHI SURFACE - GIMP 90m ****
%**********************************
if 0
[gimp,R]=geotiffread('./gimpdem_90m.tif');
 X=R.XLimWorld(1):R.DeltaX:R.XLimWorld(2);
 Y=R.YLimWorld(2):R.DeltaY:R.YLimWorld(1);

 X=0.5*(X(1:end-1)+X(2:end));
 Y=0.5*(Y(1:end-1)+Y(2:end));

else
 load GIMP_90m.mat
end

%ix=13000:17000; iy=4000:6000;
%ix=15500:16500; iy=4700:6000;
ix=15600:16800; iy=4700:6000;
gimp=gimp(ix,iy);X=X(iy);Y=Y(ix);

%X=repmat(X ,length(Y),1);
%Y=repmat(Y',1,size(X,2));

if 1
%******************************************
%****** LOAD RUNOFF w/ gimp Resolution ****
%******************************************
load RUNOFF_20120712_30000x16620

RUNOFF_ALL=RUNOFF_ALL(ix,iy);

%********************************
%****** PHI BOTTOM - MCdata *****
%********************************
%LOADING MCDATSET
MCdatsetnc='../MCdatset/MCdatset-2013-12-03.nc';

x = double(ncread(MCdatsetnc,'x'));
y = double(ncread(MCdatsetnc,'y'));
MCx=repmat(x,1,length(y));
MCy=repmat(y',size(x,1),1);
  
MCb = double(ncread(MCdatsetnc,'bed'));
%MCsurf = double(ncread(MCdatsetnc,'surface'));
%MCsource = double(ncread(MCdatsetnc,'source'));

MCb(find(MCb<=-9999))=nan; 
%MCb(find(MCsource==4 & MCsurf<500))=nan;

ix=find( x>X(1)   & x<X(end) );
iy=find( y>Y(end) & y<Y(1)   );
x=MCx(ix,iy);y=MCy(ix,iy);b=MCb(ix,iy); 

X=repmat(X ,length(Y),1);
Y=repmat(Y',1,size(X,2));

%%INTERP MULTIBEAN AND MCDATSET TO GIMP GRIDS
tic,BED=interp2(x',y',b',X,Y);toc
%figure,mypcolor(X,Y,BED),colorbar;axis equal; axis tight
end
%**********************************
%*** SMOOTH PHI, FIND QSG route ***
%**********************************

rho_water= 1000; %[kg/m3]
rho_ice  = 917;  %[kg/m3]
gravity  = 9.81; %[m/s2]

PHI=rho_ice*gravity*single(gimp); %[kg/m3]*[m/s2]*[m]=[kg/m/s2]
PHI=PHI+(rho_water-rho_ice)*gravity*BED;

PHI = mysmooth2d(PHI,3,3);

% ------ FLIPUD ------
PHI        = double(flipud(PHI));
X          = flipud(X);
Y          = flipud(Y);
RUNOFF_ALL = double(flipud(RUNOFF_ALL));

% ----- FILL SINKS -----
aaa=find(isnan(PHI));
PHI(aaa)=0;
tic,PHI=imfill(PHI,4,'holes');toc
PHI(aaa)=nan;

tic,[flowacc] = wflowacc(X,Y,PHI,'W0',RUNOFF_ALL,'exponent',5);toc %[kg/m2/s]
flowacc=flowacc*8.1; %[kg/m2/s] * 90m * 90m / 1000kg/m3 = [m3/s]
%tic,[flowacc] = wflowacc(X,Y,PHI);toc %[kg/m2/s]
flowacc(find(isnan(PHI))) = nan;

mypcolor(X,Y,log10(flowacc),0,2) %[m3/s]

