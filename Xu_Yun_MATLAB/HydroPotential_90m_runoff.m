clear

steps=[3];

%**********************************
%****** PHI SURFACE - GIMP 30m ****
%**********************************
if any(steps==1)
[gimp,R]=geotiffread('./gimpdem_90m.tif');
 X1=R.XLimWorld(1):R.DeltaX:R.XLimWorld(2);
 Y1=R.YLimWorld(2):R.DeltaY:R.YLimWorld(1);

 X=0.5*(X1(1:end-1)+X1(2:end));
 Y=0.5*(Y1(1:end-1)+Y1(2:end));

%X=repmat(X ,length(Y),1);
%Y=repmat(Y',1,size(X,2));

%PLOT TO TEST
%mypcolor(X,Y',gimp,0,2000);colorbar;axis equal;axis tight

%???
%gimp(find(gimp<35))=nan;

 save GIMP_90m X1 Y1 X Y gimp 
end

%********************************************
%****** PHI BOTTOM - MCdata + MULTIBEAM *****
%********************************************
if any(steps==2)
load GIMP_90m
X=repmat(X ,length(Y),1); %./1000;
Y=repmat(Y',1,size(X,2)); %./1000;

%LOADING MCDATSET
MCdatsetnc='../MCdatset/MCdatset-2013-12-03.nc';

x = double(ncread(MCdatsetnc,'x')); %./1000;
y = double(ncread(MCdatsetnc,'y')); %./1000;
MCx=repmat(x,1,length(y));
MCy=repmat(y',size(x,1),1);

MCb = double(ncread(MCdatsetnc,'bed'));
%MCh = double(ncread(MCdatsetnc,'thickness'));
MCsurf = double(ncread(MCdatsetnc,'surface'));
MCsource = double(ncread(MCdatsetnc,'source'));

MCb(find(MCb<=-9999))=nan;
%MCb(find(MCsource>2 & MCb<0))=nan; %Version_1  
MCb(find(MCsource==4 & MCsurf<500))=nan;

%%INTERP MULTIBEAN AND MCDATSET TO GIMP GRIDS
tic,BED=interp2(MCx',MCy',MCb',X,Y);toc
%figure,mypcolor(X,Y,BED),colorbar;axis equal; axis tight

BED=single(BED);
save MCBED_90m BED

if 0
figure(2),clf
subplot(2,1,1)
mypcolor(X,Y,double(gimp),0,2000);colorbar
axis equal; axis tight
hold on,scatter(-213,-2133,'mo','filled')

subplot(2,1,2)
mypcolor(X,Y,BED);colorbar;
axis equal; axis tight
hold on,scatter(-213,-2133,'mo','filled')
end

end %if any(steps==2)

%**********************************
%*** SMOOTH PHI, FIND QSG route ***
%**********************************

if any(steps==3)

if 0 
%****** PHI SURFACE ******
load GIMP_90m
X=repmat(X ,length(Y),1);
Y=repmat(Y',1,size(X,2));

rho_water=1000; %[kg/m3]
rho_ice  =917;  %[kg/m3]
gravity=9.81;   %[m/s2]

PHI=rho_ice*gravity*single(gimp); %[kg/m3]*[m/s2]*[m]=[kg/m/s2]

%****** PHI BOTTOM ******
load MCBED_90m
PHI=PHI+(rho_water-rho_ice)*gravity*BED;

clear gimp BED

%****** LOAD RUNOFF ******
if 0
 fname='./XGRN11_RU_daily_1980_2012/runoff.KNMI-2012.XGRN11.CLRUN.nc';
  runoffdata=ncread(fname,'runoff');
  lat = double(ncread(fname,'lat'));
  lon = double(ncread(fname,'lon'));
  [x,y]=ll2xy(lat,lon,1);
  runoff=runoffdata(:,:,1,194);

  save RUNOFF_20120712 x y runoff
%else
% load RUNOFF_20120712
%end

  RUNOFF_ALL  =single(zeros(size(gimp)));

% *** PARTITION INTO TILES ***
  [ixtile,iytile]=meshgrid(1:12,1:12);
  order=(ixtile-6.5).^2+(iytile-6.5).^2;
  tile=[order(:),ixtile(:),iytile(:)];
  tile=sortrows(tile);
  nxtile=2500; %30000/12;
  nytile=1385; %16620/12;

% *** INTERPOLATE RUNOFF ***
  for itile=1:length(tile)
  disp(itile)
    ixtile=tile(itile,2);iytile=tile(itile,3);
    ix0=(ixtile-1)*nxtile+1;	if ix0<1, ix0=1;end
    ix1=(ixtile  )*nxtile;	if ix1>30000, ix1=30000;end
    iy0=(iytile-1)*nytile+1;	if iy0<1, iy0=1;end
    iy1=(iytile  )*nytile;	if iy1>16620, iy0=16620;end

    ix=ix0:ix1;  iy=iy0:iy1;
    XX=X(ix,iy); YY=Y(ix,iy);
    tic,runoff_tile=griddata(x,y,runoff,XX,YY);toc
    aaa=find(runoff_tile~=0);
     if ~isempty(aaa) 
       RUNOFF_ALL(ix,iy)=runoff_tile;
     end
  end

  save RUNOFF_20120712_30000x16620 RUNOFF_ALL

else

  load RUNOFF_20120712_30000x16620

end

%****************************
%      FLOW ACCUMULATION
%****************************

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

X=X(1,:);Y=Y(:,1);
PHI=single(PHI);
RUNOFF_ALL=single(RUNOFF_ALL);
save READY2FLOWACC PHI X Y RUNOFF_ALL

else 

load READY2FLOWACC

end

X=repmat(X,length(Y),1);
Y=repmat(Y,1,size(X,2));

PHI=double(PHI);
RUNOFF_ALL=double(RUNOFF_ALL);

% --- FLOWACC in each tile ---
FLOWACC_ALL = zeros(size(PHI));

% PARTITION INTO TILES
[ixtile,iytile]=meshgrid(1:12,1:12);
order=(ixtile-6.5).^2+(iytile-6.5).^2;
tile=[order(:),ixtile(:),iytile(:)];
tile=sortrows(tile);
nxtile=2500; %30000/12;
nytile=1385; %16620/12;

for itile= 21 %1:length(tile)
disp(itile)
    ixtile=tile(itile,2);iytile=tile(itile,3);
    ix0=(ixtile-1)*nxtile;	if ix0<1, ix0=1;end
    ix1=(ixtile  )*nxtile+1;	if ix1>30000, ix1=30000;end
    iy0=(iytile-1)*nytile;	if iy0<1, iy0=1;end
    iy1=(iytile  )*nytile+1;	if iy1>16620, iy1=16620;end

    ix=ix0:ix1;iy=iy0:iy1;
    phi        = PHI(ix,iy);
    x_grid     = X(ix,iy);
    y_grid     = Y(ix,iy);
    runoff_grid= RUNOFF_ALL(ix,iy)+FLOWACC_ALL(ix,iy);
    
mypcolor(x_grid,y_grid,phi);
break
 
   aaa=find(runoff_grid~=0);
   if ~isempty(aaa)
     %--- FLOW ACCUMULATION ---
     tic,[flowacc] = wflowacc(x_grid,y_grid,phi,'W0',runoff_grid);toc     
     FLOWACC_ALL(ix,iy)=flowacc;
   end
end

flowacc=FLOWACC_ALL(12000:16000,4000:7000);
phi=PHI(12000:16000,4000:7000);
x_grid=X(12000:16000,4000:7000);
y_grid=Y(12000:16000,4000:7000);
aaa=find(isnan(phi));
flowacc(aaa)=nan;
break

%***** FLOW ACCUMULATION ******
if 1
%http://www.mathworks.com/matlabcentral/fileexchange/14504-flow-accumulation-upslope-area/content/wflowacc.m
tic,[flowacc] = wflowacc(x_grid,y_grid,phi,'W0',runoff_grid);toc
%tic,[flowacc] = wflowacc(x_grid,y_grid,phi);toc
flowacc(find(isnan(phi)))=nan;

%--- PLOT ---
figure,clf
%subplot(2,1,1)
%mypcolor(x_grid,y_grid,phi/0.917/9.81);colorbar
%title('Hydrostatic Potential (converted to meters of ice)','fontsize',14)
%axis equal; axis tight
%%axis([-2.7e5 -2.25e5 -1.985e6 -1.94e6])

%subplot(2,1,2)
mypcolor(x_grid,y_grid,flowacc); %??? WHAT IS THE UNIT???
hc=colorbar;
%set(h2, 'YScale', 'log')
title('Flow map of subglacial discharge','fontsize',14)
axis equal; axis tight
set(gcf,'color','w')
%axis([-2.7e5 -2.25e5 -1.985e6 -1.94e6])
%axis([-2.2e5 -1.8e5 -2.22e6 -2.16e6]) %EQIP
end

end
