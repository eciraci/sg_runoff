function h=mypcolor(x,y,data,varargin);
%MYPCOLOR - pcolor for Yun
%
%   Usage:
%      h=mypcolor(x,y,data)
%      h=mypcolor(x,y,data,minimum,maximum)
%
%   Created by: Mathieu Morlighem :)

%Get extreme values
if nargin==3,
	data_min=min(data(:));
	data_max=max(data(:));
else
	data_min=varargin{1};
	data_max=varargin{2};
	data(find(data<data_min))=data_min;
	data(find(data>data_max))=data_max;
end

%Keep track of NaN and convert mins
[data_nani data_nanj]=find(isnan(data) | data==-9999);

%colormap
map = jet(64);
map = map(17:end,:);
%map = colormap(cool);
%map = flipud(map);
%map = map(5:end,:);

%Create RGB image from data
image_rgb = ind2rgb(uint16((data - data_min)*(length(map)/(data_max-data_min))),map);

%Raplce NaN by white
if ~isempty(data_nani)
	nancolor=[1. 1. 1.];
	image_rgb(sub2ind(size(image_rgb),repmat(data_nani,1,3),repmat(data_nanj,1,3),repmat(1:3,size(data_nani,1),1))) = repmat(nancolor,size(data_nani,1),1);
end

%plot data
h=imagesc(x(1,:),y(:,1),image_rgb);
axis xy;
caxis([data_min data_max]);
