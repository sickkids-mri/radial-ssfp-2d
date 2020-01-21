function [Sy,Sx] = Enlarge_ROI(ROIy,ROIx,y,x)
% y=y/2;
% x=x/2;
Sy=ROIy(1)-floor(y):ROIy(end)+ceil(y);
Sx=ROIx(1)-floor(x):ROIx(end)+ceil(x);
if min(Sy)<0
Sy=Sy-min(Sy(:))+1;
end
if min(Sx)<0
Sx=Sx-min(Sx(:))+1;
end

end