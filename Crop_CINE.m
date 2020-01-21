function [I,y,x] = Crop_CINE(I,sz)
if ~exist('sz','var')
    sz=round(size(I,1)/2);
end
if size(I,1)>=sz && size(I,2)>=sz
x=(size(I,2)-sz)/2;
x=ceil(x+1:x+sz);
y=(size(I,1)-sz)/2;
y=ceil(y+1:y+sz);
I=I(y,x,:,:);
end
end