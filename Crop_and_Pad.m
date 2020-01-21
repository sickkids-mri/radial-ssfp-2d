function I = Crop_and_Pad(I,ROIy,ROIx)
if exist('ROIy','var')
temp=zeros(size(I));
temp(ROIy,ROIx,:,:)=I(ROIy,ROIx,:,:);
I=temp;clear temp
else
I=pad_image(Crop_CINE(I),size(I,1));        
end    
end