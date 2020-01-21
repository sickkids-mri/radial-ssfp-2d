function I = Draw_ROI(I,ROIy,ROIx)

M=max(I(:));
if ~isempty(ROIy) && ~isempty(ROIx)
if ROIy(1)>1 && ROIy(end)<size(I,1) && ROIx(1)>1 && ROIx(end)<size(I,2) 
I(ROIy(1),ROIx(1:2:end),:)=M;
I(ROIy(end),ROIx(1:2:end),:)=M;
I(ROIy(1:2:end),ROIx(1),:)=M;
I(ROIy(1:2:end),ROIx(end),:)=M;
end
end
end