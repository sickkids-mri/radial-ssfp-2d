function IMG = ISpace_Interpolation_Weights(RT,CP,MI,nFrames)
if length(CP)<size(RT,3)
   t=linspace(0,range(CP(:)),size(RT,3));
   CP = Calculate_CardiacPhases(t,CP);
end

MI=MI./max(abs(MI(:)));
if isempty(MI)
   MI=ones(length(CP),1);
end
CP=col(CP);
MI=col(MI);
MI=imresize(MI,[length(CP),1]);
IMG=zeros(size(RT,1),size(RT,2),nFrames);
CP=[CP-1;CP;CP+1];
iCP=linspace(0,1,nFrames);
for iFrame=1:nFrames
   x=mygaus(1,iCP(iFrame),mean(diff(iCP)),CP);
   x=(x(1:size(MI,1))+x(size(MI,1)+1:2*size(MI,1))+x(2*size(MI,1)+1:3*size(MI,1))).*MI;
   IMG(:,:,iFrame) = sum(bsxfun(@times,RT,permute(x,[2,3,1])),3)./sum(x(:));
end
end

function f = mygaus(h,c,s,x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function used with leastsq to fit data to
% SF=height * exp(-(x-centre)^2/(2o^2)))
% param(1) = height
% param(2) = center
% param(3) = sigma
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f= h*exp((-(x-c).^2)./(2*s.^2));
end