function I = Window_CINE(I,mn,mx)
if ~exist('mn','var')
    mn=0;
end
if ~exist('mx','var')
    mx=0.5;
end

if mx~=1
I=abs(I-min(I(:)));
I=I./abs(max(I(:)));
I(I<mn*max(I(:)))=mn*max(I(:));
% I(I<mn*min(I(:)))=mn*min(I(:));
I(I>mx*max(I(:)))=mx*max(I(:));
end
I=I./max(abs(I(:)));
end


