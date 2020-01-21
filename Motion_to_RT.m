function RT = Motion_to_RT(RT,myTransforms)
[Y,X]=meshgrid(linspace(-0.5,0.5,size(RT,1)),linspace(-0.5,0.5,size(RT,2)));
k=repmat(X+1i*Y,[1,1,size(RT,3)]);
Z=zeros(size(RT),'single');
% par
for iSlice=1:size(RT,4)
iTransforms=imresize(myTransforms(:,:,iSlice),[size(RT,3),2]);
iTransforms=permute(iTransforms,[3,2,1]);
Z(:,:,:,iSlice)=exp(-2*pi*1i*(bsxfun(@times,-iTransforms(:,2,:),real(k))-bsxfun(@times,iTransforms(:,1,:),imag(k))));
end
RT=ifft2c(Z.*fft2c(RT));
end