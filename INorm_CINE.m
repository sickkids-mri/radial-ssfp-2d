function IMG = INorm_CINE(IMG)
for j=1:size(IMG,4);
for i=1:size(IMG,3);
    I=IMG(:,:,i,j);
I=I./max(abs(I(:)));
IMG(:,:,i,j)=I;
end
end

end