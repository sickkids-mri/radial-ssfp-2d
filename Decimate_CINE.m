function img = Decimate_CINE(I,n)
for t=1:size(I,3)
    for u=1:size(I,4)
        img(:,:,t,u) = Decimate_Resolution(I(:,:,t,u),n);
    end
end
end