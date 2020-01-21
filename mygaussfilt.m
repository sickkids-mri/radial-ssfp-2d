function I = mygaussfilt(I,w)
for i=1:size(I,3)
    for j=1:size(I,4)
        for k=1:size(I,5)
            I(:,:,i,j,k)=imgaussfilt(abs(I(:,:,i,j,k)),w);
        end
    end
end
end