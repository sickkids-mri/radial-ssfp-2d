function im=pad_image(im,s)
if size(im,1)<s
    p1=s-size(im,1);
    p2=s-size(im,2);
    im=cat(1,zeros(floor(0.5*p1),size(im,2),size(im,3),size(im,4)),im,zeros(ceil(0.5*p1),size(im,2),size(im,3),size(im,4)));
    im=cat(2,zeros(size(im,1),floor(0.5*p2),size(im,3),size(im,4)),im,zeros(size(im,1),ceil(0.5*p2),size(im,3),size(im,4)));
end
end