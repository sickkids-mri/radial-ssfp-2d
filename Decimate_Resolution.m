function MDR = Decimate_Resolution(z,n)
% z=double(abs(z));
z=double(z);
kern = ones(n,n)*1/n;
MDR = conv2(z,kern);
MDR = MDR(n:n:end, n:n:end);
end