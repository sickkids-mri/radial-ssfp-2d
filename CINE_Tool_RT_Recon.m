function S = CINE_Tool_RT_Recon(param)
M = param.E'*param.y;
param.y=param.y./max(abs(M(:)));
M=M./max(abs(M(:)));

ite=0;
dt=0*M(:,:,1:end-1,:);
tic1=tic;
while(1),
    clear S
    ite=ite+1;
%     display(['iter: ',num2str(ite),' Elapsed Time ',num2str(toc(tic1))])
    dt=dt+0.25*diff(M,1,3);
    dt=dt./abs(dt).*max(min(abs(dt),param.lambda(1)/2),-param.lambda(1)/2);
    dt(isnan(dt))=0;
    S=M-cat(3,-dt(:,:,1,:),-diff(dt,1,3),dt(:,:,end,:));%-z;
    M=S-param.E'*(param.E*S-param.y);
    if (ite > param.nite), break;end
end
end
