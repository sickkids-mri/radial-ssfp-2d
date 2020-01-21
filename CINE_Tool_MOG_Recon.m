function S = CINE_Tool_MOG_Recon(param)
M = param.E'*param.y;
param.y=param.y./max(abs(M(:)));
M=M./max(abs(M(:)));

ite=0;
% dt=0*M(:,:,1:end-1,:);
dt=zeros(size(cat(3,M(:,:,1),M)));
tic1=tic;
while(1),
    clear S
    ite=ite+1;
%     dt=dt+0.25*diff(M,1,3);
    dt=dt+0.25*diff(cat(3,M(:,:,end),M,M(:,:,1)),1,3);
    dt=dt./abs(dt).*max(min(abs(dt),param.lambda(1)/2),-param.lambda(1)/2);
    dt(isnan(dt))=0;
    S=M+diff(dt,1,3);%cat(3,-dt(:,:,1,:),-diff(dt,1,3),dt(:,:,end,:));%-z;
    M=S-param.E'*(param.E*S-param.y);
%     display(['iter: ',num2str(ite),' Elapsed Time ',num2str(toc(tic1))])

    if (ite > param.nite), break;end
end
end
