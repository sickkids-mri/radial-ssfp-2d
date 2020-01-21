function RTMI = Real_Time_Mutual_Information(RT)
% X=10*permute(INorm(double(Crop_CINE(RT,size(RT,1)/2))),[3,1,2]);X=X(:,:);
RTMI=zeros(size(RT,3),size(RT,4),'single');
for iSlice=1:size(RT,4)
    X=10*permute(INorm(double(RT(:,:,:,iSlice))),[3,1,2]);X=X(:,:);
    %     MM=zeros(size(X,1),size(X,1));
    [Y,Z]=meshgrid(1:size(RT,3),1:size(RT,3));
    IND=find(triu(Y,1)>0);
    IND=[Y(IND),Z(IND)];
    mm=zeros(length(IND),1);
    parfor i=1:length(IND)
        mm(i,1) = mymutualinfo(X,IND(i,:));
    end
    MM=zeros(size(X,1),size(X,1));
    IND=find(triu(Y,1)>0);
    for i=1:length(IND)
        MM(IND(i))=mm(i);
    end
%     mm(:,iSlice)=MM+MM';
    
    MM=mean(MM+MM');
%     MM=MM./max(MM);
    RTMI(:,iSlice)=MM;
end
end

function M = mymutualinfo(X,IND)
M = mutualinfo(abs(X(IND(1),:)),abs(X(IND(2),:)));
end


% function MM = Real_Time_Mutual_Information(RT)
% Case=2;
% if Case == 1
%     ERMS=squeeze(sqrt(sum(sum(abs(Crop_CINE(RT,size(RT,1)/2)-repmat(mean(Crop_CINE(RT,size(RT,1)/2),3),[1,1,size(RT,3)])).^2))));
%     j=find(ERMS==min(ERMS(:)),1);
%     X=10*permute(INorm(double(Crop_CINE(RT,size(RT,1)/2))),[3,1,2]);X=X(:,:);
%     MM=zeros(size(X,1),1);
%     for i=1:size(X,1)
%         % for j=1:size(X,1)
%         MM(i,1) = mutualinfo(abs(X(i,:)),abs(X(j,:)));
%         %     MM(i,1) = mutualinfo(abs(X(i,:)),abs(mean(X,1)));
%         % end
%     end
% elseif Case == 2
%     X=10*permute(INorm(double(Crop_CINE(RT,size(RT,1)/2))),[3,1,2]);X=X(:,:);
%     MM=zeros(size(X,1),size(X,1));
%     for i=1:size(X,1)
%         for j=1:size(X,1)
%             MM(i,j) = mutualinfo(abs(X(i,:)),abs(X(j,:)));
%         end
%     end
%     MM=mean(MM);
% elseif Case == 3
%     X=10*permute(INorm(double(Crop_CINE(RT,size(RT,1)/2))),[3,1,2]);X=X(:,:);
%     MM=zeros(size(X,1),1);
%     for i=1:size(X,1)-1
%         MM(i,1) = mutualinfo(abs(X(i+1,:)),abs(X(i,:)));
%     end
%     MM(end,1) = mutualinfo(abs(X(end,:)),abs(X(1,:)));
%
% end
% end
