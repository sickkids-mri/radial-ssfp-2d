function [ROIy,ROIx,RRs,FIG,SIGS] = CINE_Tool_AutoROI(RT,t)
% This is a work in progress. Will Probably be replaced by a CNN version.

[~,i]=max(size(RT(:,:,1,1)));
[ROIy,ROIx] = CINE_Tool_Enlarge_ROI(round(size(RT,1)/2),round(size(RT,2)/2),round(size(RT,i)/3),round(size(RT,i)/3));

yMin=ROIy(1);yMax=ROIy(end);xMin=ROIx(1);xMax=ROIx(end);


if ~rem(size(RT,3),2)
    RT=cat(3,RT,zeros(size(RT(:,:,1,:))));
end

if length(t)==1
    f=linspace(-1/(2*t),1/(2*t),size(RT,3));
else
    f=linspace(-1/(2*mean(diff(t(~isnan(t))/1000))),1/(2*mean(diff(t(~isnan(t))/1000))),size(RT,3));
end
iF=[find(f<=-(1/0.6)&f>=-(1/0.33)),find(f>=(1/0.6)&f<=(1/0.33))];


if mean(abs(diff(1000./abs(f(iF)))))>50
    [~,i]=max(size(RT(:,:,1,1)));
    [ROIy,ROIx] = CINE_Tool_Enlarge_ROI(round(size(RT,i)/2),round(size(RT,i)/2),round(size(RT,i)/4),round(size(RT,i)/4));
    if min(ROIy)<0,ROIy=ROIy-min(ROIy(:))+1;end;if min(ROIx)<0,ROIx=ROIx-min(ROIx(:))+1;end
    if max(ROIy)>size(RT,1),ROIy=ROIy(1):size(RT,1);end;if max(ROIx)>size(RT,2),ROIx=ROIx(1):size(RT,2);end


    if isempty(ROIy), ROIy=1:size(RT,1);end
    if isempty(ROIx), ROIx=1:size(RT,2);end

    
    [Sy,Sx] = CINE_Tool_Enlarge_ROI(ROIy,ROIx,(0.5*size(RT,1)),0.5*size(RT,1));
    FIG=Draw_ROI(INorm_CINE(squeeze(mean(abs(RT),3))),ROIy,ROIx);
    if min(Sy)<0,Sy=Sy-min(Sy(:))+1;end;if min(Sx)<0,Sx=Sx-min(Sx(:))+1;end
    if max(Sy)>size(RT,1),Sy=Sy(1):size(RT,1);end;if max(Sx)>size(RT,2),Sx=Sx(1):size(RT,2);end
    if isempty(Sy), Sy=ROIy;end
    if isempty(Sx), Sx=ROIx;end
    
    clear RT
    FIG=CINE_Array(FIG(Sy,Sx,:));
    SIGS=[];
    RRs=[];
else



    iFF=1:size(RT,3);iFF(ceil(length(iFF)/2))=nan;
    iFF(iF)=nan;iFF=iFF(~isnan(iFF));
    
    % F=abs(fftc(imgaussfilt(Window_CINE(abs(RT),0,0.5),2),3));
    F=abs(fftc(RT,3));
    FF=F(:,:,iFF(1:length(iFF)/2),:)+flip(F(:,:,iFF(1+length(iFF)/2:end),:),3);
    F=F(:,:,iF(1:length(iF)/2),:)+flip(F(:,:,iF(1+length(iF)/2:end),:),3);
    
    F=bsxfun(@rdivide,F,max(FF,[],3));
    F(isnan(F(:)))=0;
    F=mygaussfilt(F,4);
    
    
    MI = CINE_Tool_Mutual_Information(Decimate_CINE(Crop_CINE(RT,size(RT,1)/2),4));
    iter=0;
    clear FIG
    while iter<50
        iter=iter+1;
        % display(iter)
        Freqs=zeros(size(RT,4),1);
        SIGS=zeros(size(squeeze(RT(:,:,1,:))));
        
        C=zeros(size(F,3),size(F,4));
        for i=1:size(F,3)
            for j=1:size(F,4)
                try C(i,j)=imagemetric(abs(F(ROIy,ROIx,i,j)),{'CINE'});end
            end
        end
        
        FreqCost=zeros(1,size(F,4));
        [~,i]=min(C);
        ii=1:size(F,3);
        for j=1:size(F,4)
            FreqCost(1,j)=range(C(:,j))./std(C(ii~=i(j),j));
        end
        FreqCost(isnan(FreqCost))=1;
        
        for iSlice=1:size(RT,4)
            Freqs(iSlice,1)=find(C(:,iSlice)==min(C(:,iSlice)),1,'first');
            SIGS(:,:,iSlice)=F(:,:,C(:,iSlice)==min(C(:,iSlice)),iSlice);
        end
        
        % MotionFilter=permute(sum(MI)./sum(MI(:)),[3,1,2]);
        % % CostFilter=permute(1./FreqCost,[3,1,2]);
        % CostFilter=permute(FreqCost,[3,1,2]);
        % CostFilter=CostFilter./sum(CostFilter(:));
        % CostFilter=ones(size(CostFilter));
        
        MotionFilter=sum(MI)./sum(MI(:));
        CostFilter=FreqCost./sum(FreqCost(:));
        
        
        % S=Crop_and_Pad(imgaussfilt(SIGS,1),ROIy,ROIx);
        S=Crop_and_Pad(SIGS,ROIy,ROIx);
        MaxPoints=zeros(4,size(S,3));
        for iSlice=1:size(S,3)
            [i,j]=find(S(:,:,iSlice)==max(col(S(:,:,iSlice))));
            % MaxPoints(:,iSlice)=[i,j,max(col(S(:,:,iSlice)))];
            MaxPoints(:,iSlice)=[i,j,max(col(S(:,:,iSlice))),std(col(S(:,:,iSlice)))];
        end
        % MaxPoints(3,:)=MaxPoints(3,:)./sum(MaxPoints(3,:));
        
        S=sum(S,3);
        [i,j]=find(S==max(S(:)),1,'first');
        
        i(2,1)=round(sum(MaxPoints(1,:).*CostFilter.*MotionFilter.*MaxPoints(3,:))/sum(CostFilter.*MotionFilter.*MaxPoints(3,:)));
        j(2,1)=round(sum(MaxPoints(2,:).*CostFilter.*MotionFilter.*MaxPoints(3,:))/sum(CostFilter.*MotionFilter.*MaxPoints(3,:)));
        
        i(find(abs(i-round(0.5*(ROIy(1)+ROIy(end))))==min(abs(i-round(0.5*(ROIy(1)+ROIy(end))))),1));
        j(find(abs(j-round(0.5*(ROIx(1)+ROIx(end))))==min(abs(j-round(0.5*(ROIx(1)+ROIx(end))))),1));
        
        if i<yMin,i=yMin;end
        if i>yMax,i=yMax;end
        if j<xMin,j=xMin;end
        if j>xMax,j=xMax;end
        
        
        % i=round(sum(MaxPoints(1,:).*CostFilter.*MotionFilter)/sum(CostFilter.*MotionFilter));
        % j=round(sum(MaxPoints(2,:).*CostFilter.*MotionFilter)/sum(CostFilter.*MotionFilter));
        
        
        % S=max(imgaussfilt(bsxfun(@times,Crop_and_Pad(SIGS,ROIy,ROIx),CostFilter.*MotionFilter),16),[],3);
        % S=sum(bsxfun(@times,Crop_and_Pad(SIGS,ROIy,ROIx),CostFilter.*MotionFilter),3);
        
        
        
        [ROIy,ROIx] = CINE_Tool_Enlarge_ROI(i,j,length(ROIy)-4,length(ROIx)-4);
        if min(ROIy)<0,ROIy=ROIy-min(ROIy(:))+1;end;if min(ROIx)<0,ROIx=ROIx-min(ROIx(:))+1;end
        if max(ROIy)>size(RT,1),ROIy=ROIy(1):size(RT,1);end;if max(ROIx)>size(RT,2),ROIx=ROIx(1):size(RT,2);end
        
        
        if length(ROIy)<45&&length(ROIx)<45
%             display(iter)
            break;
        end
    end
    clear F
    
    
    [Sy,Sx] = CINE_Tool_Enlarge_ROI(ROIy,ROIx,(0.5*size(RT,1)),0.5*size(RT,1));
    if min(Sy)<0,Sy=Sy-min(Sy(:))+1;end;if min(Sx)<0,Sx=Sx-min(Sx(:))+1;end
    if max(Sy)>size(RT,1),Sy=Sy(1):size(RT,1);end;if max(Sx)>size(RT,2),Sx=Sx(1):size(RT,2);end
    
    FIG=Draw_ROI(INorm_CINE(squeeze(mean(abs(RT),3))),ROIy,ROIx);clear RT
    
    FIG=CINE_Array(FIG(Sy,Sx,:));
    SIGS=Draw_ROI(INorm_CINE(SIGS),ROIy,ROIx);
    SIGS=CINE_Array(SIGS(Sy,Sx,:));
    
    
    RRs=1000./abs(f(iF(Freqs)))';
    WeightedMeanRRs=(sum(RRs'.*CostFilter.*MotionFilter)/sum(CostFilter.*MotionFilter));
    WeightedSTDRRs=std(1000./abs(f(iF(Freqs)))',squeeze(CostFilter.*MotionFilter));
    RRs(abs(RRs-WeightedMeanRRs)>2*WeightedSTDRRs)=WeightedMeanRRs;
    
    
end
end

