function IMG = CINE_Array(I,sz,titleflag)
if length(size(I))==3
    I=permute(I,[1,2,4,3]);
end
if size(I,4)>1
    
    if ~exist('titleflag','var')
        titleflag=0;
    end
    
    if ~exist('sz','var')
        sz=3;
    else
        if isempty(sz)
            sz=3;
        end
    end
    
    
    x=ceil(size(I,4)/sz);
    y=sz;
    IMG=zeros(size(I,1)*sz,size(I,2)*x,size(I,3),'single');
    
    Count=0;
    for j=1:y
        for i=1:x
            Count=Count+1;
            if Count<=size(I,4)
                if titleflag
                    mytitle=(1-text2im(num2str(Count)))*max(makevector(I(:,:,:,Count)));
                    
                    I(1:size(mytitle,1),1:size(mytitle,2),:,Count)=repmat(mytitle,[1,1,size(I,3)]);
                    %                 I(1+size(I,1)-size(mytitle,1):size(I,1),1+size(I,2)-size(mytitle,2):size(I,2),:,Count)=repmat(mytitle,[1,1,size(I,3)]);
                    %                 I(1+size(I,1)-size(mytitle,1):size(I,1),1:size(mytitle,2),:,Count)=repmat(mytitle,[1,1,size(I,3)]);
                end
                IMG(1+size(I,1)*(j-1):j*size(I,1),1+size(I,2)*(i-1):i*size(I,2),:)=I(:,:,:,Count);
            end
        end
    end
else
    IMG=I;
end
end