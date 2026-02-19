function [KSpace,Times,Noise,Header]=load_raw_data(d,filelist)
    Header = mapVBVD(d);
    Noise=[];
    try Noise = Header{1,1}.noise();end %#ok<*TRYNC>
    try KSpace = Header.image();end
    try KSpace = Header{1,length(Header)}.image();end
    if ~exist('KSpace','var')|| isempty(KSpace)
        KSpace = Header{1,length(Header)}.image.unsorted();
        KSpace=reshape(KSpace,[size(KSpace,1),size(KSpace,2),size(KSpace,3)/Header{1,length(Header)}.image.sqzSize(4),Header{1,length(Header)}.image.sqzSize(4)]);
    end
    KSpace=squeeze(KSpace);
    KSpace=permute(KSpace,[1,3,2,4,5]);
    
    if ~iscell(Header);temp{1,1}=Header;Header=temp;clear temp;end % If only one file is selected, change to cell array for consistency
    
    
    Times=zeros(max(Header{1,length(Header)}.image.Lin(:)),max(Header{1,length(Header)}.image.Sli(:)),max(Header{1,length(Header)}.image.Ave(:)),max(Header{1,length(Header)}.image.Rep(:)),'double');
    for i=1:length(Header{1,length(Header)}.image.Lin)
        Times(Header{1,length(Header)}.image.Lin(i),Header{1,length(Header)}.image.Sli(i),Header{1,length(Header)}.image.Ave(i),Header{1,length(Header)}.image.Rep(i))=2.5*Header{1,length(Header)}.image.timestamp(i);
    end
end
    
    