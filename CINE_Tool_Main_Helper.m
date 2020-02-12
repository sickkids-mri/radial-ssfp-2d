function CINE_Tool_Main_Helper(KSpace,Times,Noise,Header,FileName,outputDirPath,Path)

% Chris Roy 2018
    fprintf('Output path: %s \n', outputDirPath);
    outputFileNamePrefix = strrep( FileName, 'meas_', '' );
    outputFileNamePrefix = strrep( outputFileNamePrefix, '.dat', '' );
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate trajectory and ramp weighting
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     if ~contains(Header{1,length(Header)}.hdr.Config.SequenceFileName(15:end),'ty')
    if ~sum(strfind(Header{1,length(Header)}.hdr.Config.SequenceFileName(15:end),'ty')) % revert to above statement for Matlab Release 2018+
        k = CINE_Tool_Trajectory(size(KSpace,1), size(KSpace,2),1); % For standard golden angle
    else
        k = CINE_Tool_Trajectory(size(KSpace,1), size(KSpace,2),7); % For tiny golden angle
    end
    
    w = CINE_Tool_Weighting(size(KSpace,1));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Debuging: Reduce number of slices or spokes
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %{
    indSlice=round(size(KSpace,4)/2)-1+(1:2);
    indSpoke=max(1,size(KSpace,2)-1200):size(KSpace,2);
    KSpace=KSpace(:,indSpoke,:,indSlice);
    Times=Times(indSpoke,indSlice);
    k=k(:,indSpoke);
    %}
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Remove noisey channels
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~isempty(Noise)
        Noise=squeeze(mean(mean(mean(abs(Noise),1),3),6));
        Coils=1:length(Noise);
        nThresh=2;
        while sum(Noise>nThresh*mean(Noise(~isnan(Noise))))
            Coils(Noise>nThresh*mean(Noise(~isnan(Noise))))=nan;
            Noise(Noise>nThresh*mean(Noise(~isnan(Noise))))=nan;
        end
        KSpace=KSpace(:,:,Coils(~isnan(Coils)),:);
    end
    clear Noise Coils
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Perform coil compression to speed up reconstructions
    % We compress to 8 coils here somewhat arbitrarily. Could be improved
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    try % If there is enough memory this can be done in one step for all slices
        [~,~,V] = svd(reshape(permute(KSpace,[1,2,4,3]),size(KSpace,1)*size(KSpace,2)*size(KSpace,4),size(KSpace,3)),'econ');
        KSpace=permute(reshape(reshape(permute(KSpace,[1,2,4,3]),size(KSpace,1)*size(KSpace,2)*size(KSpace,4),size(KSpace,3))*V,size(KSpace,1),size(KSpace,2),size(KSpace,4),size(KSpace,3)),[1,2,4,3]);
        KSpace=KSpace(:,:,1:8,:);
        clear V
    catch % if not, perform coil compression slice by slice
        for iSlice=1:size(KSpace,4)
            [~,~,V] = svd(reshape(KSpace(:,:,:,iSlice),size(KSpace,1)*size(KSpace,2),size(KSpace,3)),'econ');
            KSpace(:,:,:,iSlice)=reshape(reshape(KSpace(:,:,:,iSlice),size(KSpace,1)*size(KSpace,2),size(KSpace,3))*V,size(KSpace,1),size(KSpace,2),size(KSpace,3));
        end
        KSpace=KSpace(:,:,1:8,:);
        clear 
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Estimate coil sensititives from image reconstructions using all acquires spokes
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
     FT=MCNUFFT(k,repmat(w,[1,size(k,2)]).^2,ones(0.5*size(KSpace,1),0.5*size(KSpace,1)));% create nufft operator
     STACKS=zeros(0.5*size(KSpace,1),0.5*size(KSpace,1),size(KSpace,3),size(KSpace,4),'single'); %Initialize variable
     for iSlice=1:size(KSpace,4)
         for iCoil=1:size(KSpace,3)
             STACKS(:,:,iCoil,iSlice)=FT'*KSpace(:,:,iCoil,iSlice); %Reconstruct images using all spokes for each slice and coil
         end
     end
     clear FT
    csm=zeros(size(STACKS),'single'); %Initialize variable
    try % use parallel processing to speed up estimations of coil sensititives
        parfor iSlice=1:size(STACKS,4)
            csm(:,:,:,iSlice) = ismrm_estimate_csm_walsh(STACKS(:,:,:,iSlice));
        end
    catch % if parallel processing failes (posibly due to memory) use a regular loop
        for iSlice=1:size(STACKS,4)
            csm(:,:,:,iSlice) = ismrm_estimate_csm_walsh(STACKS(:,:,:,iSlice));
        end
    end
    clear STACKS iSlice iCoil
    %% 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %       Step 2: Create real-time reconstructions for each slice
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [~,kRT] = CINE_Tool_RT_Sort(ones(size(k)),k,15,10); % pre compute re-sorted trajectory common to each slice
    % Note that for CINE_Tool_RT_Sort here, window size (15 spokes) and shared spokes(10) are set according to  Roy et al JCMR 2016
    
    param.E=MCNUFFT(kRT,repmat(w,[1,size(kRT,2),size(kRT,3)]).^2,csm(:,:,:,1)); % pre-compute NUFFT operator that is common to each slice
    RT=zeros(size(csm,1),size(csm,1),size(kRT,3),size(csm,4),'single'); %Initialize variable
    param.lambda=0.5; % temporal total variation regularizaler for compressed sensing reconstructions of real-time data
    
    param.nite=4; % number of iterations for compressed sensing algorithm
    
    fprintf('Step 1 of 6 - Loading Raw Data - Finished \n');
    fprintf('Step 2 of 6 - Reconstructing Real-Time Images \n');

    matFileName = strcat( 'reconstructed_RT_images_', outputFileNamePrefix, '.mat' );
    matFilePath = fullfile( outputDirPath, matFileName );
    fprintf( '   %s\n', matFilePath );

    for iSlice=1:size(csm,4)
        param.E=swap_b1(param.E,csm(:,:,:,iSlice)); % Add the correct coil-sensititivy map for each slice to the NUFFT operator
        param.y=CINE_Tool_RT_Sort(KSpace(:,:,:,iSlice),kRT,15,10); % compute re-sorted KSpace for each slice using the sliding window parameters described above
        RT(:,:,:,iSlice) = CINE_Tool_RT_Recon(param); % perform compressed sensing reconstruction for each slice
        
    end
    save( matFilePath, 'RT', '-v7.3' )
    clear param kRT iSlice
    
    %% 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %       Step 3: Perform motion correction for each slice
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    t=imresize(Times(:,1)-min(Times(:,1)),size(RT,3)/length(Times(:,1))); % create interpolated time vector corresponding to the temporal resolution of real-time reconstructions
    
    [ROIy,ROIx] = CINE_Tool_AutoROI(RT,t); % automaally detect the location of the heart from the stack of real-time images and time vector
    
    [MIy,MIx]=CINE_Tool_Enlarge_ROI(ROIy,ROIx,round(size(RT,1)/2),round(size(RT,1)/2));% % create an ROI of roughly 1/2 the field-of-view centered on the heart
    
    if max(MIy)>size(RT,1);MIy=MIy+(max(MIy)-size(RT,1));end % make sure the expanded ROI does not exceed the image size
    if max(MIx)>size(RT,1);MIx=MIx-(max(MIx)-size(RT,1));end % make sure the expanded ROI does not exceed the image size
    
    AVGFrames=round(range(t)/200); % Calculate the number of frames needed for motion correction of the real time images (temporal resolution here set to 200 somewhat arbitrarily)
    
    RTAVG=CINE_Tool_Temporal_Averaging(RT,AVGFrames); % Average real-time frames to achieve the temporal resolution described in the previous line
    
    MI = CINE_Tool_Mutual_Information(Decimate_CINE(RTAVG(MIy,MIx,:,:),4)); % Calculate the mutual information of all frames within the 1/2 FOV ROI described above
    % This measure is used to automaally detect and weigh down instances of gross fetal movement
    
    [MOCOy,MOCOx]=CINE_Tool_Enlarge_ROI(ROIy,ROIx,round(size(RT,1)/4),round(size(RT,1)/4));% % create an ROI of roughly 1/4 the field-of-view centered on the heart for motion correction
    
    if max(MOCOy)>size(RT,1);MOCOy=MOCOy+(max(MOCOy)-size(RT,1));end % make sure the expanded ROI does not exceed the image size
    if max(MOCOx)>size(RT,1);MOCOx=MOCOx-(max(MOCOx)-size(RT,1));end % make sure the expanded ROI does not exceed the image size
    
    myTransforms=zeros(AVGFrames,2,size(RT,4));
    
    fprintf('Step 2 of 6 - Reconstructing Real-Time Images - Finished \n')
    fprintf('Step 3 of 6 - Calculating Translational Motion \n');

    for iSlice=1:size(RT,4)

        Targets=repmat(RTAVG(:,:,MI(:,iSlice)==max(MI(:,iSlice)),iSlice),[1,1,AVGFrames]); % Choose frame that has the highest mutual information with all frames as a target for motion correction
        myTransforms(:,:,iSlice) = CINE_Tool_MOCO(abs(RTAVG(:,:,:,iSlice)),abs(Targets), MOCOy,MOCOx); % Estimate translational displacement from the lower temporal resolution real-time images
    end
    RT = CINE_Tool_Correct_RT(RT,-myTransforms);% Apply motion estimations to real-time reconstructions
    
    matFileName = strcat( 'calculated_translational_motion_', outputFileNamePrefix, '.mat' );
    matFilePath = fullfile( outputDirPath, matFileName );
    fprintf( '   %s\n', matFilePath );
    save( matFilePath, 'RT', '-v7.3' )
    
    clear RTAVG MOCOy MOCOx AVGFrames iSlice Targets
    %% Save Data to .mat File

    matFilePath = fullfile( outputDirPath, matFileName );

    %% 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %       Step 4: Perform metric optimized gating to estimate the heart rate for each slice
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [ROIy,ROIx,Initial_RRs] = CINE_Tool_AutoROI(RT,t);% update the location of the heart using the motion corrected real-time images and estimate the mean heart rate
    
    [MIy,MIx]=CINE_Tool_Enlarge_ROI(ROIy,ROIx,round(size(RT,1)/2),round(size(RT,1)/2));% % create an ROI of roughly 1/2 the field-of-view centered on the heart
    if max(MIy)>size(RT,1);MIy=MIy+(max(MIy)-size(RT,1));end % make sure the expanded ROI does not exceed the image size
    if max(MIx)>size(RT,1);MIx=MIx-(max(MIx)-size(RT,1));end % make sure the expanded ROI does not exceed the image size
    
    MI = CINE_Tool_Mutual_Information(Decimate_CINE(RT(MIy,MIx,:,:),4));% update mutual information of all frames using motion corrected real-time
    clear MIy MIx
    MOG=struct('RWaveTimes',[]);% Create structure that will keep the RWaveTimes determined by MOG
    fprintf('Step 3 of 6 - Calculating Translational Motion - Finished \n');
    fprintf('Step 4 of 6 - Calculating Heart Rate \n');
    for iSlice=1:size(RT,4)
        MOG(iSlice,1)=CINE_Tool_MOG(imgaussfilt(abs(RT(ROIy,ROIx,:,iSlice)),2),t,[],20,MI(:,iSlice),Initial_RRs(iSlice));% Run MOG to determine the heart-rate for each slice
    end
    clear RT Initial_RRs iSlice t ROIy ROIx
    
    %% 5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %       Step 5: Reconstruct Final motion-corrected and metric optimized gated CINE Images using compressed sensing
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    rFrames=20;
    CINEs=zeros(0.5*size(KSpace,1),0.5*size(KSpace,1),rFrames,size(KSpace,4),'single'); %Initialize variable
    fprintf('Step 4 of 6 - Calculating Heart Rate - Finished \n');
    fprintf('Step 5 of 6 - Reconstructing CINE Images \n');
   
    matFileName = strcat( 'reconstructed_CINE_images_', outputFileNamePrefix, '.mat' );
    matFilePath = fullfile( outputDirPath, matFileName );
    fprintf( '   %s\n', matFilePath );
    
    for iSlice=1:size(KSpace,4)
    
        CP = Calculate_CardiacPhases(Times(:,iSlice)-min(Times(:,iSlice)),MOG(iSlice).RWaveTimes);% Calculate caridac phase for each slice according to RWaveTimes from MOG
        MOCO=exp(-2*pi*1i*(-bsxfun(@times,imresize(-myTransforms(:,2,iSlice),[size(KSpace,2),1])',real(k))-bsxfun(@times,imresize(-myTransforms(:,1,iSlice),[size(KSpace,2),1])',imag(k)))); %Calculate
        % phase change to be applied to K-space for each slice to correct for the estimated motions
        
        W = KSpace_Interpolation_Weights(size(KSpace,1),CP,MI(:,iSlice),rFrames); % Determine weightings which incorporate the cardiac phase and mutual information for each spoke in order to create motion corrected CINEs
        
        [~,SpokeIndices]=sort(squeeze(W(1,:,:)),1,'descend');rSpokes=max(sum(squeeze(W(1,:,:))>max(W(:))/1000)); % Sort spokes according to weights
        % and calculate the total number of spokes required per frame based on exclusion of spokes that are <0.1%of the max weight
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Re-sort k-space according to the weighting matrix (W) and
        % reconstructed spokes per frame (rSpokes)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        kMOG=zeros(size(KSpace,1),rSpokes,rFrames,'single'); %Initialize variable
        param.y=zeros(size(KSpace,1),rSpokes,size(KSpace,3),rFrames,'single'); %Initialize variable
        for iFrame=1:rFrames
            param.y(:,:,:,iFrame)=bsxfun(@times,bsxfun(@times,KSpace(:,SpokeIndices(1:rSpokes,iFrame),:,iSlice),MOCO(:,SpokeIndices(1:rSpokes,iFrame))),W(:,SpokeIndices(1:rSpokes,iFrame),iFrame));
            kMOG(:,:,iFrame)=k(:,SpokeIndices(1:rSpokes,iFrame));
        end
        
        
        param.E=MCNUFFT(kMOG,repmat(w,[1,size(kMOG,2),size(kMOG,3)]).^2,csm(:,:,:,iSlice)); % create NUFFT operator for each slice
        param.lambda=0.1;
        param.nite=9;
        
        param.tol=0.0025;
        
        CINEs(:,:,:,iSlice) = CINE_Tool_MOG_Recon(param);
    end
    save( matFilePath, 'CINEs', '-v7.3' )
    clearvars -except CINEs Header FileName outputDirPath
      
    %% 6 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %       Step 6 - Generate and save DICOMs from the final CINEs
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %     StartTime=; %Keep track of internal clock time to aproximate time remaining for current step
 
   %{ 
    Patient=num2str(1);
    fileInfo = dir(FileName);
    Date = datestr( fileInfo.datenum, 'yyyyMMdd' );
    Time = datestr( fileInfo.datenum, 'hhmmss' );
    
    FrameCounter=0;
    SliceCounter=1;
    SeriesCounter=0;
    info.StudyInstanceUID=dicomuid;
    info.SeriesInstanceUID=dicomuid;
    info.Manufacturer='SIEMENS';
    info.PixelSpacing=[1;1];
    info.Modality = 'MR';
    info.PatientName=Patient;
    info.ImageComments='N/A';
    info.PatientComments='N/A';
    
    SlicePos=Header{1,2}.image.slicePos(:,1:Header{1,2}.image.NLin:end);
    
    fprintf( 'Saving results...\n' )

    outputFileNamePrefix = strrep( FileName, 'meas_', '' );
    dicomFileName = sprintf( '%s.dcm', outputFileNamePrefix );
    dicomFilePath = fullfile( outputDirPath, dicomFileName );
    
    fprintf('Step 5 of 6 - Reconstructing CINE Images - Finished \n');
    fprintf('Step 6 of 6 - Saving DICOMs - Image \n');

    for iDCM=1:size(CINEs,3)*size(CINEs,4)
        
        FrameCounter=FrameCounter+1;
        if FrameCounter>size(CINEs,3)
            FrameCounter=1;
            SliceCounter=SliceCounter+1;
        end
        
        
        I=int16(255*abs(CINEs(:,:,:,SliceCounter)./max(col(CINEs(:,:,:,SliceCounter)))));
        
        info.InstanceNumber=iDCM;
        
        info.SliceLocation=SliceCounter*Header{1,size(Header,2)}.hdr.MeasYaps.sSliceArray.asSlice{1,1}.dThickness;
        
        info.ImagePositionPatient=SlicePos(1:3,SliceCounter);
        
        info.Private_0019_1015=SlicePos(1:3,SliceCounter);
        info.Private_0019_1016=SlicePos(1:3,SliceCounter);
        
        slicePos = Header{1, size(Header,2)}.image.slicePos(4:7,1);
        a = slicePos(1); b = slicePos(2); c = slicePos(3); d = slicePos(4);
        R = [a*a+b*b-c*c-d*d, 2*b*c-2*a*d, 2*b*d+2*a*c;
            2*b*c+2*a*d, a*a-b*b+c*c-d*d, 2*c*d-2*a*b;
            2*b*d-2*a*c, 2*c*d+2*a*b, a*a-b*b-c*c+d*d  ];
        newR = R*[1 0 0; 0 1 0]';
        info.ImageOrientationPatient = newR(:);
        
        
        info.SeriesNumber=SeriesCounter;
        info.StudyDate=Date;
        info.SeriesDate=Date;
        info.AcquisitionDate=Date;
        info.ContentDate=Date;
        info.StudyTime=Time;
        
        info.AcquisitionTime=Time;
        info.SliceThickness=Header{1,size(Header,2)}.hdr.MeasYaps.sSliceArray.asSlice{1,1}.dThickness;
        info.Filename=FileName;

        dicomwrite(I(:,:,FrameCounter),[outputDirPath,'_Slice_',num2str(SliceCounter),'_Frame_',num2str(FrameCounter),'.dcm'],info);
    end
    %}
    clearvars 

    %fprintf('Step 6 of 6 - Saving DICOMs - Image - Finished \n');
    
end
