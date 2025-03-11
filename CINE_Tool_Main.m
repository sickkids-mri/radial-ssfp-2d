function CINE_Tool_Main( workDirPath, taskFileName, outputDirPath, modeFilePath )

origPath = path;
addpath(genpath('utilities'))
resetPath = onCleanup(@()path(origPath));


%% Identify and Verify Files

 % Task File
 taskFilePath = fullfile( workDirPath, taskFileName );
 assert( exist( taskFilePath, 'file' ) == 2, sprintf('File does not exist: %s',taskFilePath) )
 taskParam = yarra_read_ini_file( taskFilePath );
 
 % Scan File
 scanFileName = taskParam.task.scanfile;
 scanFilePath = fullfile( workDirPath, scanFileName );
 assert( exist( scanFilePath, 'file' ) == 2, sprintf('File does not exist: %s',scanFilePath) )
 
 % Mode File
 assert( exist( modeFilePath, 'file' ) == 2, sprintf('File does not exist: %s',modeFilePath) )
 modeParam = yarra_read_ini_file( modeFilePath );

 
%% Load Raw Data

fprintf( 'Loading %s...\n', scanFileName );

fprintf('Reading raw data %s \n', scanFilePath);
Path = scanFilePath;
[KSpace,Times,Noise,Header]=CINE_Tool_Read_Raw_Data(Path);

%% Run Reconstruction

%  Lookup Dicom StudyUID and PatientName 
studyLOID = Header{end}.hdr.Meas.StudyLOID;  % raw data study ID 
[ studyUID, patientName ] = lookup_study_info( studyLOID, taskParam.task.patientname );  % get dicom UID and patient name that are consistent across study
fprintf('PatientName       %s\n', patientName );
fprintf('StudyInstanceUID  %s\n', studyUID );

% Number of Cardiac Phases
if ( isfield( taskParam.task, 'paramvalue' ) )
    if ~isempty( taskParam.task.paramvalue )
        if ( taskParam.task.paramvalue >= 5 && taskParam.task.paramvalue <= 35 )
            numCardiacPhase = taskParam.task.paramvalue;
        else
            numCardiacPhase = 15;
        end
    end
end

CINE_Tool_Main_Helper(KSpace,Times,Noise,Header,scanFileName,numCardiacPhase,outputDirPath,patientName,studyUID);

fprintf('Reconstructions finished. \n');

end  % CINE_Tool_Main(...)
