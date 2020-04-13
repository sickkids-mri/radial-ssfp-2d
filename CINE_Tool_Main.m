function main_func( workDirPath, taskFileName, outputDirPath, modeFilePath )

origPath = path;
addpath(genpath('utilities'))
resetPath = onCleanup(@()path(origPath));


%   jfpva (joshua.vanamerom@sickkids.ca)

%
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


%dir_to_search = pwd;   %current directory

%scanFileName = 'meas_MID00047_FID836853_Radial_STACK_COR_COPYINDX.dat';
%scanFilePath = strcat(dir_to_search, '/', scanFileName);

fprintf( 'Loading %s...\n', scanFileName );

fprintf('Reading raw data %s \n', scanFilePath);
Path = scanFilePath;
[KSpace,Times,Noise,Header]=CINE_Tool_Read_Raw_Data(Path);

CINE_Tool_Main_Helper(KSpace,Times,Noise,Header,scanFileName,outputDirPath,taskParam.task.patientname);

fprintf('Reconstructions finished. \n');

end  


   