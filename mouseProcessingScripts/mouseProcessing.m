%% Read in meta data

%cd('/FUS4/data2/CJM/gitMATLAB/projects/rtwostarTemperature'); %Move to project directory
cd('/FUS4/data2/CJM/MRIMinimallyInvasiveTherapy/'); %Move to project directory


[~,~,c]=xlsread('metadata.xlsx'); % Read excel file

numdatasets=1; %Start counter for number of imported datasets
for ii=3:numel(c(:,1)) % Loop over rows starting at row 3 (first real data)
    if c{ii,1}==1 % If marked for processing loop over columns
        for jj=1:numel(c(1,:))
            switch c{2,jj} % If marked as vector metadata convert to vector otherwise read normally
                case 'vec'
                    metaDatStruct(numdatasets).(c{1,jj})=str2num(c{ii,jj});
                otherwise
                    metaDatStruct(numdatasets).(c{1,jj})=(c{ii,jj});
            end
        end
        numdatasets=numdatasets+1; % Increment valid dataset counter
    end
end

clear c numdatasets ii jj
ii=1;
%%
for ii=1:1
    
tic
    
% [mrtiData] = processMRTI( metaDatStruct(ii) );

load([metaDatStruct(ii).baseFolder '/MRTImats/' metaDatStruct(ii).studyName '.mat'])
mrtiData=processedData;

toc

writeDir1=[metaDatStruct(ii).baseFolder '/pDICOM/' metaDatStruct(ii).studyName '/MRTI'];
MGE2dcm_nii( mrtiData, 'deltaT', writeDir1);

toc

writeDir2=[metaDatStruct(ii).baseFolder '/pDICOM/' metaDatStruct(ii).studyName '/geom'];
MGE2dcm_nii( mrtiData, 'geom', writeDir2 );

toc

dicomSeg=[metaDatStruct(ii).baseFolder '/segmentations/' metaDatStruct(ii).studyName '/DCM'];
niftiSeg=[metaDatStruct(ii).baseFolder '/nifti/' metaDatStruct(ii).studyName '/Seg.nii.gz'];

eval(['!DicomSeriesReadImageWrite2 ' dicomSeg ' ' niftiSeg]);


dataFile=['/FUS4/data2/CJM/cabiGrant/20150523/segmentations/' metaDatStruct(ii).studyName,'/' metaDatStruct(ii).studyName 'Seg.nii'];
dataFile=niftiSeg;

geomFile=['/FUS4/data2/CJM/cabiGrant/20150523/nifti/' metaDatStruct(ii).studyName '/geom/ph001.nii.gz'];
writeDir=[strrep(strrep(dataFile,'segmentations','nifti'),[metaDatStruct(ii).studyName 'Seg'],[metaDatStruct(ii).studyName 'SegReformat'])]
resampleCommand='/opt/apps/ANTsR/dev//ANTsR_src/ANTsR/src/ANTS/ANTS-build//bin/WarpImageMultiTransform';
identityTrans='/FUS4/data2/CJM/MRIMinimallyInvasiveTherapy/rtwostarTemperature/identity.txt';

eval(['!' resampleCommand ' 3 ' dataFile ' ' writeDir ' ' identityTrans ' -R ' geomFile])
toc
end

%% Temperature data
% Process MRTI data using MGE object
for ii=1:7
[mrtiData] = processMRTI( metaDatStruct(ii) );
end
 
% Write processed DICOMS and nifti for temperature data
writeDir=[metaDatStruct(ii).baseFolder '/pDICOM/' metaDatStruct(ii).studyName '/post/MRTI'];
MGE2dcm_nii( mrtiData, 'deltaT', writeDir);

%% Geometric data 
% Write geometric data for reformats
writeDir=[metaDatStruct(ii).baseFolder '/pDICOM/' metaDatStruct(ii).studyName '/post/geom'];
MGE2dcm_nii( mrtiData, 'geom', writeDir );

dataFile=['/mnt/FUS4/data2/CJM/cabiGrant/20150523/Segmentations/' metaDatStruct(ii).studyName 'Seg']

%% T2w Pre
% Write OrchDICOMs from FIESTAC
processT2w( metaDatStruct(ii), 'preT2w' )
dataFile='/FUS4/data2/CJM/cabiGrant/20150404/nifti/M1/post/T2w/t2w.nii.gz';
geomFile='/FUS4/data2/CJM/cabiGrant/20150404/nifti/e231/s2793/geom/ph001.nii.gz';
resampleCommand='/opt/apps/ANTsR/dev//ANTsR_src/ANTsR/src/ANTS/ANTS-build//bin/WarpImageMultiTransform';
identityTrans='/FUS4/data2/CJM/gitMATLAB/projects/rtwostarTemperature/identity.txt';

eval(['!' resampleCommand ' 3 ' dataFile ' ' strrep(dataFile,'t2w','t2wR') ' ' identityTrans ' -R ' geomFile])

%% T1w Pre
processT1w( metaDatStruct(ii), 'preT1w' )

processT1w( metaDatStruct(ii), 'postT1w' )

%%
dataFile='/FUS4/data2/CJM/cabiGrant/20150404/segmentations/M1/post/tumorSegmentation.nii';
eval(['!' resampleCommand ' 3 ' dataFile ' ' strrep(dataFile,'segmentations','nifti') ' ' identityTrans ' -R ' geomFile])





