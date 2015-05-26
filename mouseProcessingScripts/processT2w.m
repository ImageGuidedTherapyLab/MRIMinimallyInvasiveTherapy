function [ ] = processT2w( inStruct, inString )
%UNTITLED31 Summary of this function goes here
%   Detailed explanation goes here

sepIndx=find(inString=='T');

rawDataPath=[inStruct.baseFolder '/rawData/' inStruct.studyName '/' inString(1:sepIndx-1) '/' 'T2w' ];

tmpdir=dir(rawDataPath);
pfileName=tmpdir(3).name;

%% Write DICOMS using orchestra
reconFIESTAC_CJM([rawDataPath '/' pfileName]);

%% Write nifti File
if strfind(inString,'post')

niftiPath=strrep(rawDataPath,'rawData','nifti');

inPath=strrep([rawDataPath '/' pfileName],'rawData','orchDICOM');

outPath=[niftiPath '/t2w.nii.gz'];

clearMakeDir( niftiPath )

eval(['!DicomSeriesReadImageWrite2 ' inPath ' ' outPath])

end
