function [ ] = processT1w( inStruct, inString )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%ims=readDICOMseries([inStruct.baseFolder '/DICOM/' inStruct.(inString)]);

if strfind(inString,'pre')
    niftiPath=[inStruct.baseFolder,'/nifti/' inStruct.studyName '/pre/T1w/'];
else
    niftiPath=[inStruct.baseFolder,'/nifti/' inStruct.studyName '/post/T1w/'];
end

clearMakeDir( niftiPath )

eval(['!DicomSeriesReadImageWrite2 ' [inStruct.baseFolder '/DICOM/'  inStruct.(inString)] ' ' niftiPath 'T1w.nii.gz'])

end

