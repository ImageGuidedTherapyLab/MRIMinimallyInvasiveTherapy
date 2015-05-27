function [processedData] = processTemp( inStruct )
%UNTITLED32 Summary of this function goes here
%   Detailed explanation goes here
tempDir=[inStruct.baseFolder '/DICOM/' inStruct.MRTI];
processedData=MGE(tempDir);

processedData.threshold=inStruct.MRTI_threshold;
processedData.setProcessingROI(inStruct.MRTI_pROI);
processedData.setDriftROI(inStruct.MRTI_dROI);
processedData.firstImage=inStruct.MRTI_baselineImage;

MGE2dcm_nii( tmp, 'deltaT' )




