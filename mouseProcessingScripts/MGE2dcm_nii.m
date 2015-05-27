function [ ] = MGE2dcm_nii( inMGE,  param, writePath)
%MGE2dcm_nii Creates DICOM and nifti files from MGE object. Geometry should
%be preserved for DICOM and nifti reference respectively but not between
%them
%   Detailed explanation goes here

%% Determine which data to write and extract from object
if strcmp(param,'mag');
    dataMat=squeeze(abs(inMGE.cdata(:,:,:,1,:)));
elseif strcmp(param,'geom');
    inMGE.cropTimePoints([1:1]);
    dataMat=squeeze(abs(inMGE.cdata(:,:,:,1,:)));
else
    dataMat=inMGE.(param);
end

%% Write processed DICOMS
pDICOMdir=writePath;

%pDICOMdir=[regexprep(pDICOMdir,'e\w\w\w/s\w\w\w\w',''),'examtype/' param];
clearMakeDir(pDICOMdir)

for ii=1:inMGE.numSlices;
    for jj=1:inMGE.numTimePoints;
        
        tmpim=dataMat(:,:,ii,jj);
        
        %rescaleSlope=range(tmpim(:)./2^15-1);
        %rescaleInt=min(tmpim(:));
        
        %tmpim2=int16((tmpim-rescaleInt)/rescaleSlope);
        
        metaStruct=inMGE.hdrEx;
        %metaStruct.Private_0028_1053=rescaleSlope; %Slope
        %metaStruct.Private_0028_1052=rescaleInt; %Intercept
        metaStruct.ScannerTableEntry=[];%num2str(metaStruct.ScannerTableEntry);
        metaStruct.PatientPosition=inMGE.patientPosition;
        metaStruct.ImageOrientationPatient=inMGE.sliceOrientation;
        metaStruct.ImagePositionPatient=inMGE.sliceCorner(ii,:);
        metaStruct.InstanceNumber=jj;
        metaStruct=rmfield(metaStruct,'ReferencedImageSequence');
        metaStruct=rmfield(metaStruct,'SeriesFromWhichPrescribed');
        
        dicomwrite(tmpim(:,:,1), [pDICOMdir '/sl' sprintf('%.3i',ii) 'ph' sprintf('%.3i',jj) '.dcm'],metaStruct,'WritePrivate',0,'CreateMode','Copy');
        %dicomwrite(tmpim(:,:,1), [pDICOMdir '/sl' sprintf('%.3i',ii) 'ph' sprintf('%.3i',jj) '.dcm'],metaStruct,'WritePrivate',1,'CreateMode','Copy');
    end
end

%%

niftiDir=strrep(pDICOMdir,'pDICOM','nifti');
clearMakeDir(niftiDir)
clearMakeDir([niftiDir '/tmp'])

for jj=1:inMGE.numTimePoints
    for ii=1:inMGE.numSlices
        eval(['!cp ' pDICOMdir '/sl' sprintf('%.3i',ii) 'ph' sprintf('%.3i',jj) '.dcm ' niftiDir '/tmp/']);
    end
    eval(['!DicomSeriesReadImageWrite2 ' niftiDir '/tmp ' niftiDir '/ph' sprintf('%.3i',jj) '.nii.gz']);
    
    clearMakeDir([niftiDir '/tmp'])
    
    if ~(strcmp(param,'mag')||strcmp(param,'mag'))
        tmpim=load_untouch_nii([niftiDir '/ph' sprintf('%.3i',jj) '.nii.gz']);
        tmpim.hdr.dime.datatype=64;
        tmpim.hdr.dime.bitpix=64;
        tmpim.img=double(zeros(size(tmpim.img)));
        for kk=1:inMGE.numSlices
            tmpim.img(:,:,kk)=squeeze(flipud(rot90(dataMat(:,:,kk,jj))));
        end
        save_untouch_nii(tmpim, [niftiDir '/ph' sprintf('%.3i',jj) '.nii.gz'])
    end
end
rmdir([niftiDir '/tmp'])
% 
% eval(['!DicomSeriesReadImageWrite2 ' niftiDir ' strrep(dcmFolder,'tempDICOM','tmpdcm.nii.gz') ]);


