function [ ] = writepDICOM_MGE( inMGE,pDICOM_Folder, param )
%UNTITLED33 Summary of this function goes here
%   Detailed explanation goes here

if strcmp(param,'mag');
    dataMat=squeeze(abs(inMGE.cdata(:,:,:,1,:)));
else
    dataMat=inMGE.(param);
end

for ii=1:inMGE.numSlices;
    for jj=1:inMGE.numTimePoints;
        
        tmpim=dataMat(:,:,ii,jj);
        
        
        %rescaleSlope=range(tmpim(:)./2^15-1);
        %rescaleInt=min(tmpim(:));
        
        tmpim2=int16((tmpim-rescaleInt)/rescaleSlope);
        
        metaStruct=inMGE.hdrEx;
        %metaStruct.Private_0028_1053=rescaleSlope; %Slope
        %metaStruct.Private_0028_1052=rescaleInt; %Intercept
        metaStruct.ScannerTableEntry=num2str(metaStruct.ScannerTableEntry);
        metaStruct.PatientPosition=inMGE.patientPosition;
        metaStruct.ImageOrientationPatient=inMGE.sliceOrientation;
        metaStruct.ImagePositionPatient=inMGE.sliceCorner(ii,:);
        metaStruct.InstanceNumber=jj;
        if jj==1
            dicomwrite(tmpim2(:,:,1),[pDICOM_Folder 'paramFirstTP' '/Sl' num2str(ii), 'Ph' num2str(jj) '.dcm'],metaStruct,'WritePrivate',1,'CreateMode','Copy');
        end
        dicomwrite(tmpim2(:,:,1),[pDICOM_Folder 'param' '/Sl' num2str(ii), 'Ph' num2str(jj) '.dcm'],metaStruct,'WritePrivate',1,'CreateMode','Copy');
    end
end

