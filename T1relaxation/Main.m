% loading measurement data, Note rhodata is calculated by dividing the
% corresponding data from two different flip angle, for detail see the
% following:
% MATLAB/projects/mfgreScripts/ooMFGRE
load rhodata    % data directory: /mnt/FUS4/data2/madankan/FlipAngle

% initializing variable T1ref at each spatio-temporal location: The first
% two elements specify the spatial location, and the third array assigns
% the time step!
T1ref=zeros(256,256,60);
% initializing variable slope at each spatio-temporal location: The first
% two elements specify the spatial location, and the third array assigns
% the time step!
slope=zeros(256,256,60);
ctr=0;
% looping through (x,y) spatial points, just ROI is being studied!
matlabpool open
for ix=110:180 % over x direction
    parfor iy=110:180 % over y direction
        if ~isnan(sum(Temp(ix,iy,:))) % avoid possible singularities
            % Perform minimum variance estimation to estimate coefficients
            % T1ref and slope (ZP is the concatenated vector of these), PP
            % represents corresonding covariance matrices at each time
            % step!
            [ZP,PP]=FlipAngleMinVar(5,[700,0],[1500,0.03],ix,iy);
            % assigning the corresponding values of T1ref and slope
            T1ref(ix,iy,:)=ZP(:,1); 
            slope(ix,iy,:)=ZP(:,2);           
        end
        ctr=ctr+1;
    end
end
% calculating the value of T1 from the linear model:
T1=T1ref.*(1+slope.*Temp(:,:,1:nFrames));
% saving the data
save -v7.3 estdata.mat T1 slope T1ref
% make a movie to show variations of T1 field over time:
nFrames=60;
mov(1:nFrames) = struct('cdata', [],'colormap', []);
imagesc(T1(:,:,1),[1000,2000]);
title('time = 0');
set(gca,'fontsize',13)
colorbar
set(gca,'nextplot','replacechildren');
for k = 1:nFrames 
   imagesc(T1(:,:,k),[500,2000]);colorbar
   title(['time = ',num2str(k),' Tr'])
   set(gca,'fontsize',13)
   mov(k) = getframe(gcf);
end
% save AVI file.
movie2avi(mov, 'T1vsTemp.avi', 'compression', 'None');
