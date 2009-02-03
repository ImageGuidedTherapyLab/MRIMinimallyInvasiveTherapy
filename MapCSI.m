
%First you need to load in the images you want to analyze.  Open the
%directory where the images are located.

img=readseriesDCM([1 32],'');   % note that I placed 32 here.  Replace that
                               % with the number of images in the folder.
                               % It must be a multiple of 32!

%The following three commands formats the image into complex form (real and
%imaginary numbers for each pixel.

img=double(img);
dat=cplximg(img);
dat2=reshape(dat,[256 256 16 1]);  %Replace 1 with the number of images divided
by 32.
                                  %For instance, if you have 320 images
                                  %delete "1" and type in "10"

%This command extracts the times the echoes were collected (TE values)
jj=1;
for ii=1:2:32,
  [tmp1 tmp2]=readseriesDCM([ii ii],'');
  te(jj)=tmp2.EchoTime;
  jj=jj+1;
end

%This tells Matlab how many acquisitions we did.  ONLY RUN IF THER IS MORE
%THAN ONE ACQUISITION.  IF THERE IS ONLY ONE, THEN SET s=1;
a=size(dat2);
s=a(4);

%Now we will create blank images where we will input our results into.
warning off
t2star=zeros(256,256,s);
freq=zeros(256,256,s);
amp=zeros(256,256,s);
amp0=zeros(256,256,s);

%Define gamma-B0 and the echo spacing
gB0=1.5*42.58;
esp=te(2)-te(1);


%It is time to run the spectral analysis algorithm
for z=1:s
   disp(['Processing Acquisition ' int2str(z)])
   for x=1:256;
      for y=1:256;
         if abs(dat2(x,y,2,z))>50;   % Only processes data with signal
           % Timing-correction for first echo
           dat2(x,y,1,z)=dat2(x,y,1,z).*exp(-2*pi*i*31.25e3*4e-6);  
           avg1=squeeze(dat2(x,y,:,z));
           [p q]=stmcb(avg1,1,2);
           rtz=roots(q);
           ppm=imag(log(rtz))/2/pi/(1e-3*esp*gB0);
           ppm1(x,y,z)=ppm(1);
           ppm2(x,y,z)=ppm(2);
           t2star=-1./real(log(rtz)/esp)/1000;
           amp1(x,y,z)=((p(1).*rtz(1)+p(2))/((2.*q(1).*rtz(1))+q(2)));
           amp2(x,y,z)=((p(1).*rtz(2)+p(2))/((2.*q(1).*rtz(2))+q(2)));
           t2star1(x,y,z)=(t2star(1));
           t2star2(x,y,z)=(t2star(2));
         end
      end
   end
end

