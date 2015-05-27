clear all
close all
format shortg

%% example output:
%% 
%% idx=0 alpha[0]=(1.017979,-3.893608) alpha[2]=(-0.238646,-0.055135) beta[1]=(-2.339487,0.046365) beta[2]=(1.398362,-0.060239)
%% idx=0 Lambda[0]=( 4.93779e-01,-1.12247e-02) Lambda[1]=(-2.83327e+00, 5.75897e-02)
%% idx=0 Amp[0]=(-2.31687e+16, 1.43808e+14) Amp[1]=(-7.69339e+03, 2.22746e+02)
%% Prony Matlab:
%%    2.4507e+05 +          0i   2.0947e+05 +       1106i  -2.8365e+05 +       1344i
%%    2.0947e+05 -       1106i   1.7933e+05 +          0i  -2.4209e+05 +       2422i
%% 
%% 
%% betanormaleqn =
%% 
%%       -2.2229 +   0.020647i
%%        1.2464 -    0.02432i
%% 
%% Prony;Stmcb
%%             1 -          4i     -0.14026 -    0.08795i            1 +          0i      -2.2229 +   0.020647i       1.2464 -    0.02432i
%%        1.1078 -     2.6738i      -0.4314 -     2.9822i            1 +          0i      -2.2278 +    0.01933i       1.2537 -   0.023188i
%% 

% % GPU must be reset on out of bounds errors
% reset(gpuDevice(1))

%%diary
npixel   = 256*256;
nslice   = 5;
necho    = 16;
nspecies = 2;
echospacing = 2.;
imagingfreq = 128.;

% inialized data
h_fid     = complex(repmat(linspace(1,necho,necho)',1,npixel),-repmat((linspace(2,necho+1,necho).^2)',1,npixel));

% perform ONE transfer from host to device
d_fid     = gpuArray( h_fid );

% build the matrix on the device AFTER the data has been transferred
PronyRHS    = reshape( d_fid(nspecies+1:necho  ,:)                             ,necho-nspecies,   1    ,npixel);
PronyMatrix = reshape([d_fid(nspecies  :necho-1,:),d_fid(nspecies-1:necho-2,:)],necho-nspecies,nspecies,npixel);

% batch compute initial beta

d_ppm     = gpuArray( zeros(npixel,nslice,nspecies) );
d_t2star  = gpuArray( zeros(npixel,nslice,nspecies) );
d_T1map   = gpuArray( zeros(npixel,nslice,nspecies) );
d_phase   = gpuArray( zeros(npixel,nslice,nspecies) );

rhs = [
  h_fid( 3);
  h_fid( 4);
  h_fid( 5);
  h_fid( 6);
  h_fid( 7);
  h_fid( 8);
  h_fid( 9);
  h_fid(10);
  h_fid(11);
  h_fid(12);
  h_fid(13);
  h_fid(14);
  h_fid(15);
  h_fid(16);
      ];

matrix = [
 h_fid( 2) ,h_fid( 1);
 h_fid( 3) ,h_fid( 2);
 h_fid( 4) ,h_fid( 3);
 h_fid( 5) ,h_fid( 4);
 h_fid( 6) ,h_fid( 5);
 h_fid( 7) ,h_fid( 6);
 h_fid( 8) ,h_fid( 7);
 h_fid( 9) ,h_fid( 8);
 h_fid(10) ,h_fid( 9);
 h_fid(11) ,h_fid(10);
 h_fid(12) ,h_fid(11);
 h_fid(13) ,h_fid(12);
 h_fid(14) ,h_fid(13);
 h_fid(15) ,h_fid(14);
          ];

normaleqnmatrix =  matrix'*matrix;
normaleqnrhs    = -matrix'*rhs;
disp('Prony Matlab:');
disp([normaleqnmatrix ,normaleqnrhs ]);

betanormaleqn  = normaleqnmatrix \normaleqnrhs    

% get matlab solution
u_in = zeros(necho,1 );
u_in(1) = 1; % make a unit impulse whose length is same as x
[alphaprony betaprony] = prony(h_fid(:,1),nspecies-1,nspecies );
[alphastmcb betastmcb] = stmcbdf(h_fid(:,1),u_in,nspecies-1,nspecies,5,betaprony );
disp('Prony;Stmcb');
disp([alphaprony,betaprony;alphastmcb,betastmcb] );

[alphapronyzero betapronyzero] = prony(h_fid(:,1),0,nspecies );
[ filter( 1, betapronyzero     , h_fid(:,1) ), filter( 1, betaprony         , h_fid(:,1) )]

lambda  = roots([1; betanormaleqn  ]  )
%%  compile and run
StmcbKernelptx = parallel.gpu.CUDAKernel('StmcbKernel.ptx', 'StmcbKernel.cu');
threadsPerBlock = 256;
StmcbKernelptx.ThreadBlockSize=[threadsPerBlock  1];
blocksPerGrid = (npixel  + threadsPerBlock - 1) / threadsPerBlock;
StmcbKernelptx.GridSize=[ceil(blocksPerGrid)  1];
[d_ppm, d_t2star, d_T1map,d_phase, ] = feval(StmcbKernelptx,real(d_fid),imag(d_fid), d_ppm, d_t2star, d_T1map,d_phase, echospacing,imagingfreq,necho,nspecies,npixel);
%% [d_ppm(1:40)', d_t2star(1:40)', d_T1map(1:40)' ,d_phase(1:40)'] 

%% TODO  add total variation regression
%% TODO  %% build CS Kernel
%% TODO  StmcbTVKernelptx = parallel.gpu.CUDAKernel('StmcbTVKernel.ptx', 'StmcbTVKernel.cu');
%% TODO  
%% TODO  for iii = 1:Niter
%% TODO    % perform a least square solve for the alpha and beta
%% TODO    [d_alpha, d_beta ] = feval(StmcbTVKernelptx,real(d_fid),imag(d_fid), d_alpha, d_beta, echospacing,imagingfreq,necho,nspecies,npixel);
%% TODO    % solve L1 subproblem for spatial regulariation of alpha beta 
%% TODO    w = l1subproblem([d_alpha, d_beta ],mu,lambda);
%% TODO    % update lagrange multipliers 
%% TODO   lambda =  lambda - 1/mu*(x-w);
%% TODO  end

% exit in order for NVVP to complete profiling
% exit

