import vtk
# echo vtk version info
print "using vtk version", vtk.vtkVersion.GetVTKVersion()
import vtk.util.numpy_support as vtkNumPy 
import numpy
import os
import time
import dicom
import scipy.io as scipyio
# cuda environment
import pycuda.driver as cuda
import pycuda.autoinit
from pycuda.compiler import SourceModule

with open ("StmcbKernel.cu", "r") as cudafile:
    codetemplate=cudafile.read()
codetemplate = """
/*
 * Example of how to use the mxGPUArray API in a MEX file.  This example shows
 * how to write a MEX function that takes a gpuArray input and returns a
 * gpuArray output, e.g. B=mexFunction(A).
 *
 * Copyright 2012 The MathWorks, Inc.
 */
#include <cusp/complex.h>

// borrow petsc types
typedef double PetscScalar;
typedef cusp::complex<PetscScalar> PetscComplex;
typedef int PetscInt;

// FIXME global var used to pass array of data
int const MaxSpecies = 2;
int const MaxSize    = 2*MaxSpecies ;
int const tol        = 1.E-4;
int const upperbound = 100;


extern "C"{
/*
  filter operation
   \hat{Y} =  Y/D
*/
__device__ 
void  signalfilter(int Necho,PetscComplex *y,int Nspecies,PetscComplex beta[])
{
     for (int iii=0;iii<Necho;iii++) 
        for (int jjj=0;jjj<Nspecies;jjj++) 
           y[iii] = y[iii] - beta[jjj] *y[iii-(jjj+1)];
}
/*
 * Device code
   Solve a system of n equations in n unknowns using Gaussian Elimination
   Solve an equation in matrix form Ax = b
   The 2D array a is the matrix A with an additional column b.
   This is often written (A:b)

   TODO notice that the first index is the largest dimension
   A0,0    A1,0    A2,0    ....  An-1,0     b0
   A0,1    A1,1    A2,1    ....  An-1,1     b1
   A0,2    A1,2    A2,2    ....  An-1,2     b2
   :       :       :             :          :
   :       :       :             :          :
   A0,n-1  A1,n-1  A2,n-1  ....  An-1,n-1   bn-1
 */
__device__ 
void  GSolve(PetscComplex a[][MaxSize],int n,PetscComplex x[])
{
   int i,j,k; //,maxrow;
   PetscComplex tmp;

   for (i=0;i<n;i++) {

      ///* Find the row with the largest first value */
      //maxrow = i;
      //for (j=i+1;j<n;j++) {
      //   if (ABS(a[i][j]) > ABS(a[i][maxrow]))
      //      maxrow = j;
      //}

      ///* Swap the maxrow and ith row */
      //for (k=i;k<n+1;k++) {
      //   tmp = a[k][i];
      //   a[k][i] = a[k][maxrow];
      //   a[k][maxrow] = tmp;
      //}

      ///* Singular matrix? */
      //if (ABS(a[i][i]) < EPS)
      //   return(FALSE);

      /* Eliminate the ith element of the jth row */
      for (j=i+1;j<n;j++) {
         for (k=n;k>=i;k--) {
            a[k][j] -= a[k][i] * a[i][j] / a[i][i];
         }
      }
   }

   /* Do the back substitution */
   for (j=n-1;j>=0;j--) {
      tmp = 0;
      for (k=j+1;k<n;k++)
         tmp += a[k][j] * x[k];
      x[j] = (a[n][j] - tmp) / a[j][j];
   }

   return;
}

/*************QR Root Solve*************/
__device__
PetscComplex dotprod(PetscComplex *vec1, PetscComplex *vec2, int nDim)
{
	PetscComplex x = 0;
	PetscComplex tmp = 0;
	for (int i = 0; i < nDim; ++i)
	{
		tmp.real(vec1[i].real());
		tmp.imag(-vec1[i].imag());
		x += tmp * vec2[i];
	}
	return x;
}

__device__
PetscComplex* matmult(PetscComplex *mat1, PetscComplex *mat2, int nDim)
{
	PetscComplex *x = new PetscComplex[nDim * nDim];
	for (int i = 0; i < nDim * nDim; ++i)
		x[i] = 0;
	for (int k = 0; k < nDim; ++k)
		for (int j = 0; j < nDim; ++j)
			for (int i = 0; i < nDim; ++i)
				x[j + k * nDim] += mat1[j + i * nDim] * mat2[i + k * nDim];
	return x;
}

__device__
PetscComplex l2norm(PetscComplex *vec, int nDim)
{
	PetscComplex x = 0;
	for (int i = 0; i < nDim; ++i)
		x += vec[i].real() * vec[i].real() + vec[i].imag() * vec[i].imag();
	x = sqrt(x);
	return x;
}

__device__
void make_comp_mat(PetscComplex *polynomial, PetscComplex *companion, int nDim)
{
	for (int i = 0; i < nDim * nDim; ++i)
		companion[i] = 0;
	for (int i = 0; i < nDim; ++i)
		companion[i * nDim] = -polynomial[i + 1] / polynomial[0];
	for (int i = 0; i < nDim - 1; ++i)
		companion[i * nDim + i + 1] = 1;
}

__device__
void select_diag(PetscComplex *vector, PetscComplex *matrix, int nDim)
{
	for (int i = 0; i < nDim; ++i)
		vector[i] = matrix[i * nDim + i];
}

__device__
void modified_gram_schmidt(PetscComplex *a, PetscComplex *Q, PetscComplex *R, int nDim)
{
	PetscComplex *u = new PetscComplex[nDim];
	PetscComplex *v = new PetscComplex[nDim];
	PetscComplex prj = 0;
	PetscComplex l2 = 0;
	for (int i = 0; i < nDim * nDim; ++i)
		Q[i] = R[i] = 0;
	for (int i = 0; i < nDim; ++i)
		u[i] = v[i] = 0;

	for (int k = 0; k < nDim; ++k)
	{
		for (int i = 0; i < nDim; ++i)
			u[i] = a[i + k * nDim];
		for (int j = 0; j < k; ++j)
		{
			for (int i = 0; i < nDim; ++i)
				v[i] = Q[i + j * nDim];
			prj = dotprod(v, u, nDim);
			for (int i = 0; i < nDim; ++i)
				u[i] -= prj * v[i];
		}
		l2 = l2norm(u, nDim);
		for (int i = 0; i < nDim; ++i)
			Q[i + k * nDim] = u[i] / l2;
		for (int j = k; j < nDim; ++j)
		{
			for (int i = 0; i < nDim; ++i)
			{
				u[i] = a[i + j * nDim];
				v[i] = Q[i + k * nDim];
			}
			R[k + j * nDim] = dotprod(u, v, nDim);
		}
	}

	delete[] u;
	delete[] v;
}

__device__
void roots(
	PetscComplex *polynomial,
	PetscComplex *root,
	int nDim_in,
	double tolerance,
	int upperbound)
{
	int nDim = nDim_in - 1;
	PetscComplex *a = new PetscComplex[nDim * nDim];
	PetscComplex *Q = new PetscComplex[nDim * nDim];
	PetscComplex *R = new PetscComplex[nDim * nDim];
	int nTol = 0;
	for (int i = 0; i < nDim; ++i)
		root[i] = 0;

	make_comp_mat(polynomial, a, nDim);

	for (int k = 0; k < upperbound; ++k)
	{
		modified_gram_schmidt(a, Q, R, nDim);
		a = matmult(R, Q, nDim);
		nTol = 0;
		for (int j = 0; j < nDim; ++j)
			for (int i = 0; i < nDim; ++i) {
				if (i > j && sqrt(a[i + j * nDim].real() * a[i + j * nDim].real() + 
					a[i + j * nDim].imag() * a[i + j * nDim].imag()) > tolerance) 
					++nTol; }
		if (nTol == 0) break;
	}

	select_diag(root, a, nDim);

	delete[] a;
	delete[] Q;
	delete[] R;
}
/*************End QR Root Solve*************/

/*
 * Device code
 */
__global__ 
void StmcbKernel(
         const float* d_RealDataArray,
         const float* d_ImagDataArray,
               float* d_Ppm,
               float* d_T2star,
               float* d_Amplitude,
               float* d_Phase,
         float const EchoSpacing,
         float const ImagingFreq,
         float const ThresholdSignal,
         int const Necho,
         int const Nspecies,
         int const Npixel,
        const int debugthread   , 
        const int displaythread , 
        const int debugverbose  ) 
{
    /*
      grid stride loop design pattern, 1-d grid
      http://devblogs.nvidia.com/parallelforall/cuda-pro-tip-write-flexible-kernels-grid-stride-loops/
         - By using a loop, you can support any problem size even if it exceeds the largest grid size your CUDA device supports. Moreover, you can limit the number of blocks you use to tune performance.
    */
    for (int idx = blockIdx.x * blockDim.x + threadIdx.x; 
         idx < Npixel;
         idx += blockDim.x * gridDim.x) 
      {
        /* define temporary data structures in register memory */
        int iii,jjj,kkk;
        int const KernelMaxEcho = 16;
        // array for original data
        PetscComplex  dataworkfull[KernelMaxEcho+MaxSpecies];
        PetscComplex  sysinputfull[KernelMaxEcho+MaxSpecies];
        // initialize
        for(iii = 0; iii< KernelMaxEcho+MaxSpecies; iii++)
         {
           dataworkfull[iii] = 0.0;
           sysinputfull[iii] = 0.0;
         }
        // setup pointer to allow negative index during assembly
        PetscComplex  *datawork = &dataworkfull[MaxSpecies];
        PetscComplex  *sysinput = &sysinputfull[MaxSpecies];

        // storage for augmented matrix 
        PetscComplex  wrkMatrix[MaxSize+1][MaxSize];
        // storage for solution = [beta_1 ... beta_Nspecies alpha_0 ... alpha_{Nspecies-1}]
        PetscComplex  slnVector[MaxSize];
        
              
        if (Necho > KernelMaxEcho)  // error check
         {
          iii = 0;
         }
        else if (Nspecies > MaxSpecies)  // error check
         {
          iii = 0;
         }
        else 
         {
           /* Copy Global data to register memory (ie Opencl Private) */
           double signalmagnitude = 0.0;
           for(iii = 0; iii< Necho; iii++)
            {
              datawork[iii] = PetscComplex(d_RealDataArray[idx * Necho+iii],d_ImagDataArray[idx * Necho+iii]);
              signalmagnitude  = signalmagnitude  + cusp::abs(datawork[iii]) ;
              sysinput[iii] = 0.0;
            }
           // only compute for sufficient signal
           if(signalmagnitude > ThresholdSignal)
            {
              // system input is deltat function
              sysinput[0] = 1.0;

              /* initialize */
              for(iii = 0; iii< Nspecies+1; iii++)
                  for(jjj = 0; jjj< Nspecies; jjj++)
                    wrkMatrix[iii][jjj] = 0.0;

              /* Build Matrix and RHS for prony solve  */
              for(iii = 0; iii< Nspecies; iii++)
               {
                for(jjj = 0; jjj< Nspecies; jjj++)
                  {
                    for(kkk = Nspecies; kkk< Necho; kkk++)
                      {
                       // TODO Notice Gauss Solver is Row MAJOR
                       //wrkMatrix[iii][jjj] = wrkMatrix[iii][jjj] + cusp::conj(datawork[kkk-iii-1]) * datawork[kkk-jjj-1];
                       wrkMatrix[jjj][iii] = wrkMatrix[jjj][iii] + cusp::conj(datawork[kkk-iii-1]) * datawork[kkk-jjj-1];
                      }
                  }
                for(kkk = Nspecies; kkk< Necho; kkk++)
                  {
                    wrkMatrix[Nspecies][iii] = wrkMatrix[Nspecies][iii]  - cusp::conj(datawork[kkk-iii-1]) * datawork[kkk];
                  }
               }

              // initialize solution
              for(iii = 0; iii< 2*Nspecies; iii++) slnVector[iii] = 0.0;

              /* solve prony linear system */
              //if (idx == debugthread ) WriteSolution(wrkMatrix,Nspecies,slnVector);
              GSolve(wrkMatrix,Nspecies,slnVector);
              //if (idx == debugthread ) WriteSolution(wrkMatrix,Nspecies,slnVector);


              // buffer for stmcb iteration
              PetscComplex  dataworkstmcb[KernelMaxEcho+MaxSpecies];
              PetscComplex  sysinputstmcb[KernelMaxEcho+MaxSpecies];

              /* steiglitz iteration */
              for(int isteig = 0 ; isteig <5 ; isteig++)
                {
                   // initialize - restore pre-filter data
                   for(iii = 0; iii< KernelMaxEcho+MaxSpecies; iii++)
                    {
                      dataworkstmcb[iii] = dataworkfull[iii] ;
                      sysinputstmcb[iii] = sysinputfull[iii] ;
                    }
                   // initialize
                   for(iii = 0; iii< 2* Nspecies+1; iii++)
                       for(jjj = 0; jjj< 2* Nspecies; jjj++)
                         wrkMatrix[iii][jjj] = 0.0;
                   // FIXME - bad practice - difficult to follow
                   // setup pointer to allow negative index during assembly
                   datawork = &dataworkstmcb[MaxSpecies];
                   sysinput = &sysinputstmcb[MaxSpecies];

                     if (idx == debugthread )
                       {
                       PetscComplex tmpcheck=(datawork[2]-slnVector[0]*(datawork[1]-slnVector[0]*datawork[0])-slnVector[1]*datawork[0]);
                       }
                  signalfilter(Necho,datawork,Nspecies,slnVector);
                  signalfilter(Necho,sysinput,Nspecies,slnVector);
                  //for(iii = 0; iii< Necho; iii++)
                  /* Build Matrix and RHS for steiglitz  solve  */
                  for(iii = 0; iii< Nspecies; iii++)
                   {
                    // matrix
                    for(jjj = 0; jjj< Nspecies; jjj++)
                      {
                        for(kkk = 0; kkk< Necho; kkk++)
                          {
                           // TODO Notice Gauss Solver is Row MAJOR
                           wrkMatrix[jjj         ][iii         ] = wrkMatrix[jjj         ][iii         ] + cusp::conj(-datawork[kkk-iii-1]) *-datawork[kkk-jjj-1];
                           wrkMatrix[jjj+Nspecies][iii         ] = wrkMatrix[jjj+Nspecies][iii         ] + cusp::conj(-datawork[kkk-iii-1]) * sysinput[kkk-jjj  ];
                           wrkMatrix[jjj         ][iii+Nspecies] = wrkMatrix[jjj         ][iii+Nspecies] + cusp::conj( sysinput[kkk-iii  ]) *-datawork[kkk-jjj-1];
                           wrkMatrix[jjj+Nspecies][iii+Nspecies] = wrkMatrix[jjj+Nspecies][iii+Nspecies] + cusp::conj( sysinput[kkk-iii  ]) * sysinput[kkk-jjj  ];
                          }
                      }
                    // rhs
                    for(kkk = 0; kkk< Necho; kkk++)
                      {
                        wrkMatrix[2*Nspecies][iii         ] = wrkMatrix[2*Nspecies][iii         ] + cusp::conj(-datawork[kkk-iii-1]) * datawork[kkk];
                        // offest by Nspecies to for the signal input
                        wrkMatrix[2*Nspecies][iii+Nspecies] = wrkMatrix[2*Nspecies][iii+Nspecies] + cusp::conj( sysinput[kkk-iii  ]) * datawork[kkk];
                      }
                   }

                  /* solve */
                  //if (idx == debugverbose ) WriteSolution(wrkMatrix,2*Nspecies,slnVector);
                  GSolve(wrkMatrix,2*Nspecies,slnVector);
                  //if (idx == debugthread ) WriteSolution(wrkMatrix,2*Nspecies,slnVector);
                }

              // analytic 1 peak 
              PetscComplex  Lambda[MaxSpecies];
              PetscComplex  amplitude[MaxSpecies];
              // initialize amplitude
              for(iii = 0; iii< Nspecies; iii++) amplitude[iii] = 0.0;

              if ( Nspecies ==  1 )  
                {
                  /* compute roots */
                  Lambda[0] = - slnVector[0] ;
                  /* compute amplitude from residue */
                  amplitude[0] = (slnVector[1])/(slnVector[0]*Lambda[0]);
                }
              // analytic 2 peak 
              else if ( Nspecies ==  2 )  
                {
                  /* compute roots */
                  Lambda[0] = 0.5 * ( slnVector[0] + sqrt(slnVector[0]*slnVector[0] + 4.0 * slnVector[1]) ); 
                  Lambda[1] = 0.5 * ( slnVector[0] - sqrt(slnVector[0]*slnVector[0] + 4.0 * slnVector[1]) ); 
               
                  /* compute amplitude from initial conditions*/
                  wrkMatrix[0][0] = 1;
                  wrkMatrix[1][0] = 1;
                  wrkMatrix[0][1] = Lambda[0]+slnVector[0];
                  wrkMatrix[1][1] = Lambda[1]+slnVector[0];
                  wrkMatrix[2][0] = slnVector[2];
                  wrkMatrix[2][1] = slnVector[3];

                  //if (idx == debugthread ) WriteSolution(wrkMatrix,Nspecies,amplitude);
                  GSolve(wrkMatrix,2*Nspecies,amplitude);
                  //if (idx == debugthread ) WriteSolution(wrkMatrix,Nspecies,amplitude);

                  /* compute amplitude from residue */
                  /* FIXME: extend to multiple species */
                  for(iii = 0; iii< Nspecies; iii++)
                    amplitude[iii] = (slnVector[2]+slnVector[3]*Lambda[iii])/(2.0*slnVector[0]+slnVector[1]*Lambda[iii])/Lambda[iii];

                }
              // TODO Npeaks > 2 needed
              else if (Nspecies > MaxSpecies)  // error check
               	roots(slnVector, Lambda, Nspecies+1, tol, upperbound);

              // compute amplitudes, frequency, and t2star
              for(iii = 0; iii< Nspecies; iii++)
                {
                  PetscComplex  logroot = cusp::log(Lambda[iii]);
                  d_Ppm[      idx*Nspecies+iii] = logroot.imag()/2.0/M_PI/(EchoSpacing*ImagingFreq) * 1.e-3;
                  d_T2star[   idx*Nspecies+iii] = -EchoSpacing/logroot.real() ;

                  d_Amplitude[idx*Nspecies+iii] = cusp::abs(amplitude[iii]) ;
                  d_Phase[    idx*Nspecies+iii] = cusp::arg(amplitude[iii]) ;
                  
                }

            }
           else
            {
              // return default values for no signal
              for(iii = 0; iii< Nspecies; iii++)
                {
                  d_Ppm[      idx*Nspecies+iii] = 0.0;
                  d_T2star[   idx*Nspecies+iii] = 0.0;
                  d_Amplitude[idx*Nspecies+iii] = 0.0;
                  d_Phase[    idx*Nspecies+iii] = 0.0;
                }
            }
         } 
      } // end grid stride loop design pattern, 1-d grid
}

}
"""
mod = SourceModule(codetemplate,include_dirs = ['/opt/apps/cuda/4.2/cuda/include/'],no_extern_c=True)
gpukernel = mod.get_function("StmcbKernel")



# FIXME - hack to reduce file usage
SetFalseToReduceFileSystemUsage  = False

class RealTimeDicomFileRead:
  """ Base Class for realtime image header parsing...  """
  def __init__(self,rootDirectory,ExpectedFileSize,DefaultNstep,DefaultOffset,RemoteServer,RemoteRsync):
    print " base class constructor called \n\n" 
    # dictionary key template = time echo slice type
    self.keyTemplate   = "%04d_%03d_%03d_%02d" 
    self.FileSize = ExpectedFileSize
    self.DicomDataDictionary = {}
    self.FilesReadIn = set()
    self.NumTimeStep = DefaultNstep
    self.MinRawDataNumber = 100000000
    self.TimeOffset = DefaultOffset
    self.SignalThreshold = 5.e2
    # Debug flags
    self.Debug = 0
    # assume local directory
    self.dataDirectory = rootDirectory
    # remote server
    self.RemoteServer = RemoteServer 
    self.SocketFile = "/tmp/%r@%h:%p "
    self.CheckConnectionCMD = None
    self.KillConnectionCMD  = None
    self.RemoteDirectory    = None
    self.RsyncCMD           = None
    if (self.RemoteServer != None):
       self.RemoteDirectory    = rootDirectory
       self.dataDirectory      = "RawData/%s" % rootDirectory.split('/').pop()
       os.system( "mkdir -p %s" % self.dataDirectory)
       self.CheckConnectionCMD = "ssh -O check -S %s %s" % (self.SocketFile,self.RemoteServer) 
       self.CreateSocketCMD    = "ssh -MNf     -S %s %s" % (self.SocketFile,self.RemoteServer) 
       self.KillConnectionCMD  = "ssh -O exit  -S %s %s" % (self.SocketFile,self.RemoteServer) 
       # setup rsync command
       self.RsyncCMD = 'rsync -e "ssh -S %s"  ' %  (self.SocketFile)
       if(RemoteRsync != None):
         self.RsyncCMD = self.RsyncCMD + ' --rsync-path=%s ' % (RemoteRsync)
       self.RsyncCMD = self.RsyncCMD + ' -avz %s:%s/ %s/ ' %  (self.RemoteServer,self.RemoteDirectory,self.dataDirectory)
       # setup connection
       if(os.system(self.CheckConnectionCMD  ) ):
          print "\n\n   Creating Persisent Connection  %s!!!! \n\n " % self.RemoteServer
          print self.CreateSocketCMD  
          os.system(self.CreateSocketCMD )
       else:
          print "\n\n   Found Persisent Connection on %s!!!! \n\n " % self.RemoteServer
  def __del__(self):
    if( self.KillConnectionCMD != None):
      import os
      print self.KillConnectionCMD 
      os.system(self.KillConnectionCMD)

  def GetHeaderInfo(self):
    """ get initial header info"""
    # get initial file data
    dcmHeaderFileName =  self.QueryDictionary(0,1,0,2)

    headerData = dicom.read_file( "%s/%s" % (self.dataDirectory,dcmHeaderFileName) )

    # get number of slices
    self.nslice = headerData[0x0021,0x104f].value
    # set number of species
    self.NSpecies    = 2
    # get number of echos
    self.NumberEcho  = headerData[0x0019,0x107e].value
    self.Npixel         =  headerData.Rows*headerData.Rows*self.nslice
    self.FullSizeRaw    =  headerData.Rows*headerData.Rows*self.nslice*self.NumberEcho 
    self.RawDimensions  = [headerData.Rows,headerData.Rows,self.nslice,self.NumberEcho]
    self.FullSizeMap    =  headerData.Rows*headerData.Rows*self.nslice*self.NSpecies 
    self.MapDimensions  = [headerData.Rows,headerData.Rows,self.nslice,self.NSpecies ]
    spacing_mm = headerData.PixelSpacing
    spacing_mm.append(headerData.SpacingBetweenSlices) 
    origin_mm  = headerData.ImagePositionPatient
    #convert to meter
    self.spacing = [ 0.001 * float(dXi) for dXi in spacing_mm  ]
    self.origin  = [ 0.001 * float(Xi) for Xi in origin_mm   ]
    print self.RawDimensions, self.NSpecies, self.spacing, self.origin
    
    # FIXME should be negative but phase messed up somewhere
    self.alpha = -0.0097   
    # FIXME need to read in multiple echo times to process MFGRE should 
    # FIXME still be fine for 1st echo
    echoTime = float(headerData.EchoTime)
    self.imagFreq    = float(headerData.ImagingFrequency)
    # FIXME read from header
    self.EchoSpacing = 1.e-3
    # temperature map factor
    self.tmap_factor = 1.0 / (2.0 * numpy.pi * self.imagFreq * self.alpha * echoTime * 1.e-3)

    # should be equal and imaginary data
    expectedNtime = int(headerData.ImagesinAcquisition)/self.nslice/self.NumberEcho/2
    if( self.NumTimeStep != expectedNtime ):
       print headerData.ImagesinAcquisition,self.nslice,self.NumberEcho
       #raise RuntimeError("expecting %d total time points" % expectedNtime )
    # end GetHeaderInfo(self):
    return

  # return image data from raw file names
  def QueryDictionary(self,timeInstance,echo_id,slice_id,imagetype):
    "infinite loop until file is available and of proper size "
    localFileKey = self.keyTemplate % (timeInstance,echo_id,slice_id,imagetype) 
    while ( True ):
      try: 
        # get current directory list 
        directoryList = os.listdir(self.dataDirectory)
        # check if this has already been read in
        filename = self.DicomDataDictionary[localFileKey][0]
        return filename 
      except OSError: 
        print "waiting for directory %s ... " % ( self.dataDirectory) 
        time.sleep(1)
      except KeyError: 
        print "waiting ... min raw %d offset %d time %04d Echo %03d slice %03d type %02d" % (self.MinRawDataNumber,self.TimeOffset  ,  timeInstance,echo_id,slice_id,imagetype) 
        # rsync remote directory
        if(self.RemoteServer != None):
           print self.RsyncCMD 
           os.system( self.RsyncCMD )
        else:
           time.sleep(1)
        files = set( filter(lambda x:os.path.isfile("%s/%s" % (self.dataDirectory,x) ) ,directoryList) )
        ### filestmp = filter(lambda x:os.path.isfile("%s/%s" % (self.dataDirectory,x) ) ,directoryList) 
        ### files = set (filter(lambda x:int(x.split(".").pop()) < 50 ,filestmp ) )

        # we will only read files that have not been read
        FilesNotYetRead = files - self.FilesReadIn
        for filename in FilesNotYetRead :
          FullPathToFile = "%s/%s" %(self.dataDirectory,filename) 
          if ( os.path.getsize( FullPathToFile ) > self.FileSize ):
            print "found", FullPathToFile
            dcmimage = dicom.read_file( FullPathToFile )
            #deltat = dcmimage[0x0019,0x105a].value/dcmimage[0x0019,0x10f2].value*1.e-6
            deltat = 5.0
            rawdataNumber = dcmimage[0x0019,0x10a2].value
            numEchoes = dcmimage[0x0019,0x107e].value
            #check if default ntime not set
            if (self.NumTimeStep == None):
              try:
                self.NumTimeStep = dcmimage.NumberofTemporalPositions
              except AttributeError:
                raise RuntimeError("NumberofTemporalPositions Not found try setting directly\n\t\t--nstep=...")
            #compute timeIntID
            rawdataNumber = dcmimage[0x0019,0x10a2].value
            if( rawdataNumber < self.MinRawDataNumber ):
              self.MinRawDataNumber = rawdataNumber
              #need to start over
              self.FilesReadIn.clear()
            else:
              #sliceIntID = int(abs(round((dcmimage.SliceLocation - dcmimage[0x0019,0x1019].value)/ dcmimage.SpacingBetweenSlices)))
              # slice location grouped by temporal position
              #  ie slice one goes first, then slice two , three, etc
              sliceIntID = (rawdataNumber - self.MinRawDataNumber )/ dcmimage.NumberofTemporalPositions

              #if(sliceIntID < 0  or sliceIntID >= self.nslice):
              #  print "SliceLocation", dcmimage.SliceLocation , "spacing between slices",dcmimage.SpacingBetweenSlices, "first scan location " , dcmimage[0x0019,0x1019].value, "nslice", self.nslice
              if(sliceIntID < 0  ):
                print "SliceLocation", dcmimage.SliceLocation , "spacing between slices",dcmimage.SpacingBetweenSlices, "first scan location " , dcmimage[0x0019,0x1019].value
                raise RuntimeError("slice integer %d ? " % sliceIntID )
              if ( numEchoes == 1 ) : 
                 # for 1 echo assume CPD
                 numberSlice = dcmimage[0x0021,0x104f].value
                 timeIntID = int(dcmimage.InstanceNumber - 1)/int(numberSlice*2)
              elif( numEchoes > 1 ) :
                 # for multiple echo assume MFGRE
                 tmptimeID = dcmimage.InstanceNumber - 1 - self.NumTimeStep * numEchoes * sliceIntID * 2
                 timeIntID = tmptimeID /numEchoes / 2
                 timeIntID = rawdataNumber - self.TimeOffset
                 timeIntID = rawdataNumber - sliceIntID*dcmimage.NumberofTemporalPositions  - self.MinRawDataNumber 
              else :
                 raise RuntimeError("unknown sequence ")
              #error check
              if( timeIntID < 0  or timeIntID >= self.NumTimeStep ):
                 print 'timeIntID', timeIntID ,"numEchoes ", numEchoes 
                 print "InstanceNumber ", dcmimage.InstanceNumber, "NumberofTemporalPositions", self.NumTimeStep, "sliceIntID" ,sliceIntID 
                 #print "TriggerTime", dcmimage.TriggerTime, "deltat", deltat, "number slice ", dcmimage[0x0021,0x104f].value
                 print "time error: %d not \\notin [0,%d) " % (timeIntID,self.NumTimeStep) 
              datatype = int(dcmimage[0x0043,0x102f].value)
              #error check
              if ( datatype == 2 or datatype == 3 ) : 
                keyID = self.keyTemplate % ( timeIntID, int(dcmimage.EchoNumbers), 
                                   sliceIntID, datatype ) 
              else :
                 raise RuntimeError("\n\n\t unknown datatype %d : expecting real and imaginary data" % datatype)
              #error check key
              if ( keyID in self.DicomDataDictionary) : 
                 "\n\n\t overwriting keyID %s ..." % keyID 
              self.DicomDataDictionary[keyID]=(filename,dcmimage.EchoTime)
              print "deltat", deltat," min raw", self.MinRawDataNumber,"raw data", rawdataNumber, "key", keyID
              # not all headers have this ? 
              try:
                print "trigger ",dcmimage.TriggerTime,"temporal ID ",dcmimage.TemporalPositionIdentifier
              except :
                pass
              # ensure we do not read this in again
              self.FilesReadIn.add( filename )
          else: 
            print "filesize too small", os.path.getsize( FullPathToFile )
        print "read in ", len(self.FilesReadIn), "files"

  # return image data from raw file names
  def GetRawDICOMData(self,idtime,outDirectoryID):
    
    # loop until files are ready to be read in
    realImageFilenames = []
    imagImageFilenames = []
    # FIXME: index fun, slice start from 0, echo start from 1
    for idslice in range(self.nslice):
      for idecho in range(1,self.NumberEcho+1):
        realImageFilenames.append( self.QueryDictionary( idtime,idecho,idslice,2 ) )  
        imagImageFilenames.append( self.QueryDictionary( idtime,idecho,idslice,3 ) )  
  
    #create local vars
    rootdir = self.dataDirectory
    RawDim  = self.RawDimensions
    real_array=numpy.zeros(self.FullSizeRaw,dtype=numpy.float32)
    imag_array=numpy.zeros(self.FullSizeRaw,dtype=numpy.float32) 

    vtkAppendReal = vtk.vtkImageAppendComponents()
    vtkAppendImag = vtk.vtkImageAppendComponents()

    for idEchoLoc,(fileNameReal,fileNameImag) in enumerate(zip(realImageFilenames,imagImageFilenames)):
      # FIXME: index nightmare
      # FIXME: will be wrong for different ordering
      # arrange such that echo varies fast then x, then y, then z
      # example of how slicing behaves
      #>>> x = range(100)
      #>>> x
      #[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99]
      #>>> x[0:100:10]
      #[0, 10, 20, 30, 40, 50, 60, 70, 80, 90]
      #>>> x[1:100:10]
      #[1, 11, 21, 31, 41, 51, 61, 71, 81, 91]
      #>>> x[2:100:10]
      #[2, 12, 22, 32, 42, 52, 62, 72, 82, 92]
      #>>> x[3:100:10]
      #[3, 13, 23, 33, 43, 53, 63, 73, 83, 93]
      idEcho  = idEchoLoc % RawDim[3]
      idSlice = idEchoLoc / RawDim[3]
      beginIndex = RawDim[0]*RawDim[1]*RawDim[3]* idSlice   +idEcho
      finalIndex = RawDim[0]*RawDim[1]*RawDim[3]*(idSlice+1)
      stepIndex  = RawDim[3]
  
      ## realds = dicom.read_file( "%s/%s"%(rootdir,fileNameReal) )
      ## imagds = dicom.read_file( "%s/%s"%(rootdir,fileNameImag) )
      ## realsliceID = int( round((float(realds.SliceLocation) - float(realds[0x0019,0x1019].value))/ realds.SliceThickness))
      ## imagsliceID = int( round((float(imagds.SliceLocation) - float(imagds[0x0019,0x1019].value))/ imagds.SliceThickness))
  
      ## print "%03d echo %03d slice %03d slice [%d:%d:%d] %03d %s %03d %s "% (idEchoLoc,idEcho,idSlice,beginIndex,finalIndex,stepIndex,realsliceID,fileNameReal,imagsliceID,fileNameImag )
  
      vtkRealDcmReader = vtk.vtkDICOMImageReader()
      vtkRealDcmReader.SetFileName("%s/%s"%(rootdir,fileNameReal) )
      vtkRealDcmReader.Update()
      vtkRealData = vtk.vtkImageCast()
      vtkRealData.SetOutputScalarTypeToFloat()
      vtkRealData.SetInput( vtkRealDcmReader.GetOutput() )
      vtkRealData.Update( )
      real_image = vtkRealData.GetOutput().GetPointData() 
      real_array[ beginIndex: finalIndex : stepIndex ] = vtkNumPy.vtk_to_numpy(real_image.GetArray(0)) 
  
      vtkImagDcmReader = vtk.vtkDICOMImageReader()
      vtkImagDcmReader.SetFileName("%s/%s"%(rootdir,fileNameImag) )
      vtkImagDcmReader.Update()
      vtkImagData = vtk.vtkImageCast()
      vtkImagData.SetOutputScalarTypeToFloat()
      vtkImagData.SetInput( vtkImagDcmReader.GetOutput() )
      vtkImagData.Update( )
      imag_image = vtkImagData.GetOutput().GetPointData() 
      imag_array[ beginIndex: finalIndex : stepIndex ] = vtkNumPy.vtk_to_numpy(imag_image.GetArray(0)) 
  
      vtkAppendReal.SetInput( idEchoLoc ,vtkRealDcmReader.GetOutput() )
      vtkAppendImag.SetInput( idEchoLoc ,vtkImagDcmReader.GetOutput() )
      vtkAppendReal.Update( )
      vtkAppendImag.Update( )
    if (SetFalseToReduceFileSystemUsage):
      vtkRealDcmWriter = vtk.vtkDataSetWriter()
      vtkRealDcmWriter.SetFileName("Processed/%s/realrawdata.%04d.vtk" % (outDirectoryID,idtime) )
      vtkRealDcmWriter.SetInput(vtkAppendReal.GetOutput())
      vtkRealDcmWriter.Update()
  
    if (SetFalseToReduceFileSystemUsage):
      vtkImagDcmWriter = vtk.vtkDataSetWriter()
      vtkImagDcmWriter.SetFileName("Processed/%s/imagrawdata.%04d.vtk" % (outDirectoryID,idtime) )
      vtkImagDcmWriter.SetInput(vtkAppendImag.GetOutput())
      vtkImagDcmWriter.Update()
  
    # write numpy to disk in matlab
    echoTimes = []
    for idecho in range(1,self.NumberEcho+1):
       localKey = self.keyTemplate % ( idtime,idecho,0,2 )
       echoTimes.append(self.DicomDataDictionary[localKey][1])
    if (SetFalseToReduceFileSystemUsage):
      scipyio.savemat("Processed/%s/rawdata.%04d.mat"%(outDirectoryID,idtime), {'dimensions':RawDim,'echoTimes':echoTimes,'real':real_array,'imag':imag_array})
  
    # end GetRawDICOMData
    return (real_array,imag_array)

  # write a numpy data to disk in vtk format
  def ConvertNumpyVTKImage(self,NumpyImageData):
    # Create initial image
    MapDim = self.MapDimensions
    # imports raw data and stores it.
    dataImporter = vtk.vtkImageImport()
    # array is converted to a string of chars and imported.
    data_string = NumpyImageData.tostring()
    dataImporter.CopyImportVoidPointer(data_string, len(data_string))
    # The type of the newly imported data is set to unsigned char (uint8)
    dataImporter.SetDataScalarTypeToFloat()
    # Because the data that is imported only contains an intensity value (it isnt RGB-coded or someting similar), the importer
    # must be told this is the case.
    dataImporter.SetNumberOfScalarComponents(MapDim[3])
    # The following two functions describe how the data is stored and the dimensions of the array it is stored in. For this
    # simple case, all axes are of length 75 and begins with the first element. For other data, this is probably not the case.
    # I have to admit however, that I honestly dont know the difference between SetDataExtent() and SetWholeExtent() although
    # VTK complains if not both are used.
    dataImporter.SetDataExtent( 0, MapDim[0]-1, 0, MapDim[1]-1, 0, MapDim[2]-1)
    dataImporter.SetWholeExtent(0, MapDim[0]-1, 0, MapDim[1]-1, 0, MapDim[2]-1)
    dataImporter.SetDataSpacing( self.spacing )
    dataImporter.SetDataOrigin(  self.origin )
    dataImporter.Update()
    return dataImporter.GetOutput()
  
# setup command line parser to control execution
from optparse import OptionParser
parser = OptionParser()
parser.add_option("--datadir", 
                  action="store", dest="datadir", default=None,
                  help="[REQUIRED] full path to data directory", metavar="DIR")
parser.add_option("--remoteserver", 
                  action="store", dest="remoteserver", default=None,
                  help="[OPTIONAL] transfer files from remoteserver:datadir", metavar="IP")
parser.add_option("--remotersync", 
                  action="store", dest="remotersync", default=None,
                  help="[OPTIONAL] transfer files using remoteserver:remotesrsync", metavar="PATH")
parser.add_option("--speciesID", 
                  action="store", dest="speciesID", type="int", default=0,
                  help="[OPTIONAL] species # to display ", metavar="INT")
parser.add_option("--sliceID", 
                  action="store", dest="sliceID", type="int", default=None,
                  help="[OPTIONAL] slice to display", metavar="INT")
parser.add_option("--nstep", 
                  action="store", dest="nstep", type="int", default=None,
                  help="[OPTIONAL] # of expected time steps ", metavar="INT")
parser.add_option("--offset", 
                  action="store", dest="offset", type="int", default=0,
                  help="[OPTIONAL] time offset of raw data number", metavar="INT")
parser.add_option("--baseline", 
                  action="store", dest="baseline", type="float", default=0.0,
                  help="[OPTIONAL] initial temperature ", metavar="FLOAT")
parser.add_option("-q", "--quiet",
                  action="store_false", dest="verbose", default=True,
                  help="don't print status messages to stdout")
(options, args) = parser.parse_args()
if (options.datadir != None):
  
  # instantiate helper class
  fileHelper = RealTimeDicomFileRead( options.datadir, 256*256,options.nstep ,options.offset
                                                     ,options.remoteserver ,options.remotersync)
  outputDirID = filter(len,options.datadir.split("/")).pop()
  #os.system( "mkdir -p Processed/%s" % outputDirID )
  try:
    os.mkdir( "Processed/%s" % outputDirID )
  except OSError:
    print "Processed/%s exists..." % outputDirID 

  # Get Header data
  fileHelper.GetHeaderInfo( )
  print fileHelper.tmap_factor

  # display the center slice by default
  # FIXME should we display more than one slice ? 
  # FIXME should we display the max of all slice ? 
  DisplaySlice = fileHelper.nslice/2
  if ( options.sliceID!= None ) : 
     DisplaySlice = options.sliceID
  
  print "# time steps %d, display slice %d " % ( fileHelper.NumTimeStep,DisplaySlice )

  deltat = 6.0
  pvd=open("Processed/%s/temperature.pvd" % outputDirID ,"w")
  pvd.write('<?xml version="1.0"?>\n')
  pvd.write('<VTKFile type="Collection" version="0.1" byte_order="LittleEndian" compressor="vtkZLibDataCompressor">\n')
  pvd.write('  <Collection>\n')
  for idtime in range(fileHelper.NumTimeStep):
       pvd.write('   <DataSet timestep="%f" part="0" file="%s.%04d.vti"/>\n' % (idtime*deltat,"temperature",idtime) )
  pvd.write('  </Collection>\n')
  pvd.write('</VTKFile>\n')
  
  # create initial image as 1d array
  absTemp = numpy.zeros(fileHelper.FullSizeMap,
                         dtype=numpy.float32) + options.baseline
  vtkTempImage = fileHelper.ConvertNumpyVTKImage(absTemp)
  vtkTempWriter = vtk.vtkXMLImageDataWriter()
  vtkTempWriter.SetFileName( "Processed/%s/temperature.%04d.vti" % (outputDirID,0))
  vtkTempWriter.SetInput( vtkTempImage )
  vtkTempWriter.Update()
  
  # create a rendering window and renderer
  ren = vtk.vtkRenderer()
  renWin = vtk.vtkRenderWindow()
  renWin.AddRenderer(ren)
   
  # create a renderwindowinteractor
  iren = vtk.vtkRenderWindowInteractor()
  iren.SetRenderWindow(renWin)
  
  # try image viewer for multiple windows
  imageViewer = vtk.vtkImageViewer2()

  # setup storage
  ppm_array      =numpy.zeros(fileHelper.FullSizeMap,dtype=numpy.float32)
  t2star_array   =numpy.zeros(fileHelper.FullSizeMap,dtype=numpy.float32)
  amplitude_array=numpy.zeros(fileHelper.FullSizeMap,dtype=numpy.float32)
  phase_array    =numpy.zeros(fileHelper.FullSizeMap,dtype=numpy.float32)

  # get initial data
  vtkPreviousImage = fileHelper.GetRawDICOMData( 0, outputDirID )

  # configure kernel
  threadsPerBlock = 256;
  blocksPerGrid = (fileHelper.Npixel  + threadsPerBlock - 1) / threadsPerBlock;

  # run kernel
  gpukernel(cuda.In(vtkPreviousImage[0] ),cuda.In(vtkPreviousImage[1] ),
            cuda.Out(ppm_array       ),cuda.Out(t2star_array   ),
            cuda.Out(amplitude_array ),cuda.Out(phase_array    ),
            numpy.array(fileHelper.EchoSpacing    ,dtype=numpy.float32),
            numpy.array(fileHelper.imagFreq       ,dtype=numpy.float32),
            numpy.array(fileHelper.SignalThreshold,dtype=numpy.float32),
            numpy.array(fileHelper.NumberEcho     ,dtype=numpy.int32),
            numpy.array(fileHelper.NSpecies       ,dtype=numpy.int32),
            numpy.array(fileHelper.Npixel         ,dtype=numpy.int32),
            numpy.array(fileHelper.Debug          ,dtype=numpy.int32),
            numpy.array(fileHelper.Debug          ,dtype=numpy.int32),
            numpy.array(fileHelper.Debug          ,dtype=numpy.int32),
                  block=(threadsPerBlock,1, 1), grid=(blocksPerGrid,1) )

  # do not finish until all files processed
  for idfile in range(fileHelper.NumTimeStep): 
    try:
      # get current data set
      vtkCurrent_Image = fileHelper.GetRawDICOMData( idfile, outputDirID )

      # save previous ppm info
      prevppm_array    = numpy.copy(ppm_array)

      # run kernel
      gpukernel(cuda.In(vtkCurrent_Image[0] ),cuda.In(vtkCurrent_Image[1] ),
                cuda.Out(ppm_array       ),cuda.Out(t2star_array   ),
                cuda.Out(amplitude_array ),cuda.Out(phase_array    ),
                numpy.array(fileHelper.EchoSpacing    ,dtype=numpy.float32),
                numpy.array(fileHelper.imagFreq       ,dtype=numpy.float32),
                numpy.array(fileHelper.SignalThreshold,dtype=numpy.float32),
                numpy.array(fileHelper.NumberEcho     ,dtype=numpy.int32),
                numpy.array(fileHelper.NSpecies       ,dtype=numpy.int32),
                numpy.array(fileHelper.Npixel         ,dtype=numpy.int32),
                numpy.array(fileHelper.Debug          ,dtype=numpy.int32),
                numpy.array(fileHelper.Debug          ,dtype=numpy.int32),
                numpy.array(fileHelper.Debug          ,dtype=numpy.int32),
                      block=(threadsPerBlock,1, 1), grid=(blocksPerGrid,1) )
                            
      deltaTemp = (ppm_array - prevppm_array )/fileHelper.alpha 
      absTemp  = absTemp + deltaTemp 
  
      # write numpy to disk in vtk format
      print "writing timeID %d " % (idfile)
      vtkTempImage = fileHelper.ConvertNumpyVTKImage(absTemp)
      vtkTempWriter = vtk.vtkXMLImageDataWriter()
      vtkTempWriter.SetFileName( "Processed/%s/temperature.%04d.vti" % (outputDirID,idfile))
      vtkTempWriter.SetInput( vtkTempImage )
      vtkTempWriter.Update()
  
      # write numpy to disk in matlab
      if (SetFalseToReduceFileSystemUsage):
        scipyio.savemat("Processed/%s/temperature.%05d.mat"%(outputDirID,idfile), {'temp':absTemp})
  
      # update for next time step
      vtkPreviousImage = vtkCurrent_Image 
  
      # color table
      # http://www.vtk.org/doc/release/5.8/html/c2_vtk_e_3.html#c2_vtk_e_vtkLookupTable
      # http://vtk.org/gitweb?p=VTK.git;a=blob;f=Examples/ImageProcessing/Python/ImageSlicing.py
      hueLut = vtk.vtkLookupTable()
      hueLut.SetNumberOfColors (256)
      #FIXME: adjust here to change color  range
      hueLut.SetRange (-5.0, 20.0)  
      #hueLut.SetSaturationRange (0.0, 1.0)
      #hueLut.SetValueRange (0.0, 1.0)
      hueLut.SetHueRange (0.667, 0.0)
      hueLut.SetRampToLinear ()
      hueLut.Build()
  
      # colorbar
      # http://www.vtk.org/doc/release/5.8/html/c2_vtk_e_3.html#c2_vtk_e_vtkLookupTable
      scalarBar = vtk.vtkScalarBarActor()
      scalarBar.SetTitle("Temperature")
      scalarBar.SetNumberOfLabels(4)
      scalarBar.SetLookupTable(hueLut)
  
      # image viewer
      imageViewer.SetInput(vtkTempImage)
      imageViewer.SetSize(512,512)
      imageViewer.SetSlice(DisplaySlice)
      imageViewer.SetPosition(512,0)
      imageViewer.Render()

      # extract VOI to display
      extractVOI = vtk.vtkExtractVOI()
      extractVOI.SetInput(vtkTempImage) 
      extractVOI.SetVOI([0,fileHelper.MapDimensions[0],0,fileHelper.MapDimensions[1],DisplaySlice,DisplaySlice]) 
      extractVOI.Update()
      
      # mapper
      #mapper = vtk.vtkDataSetMapper()
      mapper = vtk.vtkImageMapToColors()
      mapper.SetInput(extractVOI.GetOutput())
      # set echo to display
      mapper.SetActiveComponent( options.speciesID )
      mapper.SetLookupTable(hueLut)
  
      # actor
      actor = vtk.vtkImageActor()
      actor.SetInput(mapper.GetOutput())
       
      # assign actor to the renderer
      ren.AddActor(actor)
      ren.AddActor2D(scalarBar)
       
      # uncomment to enable user interface interactor
      #iren.Initialize()
      renWin.SetSize(512,512)
      renWin.Render()
      #iren.Start()

      # save as a movie for animation
      windowToImage = vtk.vtkWindowToImageFilter() 
      windowToImage.SetInput(renWin)
      windowToImage.Update()
      if (SetFalseToReduceFileSystemUsage):
        jpgWriter     = vtk.vtkJPEGWriter() 
        jpgWriter.SetFileName( "Processed/%s/temperature.%04d.jpg" % (outputDirID,idfile))
        #jpgWriter.SetInput(extractVOI.GetOutput())
        jpgWriter.SetInput(windowToImage.GetOutput())
        jpgWriter.Write()

    except KeyboardInterrupt:
      #reset reference phase
      print "reseting base phase image at time ", idfile
      time.sleep(1)
      vtkPreviousImage = fileHelper.GetRawDICOMData( idfile, outputDirID )
      absTemp = numpy.zeros(fileHelper.FullSizeMap,
                            dtype=numpy.float32) + options.baseline
  # save gif animations
  print "saving animations... "
  os.system( "convert -delay 30 -resize 80%% -loop 0 Processed/%s/temperature.*.jpg Processed/%s/temperature.gif  " % (outputDirID,outputDirID) ) 
else:
  parser.print_help()
  print options
