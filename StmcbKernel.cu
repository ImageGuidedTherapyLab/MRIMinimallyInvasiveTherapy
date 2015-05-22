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
__device__
void WriteSolution(PetscComplex a[][MaxSize],int n,PetscComplex *x)
{
   int j,k;

   for (j=0;j<n;j++) {
      for (k=0;k<n+1;k++) {
         printf("%d %d %12.5e %12.5e ",k,j,a[k][j].real(),a[k][j].imag());
      }
      printf(" | %d  %12.5e %12.5e \n",j,x[j].real(),x[j].imag());
   }
   printf("\n");
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
         const double* d_RealDataArray,
         const double* d_ImagDataArray,
               double* d_Ppm,
               double* d_T2star,
               double* d_Amplitude,
               double* d_Phase,
         double const EchoSpacing,
         double const ImagingFreq,
         double const ThresholdSignal,
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
        
        if (idx == debugthread ) printf("idx=%d Necho=%d Nspecies=%d Npixel=%d\n",idx,Necho, Nspecies, Npixel);
              
        if (Necho > KernelMaxEcho)  // error check
          printf("Necho=%d > %d not supported \n",Necho,KernelMaxEcho);
        else if (Nspecies > MaxSpecies)  // error check
          printf("NSpecies=%d > %d not supported \n",Nspecies,MaxSpecies);
        else 
         {
           /* Copy Global data to register memory (ie Opencl Private) */
           double signalmagnitude = 0.0;
           for(iii = 0; iii< Necho; iii++)
            {
              datawork[iii] = PetscComplex(d_RealDataArray[idx * Necho+iii],d_ImagDataArray[idx * Necho+iii]);
              signalmagnitude  = signalmagnitude  + cusp::abs(datawork[iii]) ;
              sysinput[iii] = 0.0;
              if (idx == debugthread ) printf("idx=%d datawork[%d]=(%f,%f) sysinput[%d]=(%f,%f), mag = %f\n",idx,iii,datawork[iii].real(),datawork[iii].imag(),iii,sysinput[iii].real(),sysinput[iii].imag(),signalmagnitude  );
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
                    //if (idx == debugthread ) printf("idx=%d %d, %d,%d \n",idx,iii,kkk-iii-1,kkk);
                  }
               }

              // initialize solution
              for(iii = 0; iii< 2*Nspecies; iii++) slnVector[iii] = 0.0;

              /* solve prony linear system */
              if (idx == debugthread ) printf("Prony:\n");
              if (idx == debugthread ) WriteSolution(wrkMatrix,Nspecies,slnVector);
              GSolve(wrkMatrix,Nspecies,slnVector);
              if (idx == debugthread ) WriteSolution(wrkMatrix,Nspecies,slnVector);


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
                       printf("expected filter w[2]=(%12.5e,%12.5e) beta[0]=(%12.5e,%12.5e)\n",tmpcheck.real(),tmpcheck.imag(),slnVector[0].real(), slnVector[0].imag() );
                       }
                  signalfilter(Necho,datawork,Nspecies,slnVector);
                  signalfilter(Necho,sysinput,Nspecies,slnVector);
                  for(iii = 0; iii< Necho; iii++)
                     if (idx == debugverbose ) printf("idx=%d datawork[%d]=(%12.5e,%12.5e) sysinput[%d]=(%12.5e,%12.5e)\n",idx,iii,datawork[iii].real(),datawork[iii].imag(),iii,sysinput[iii].real(),sysinput[iii].imag());
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
                        //if (idx == debugthread ) printf("idx=%d %d, %d,%d \n",idx,iii,kkk-iii-1,kkk);
                      }
                   }

                  /* solve */
                  if (idx == debugverbose ) WriteSolution(wrkMatrix,2*Nspecies,slnVector);
                  GSolve(wrkMatrix,2*Nspecies,slnVector);
                  if (idx == debugthread ) WriteSolution(wrkMatrix,2*Nspecies,slnVector);
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
                  if (idx == displaythread ) printf("idx=%d alpha[0]=(%f,%f) alpha[1]=(%f,%f) beta[1]=(%f,%f) beta[2]=(%f,%f)\n",idx,
                                         slnVector[2].real(),slnVector[2].imag(),
                                         slnVector[3].real(),slnVector[3].imag(),
                                         slnVector[0].real(),slnVector[0].imag(),
                                         slnVector[1].real(),slnVector[1].imag()
                                       );
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

                  if (idx == debugthread ) WriteSolution(wrkMatrix,Nspecies,amplitude);
                  GSolve(wrkMatrix,2*Nspecies,amplitude);
                  if (idx == debugthread ) WriteSolution(wrkMatrix,Nspecies,amplitude);
                  if (idx == displaythread ) printf("idx=%d Lambda[0]=(%12.5e,%12.5e) Lambda[1]=(%12.5e,%12.5e) \n",idx,Lambda[0].real(),Lambda[0].imag(),Lambda[1].real(),Lambda[1].imag());
                  if (idx == displaythread ) printf("idx=%d Amp[0]=(%12.5e,%12.5e) Amp[1]=(%12.5e,%12.5e) \n",idx,amplitude[0].real(),amplitude[0].imag(),amplitude[1].real(),amplitude[1].imag());

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
                  
                  if (idx == displaythread ) printf("idx=%d ppm[%d]=%12.5e t2star[%d]=%12.5e amplitude[%d]=%12.5e phase[%d]=%12.5e  \n",idx,iii,d_Ppm[idx*Nspecies+iii], iii,d_T2star[   idx*Nspecies+iii], iii,d_Amplitude[idx*Nspecies+iii], iii,d_Phase[    idx*Nspecies+iii] );
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


