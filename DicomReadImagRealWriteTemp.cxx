/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: DicomSeriesReadSeriesWrite.cxx,v $
  Language:  C++
  Date:      $Date: 2005-11-19 16:31:50 $
  Version:   $Revision: 1.9 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

//  Software Guide : BeginLatex
//
//  This example illustrates how to read a DICOM series into a volume and then
//  save this volume into another DICOM series using the exact same header
//  information. It makes use of the GDCM library.
//
//  The main purpose of this example is to show how to properly propagate the
//  DICOM specific information along the pipeline to be able to correctly write
//  back the image using the information from the input DICOM files.
//
//  Please note that writing DICOM files is quite a delicate operation since we
//  are dealing with a significant amount of patient specific data. It is your
//  responsibility to verify that the DICOM headers generated from this code
//  are not introducing risks in the diagnosis or treatment of patients. It is
//  as well your responsibility to make sure that the privacy of the patient is
//  respected when you process data sets that contain personal information.
//  Privacy issues are regulated in the United States by the HIPAA
//  norms\footnote{The Health Insurance Portability and Accountability Act of
//  1996. \url{http://www.cms.hhs.gov/hipaa/}}. You would probably find similar
//  legislation in every country.
//
//  \index{HIPAA!Privacy}
//  \index{HIPAA!Dicom}
//  \index{Dicom!HIPPA}
//
//  When saving datasets in DICOM format it must be made clear whether this
//  datasets have been processed in any way, and if so, you should inform the
//  recipients of the data about the purpose and potential consequences of the
//  processing. This is fundamental if the datasets are intended to be used for
//  diagnosis, treatment or follow-up of patients. For example, the simple
//  reduction of a dataset form a 16-bits/pixel to a 8-bits/pixel
//  representation may make impossible to detect certain pathologies and as a
//  result will expose the patient to the risk or remaining untreated for a
//  long period of time while her/his pathology progresses.
//  
//  You are strongly encouraged to get familiar with the report on medical
//  errors ``To Err is Human'', produced by the U.S. Institute of
//  Medicine~\cite{ToErrIsHuman2001}. Raising awareness about the high
//  frequency of medical errors is a first step in reducing their occurrence.
//  
//  \index{Medical Errors}
//
//  Software Guide : EndLatex 

// Software Guide : BeginLatex
// 
// After all these warnings, let us now go back to the code and get familiar
// with the use of ITK and GDCM for writing DICOM Series. The first step that
// we must take is to include the header files of the relevant classes. We
// include the GDCM image IO class, the GDCM filenames generator, the series
// reader and writer.
//
// Software Guide : EndLatex

// Software Guide : BeginCodeSnippet
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include "itkExtractImageFilter.h"
#include "itkMedianImageFilter.h"
#include "itkAtan2ImageFilter.h"
#include "itkScalarToArrayCastImageFilter.h"
#include "itkGDCMImageIO.h"
#include "itkGDCMSeriesFileNames.h"
#include "itkImageSeriesReader.h"
#include "itkImageFileWriter.h"
#include "itkImageSliceConstIteratorWithIndex.h"
#include "itkImageSliceIteratorWithIndex.h"
// Software Guide : EndCodeSnippet

#include <vector>
#include <itksys/SystemTools.hxx>
#include <boost/lexical_cast.hpp>

// petsc includes
#include "petscda.h"
#include "petscsys.h"

# define OSSRealzeroright(o,v,p,d)  (o).width(v);  (o).precision(p); (o).fill('0'); (o) << std::right << (d)



// globals 
static char help[] = "Tests DAGetInterpolation for nonuniform DA coordinates.\n\n";
const unsigned int          Dimension = 3;
const double                       pi = 3.1415926535897931;
const PetscInt maxCSIecho=16,nvarplot=6;
// useful typedefs
typedef PetscScalar                                            InputPixelType;
typedef float                                                 OutputPixelType;
typedef itk::Image< InputPixelType, Dimension >               InputImageType;
typedef itk::Image< itk::Vector< PetscScalar , maxCSIecho >,
                                                 Dimension >    CSIImageType;
typedef itk::Image< OutputPixelType, Dimension >              OutputImageType;
typedef itk::Image< itk::Vector< OutputPixelType , nvarplot > , Dimension > 
                                                           VecOutputImageType;
typedef itk::ExtractImageFilter< InputImageType, InputImageType> 
                                                               ProcFilterType;
typedef itk::Atan2ImageFilter<InputImageType,InputImageType,InputImageType>
                                                              PhaseFilterType;
typedef itk::ScalarToArrayCastImageFilter< InputImageType, CSIImageType >
                                                                CSIFilterType;
typedef itk::ImageSeriesReader<      InputImageType >              ReaderType;
typedef itk::ImageFileWriter<       OutputImageType >              WriterType;
typedef itk::ImageFileWriter<    VecOutputImageType >           VecWriterType;
typedef itk::GDCMImageIO                                          ImageIOType;
typedef itk::ImageSliceConstIteratorWithIndex< CSIImageType > CSIIteratorType;
typedef itk::ImageSliceIteratorWithIndex< VecOutputImageType > VecOutIterType;
typedef itk::ImageSliceConstIteratorWithIndex< InputImageType >  
                                                            InputIteratorType;
typedef itk::ImageSliceIteratorWithIndex< InputImageType > BufferIteratorType;
typedef itk::ImageSliceIteratorWithIndex<OutputImageType > OutputIteratorType;

// subroutine to read the dicom data in real-time
template<class T>
InputImageType::Pointer GetImageData(T ITKPointer, const PetscInt rank,  std::vector<std::string>  &filenames )
{
  // Software Guide : BeginLatex
  // 
  // We can trigger the reading process by calling the \code{Update()}
  // method on the series reader. It is wise to put this invocation inside
  // a \code{try/catch} block since the process may eventually throw
  // exceptions.
  //
  // Software Guide : EndLatex
     
  bool DataNotRead = true;
  while(DataNotRead)
    {
     for (unsigned int i = 0; i < filenames.size(); i++)
      {
       if(access(filenames[i].c_str(),R_OK))
        {
         if(!rank) std::cout << filenames[i] << " NOT FOUND" << std::endl;
        }
       else 
        {
         if(!rank) std::cout << filenames[i].substr( 
                              filenames[i].rfind('/')+1, std::string::npos ) 
                   << " found...." << std::endl;
        }
      }
     try
       {
       ITKPointer->Update();
       DataNotRead = false;
       }
     catch (itk::ExceptionObject &excp)
       {
       std::cout << "Exception thrown" << excp << std::endl;
       sleep(1);
       }
    }
  std::cout << std::flush ;
  return ITKPointer->GetOutput();
}
// subroutine to write the image to disk
#undef __FUNCT__
#define __FUNCT__ "WriteImage"
void WriteImage(InputImageType::Pointer Image,std::ostringstream &filename,
                                        InputImageType::SizeType &filterRadius)
{
  PetscFunctionBegin;

  // setup writer
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( filename.str() );
  std::cout << "writing " << filename.str() << std::endl;

  //  Software Guide : BeginLatex
  //
  //  Using the image types, it is now possible to define the filter type
  //  and create the filter object.
  //
  //  \index{itk::MedianImageFilter!instantiation}
  //  \index{itk::MedianImageFilter!New()}
  //  \index{itk::MedianImageFilter!Pointer}
  // 
  //  Software Guide : EndLatex 

  // Software Guide : BeginCodeSnippet
  typedef itk::MedianImageFilter<
               InputImageType, OutputImageType >  OutputFilterType;

  OutputFilterType::Pointer filter = OutputFilterType::New();
  // Software Guide : EndCodeSnippet

  // Software Guide : BeginCodeSnippet
  
  filter->SetRadius( filterRadius );
  // Software Guide : EndCodeSnippet

  //  Software Guide : BeginLatex
  //
  //  The input to the filter can be taken from any other filter, for example
  //  a reader. The output can be passed down the pipeline to other filters,
  //  for example, a writer. An update call on any downstream filter will
  //  trigger the execution of the median filter.
  //
  //  \index{itk::MedianImageFilter!SetInput()}
  //  \index{itk::MedianImageFilter!GetOutput()}
  //
  //  Software Guide : EndLatex 


  // Software Guide : BeginCodeSnippet
  filter->SetInput( Image );
  writer->SetInput( filter->GetOutput() );
  // Software Guide : EndCodeSnippet
  try
    {
    writer->Update();
    }
  catch (itk::ExceptionObject &ex)
    {
    std::cout << ex; abort();
    }
  PetscFunctionReturnVoid();
}
// main data structure for thermal imaging
class RealTimeThermalImaging 
{
public:
  RealTimeThermalImaging(); //default constructor
  void DebugOn(); // enable debugging
  //setup data structures from command line arguements
  PetscErrorCode GetCommandLine(int ,char**);
  PetscErrorCode GetHeaderData(); //get dicom header info
  PetscErrorCode SetupDA();       //setup structured grid infrastructure
  PetscErrorCode FinalizeDA();    //close structured grid infrastructure
  PetscErrorCode GenerateCSITmap(); // generate CSI TMAP
  PetscErrorCode GeneratePRFTmap(); // generate PRF TMAP
  // generate phase images from real imaginary
  InputImageType::Pointer GetPhaseImage(PhaseFilterType::Pointer ,
                                        std::vector<std::string>  &,
                                        std::vector<std::string>  &); 
  // write an ini file containing dicom header info
  PetscErrorCode WriteIni();
  // get number of echo times
  int get_necho () const { return necho; }
private:
  // generate the expected file of a given image acquisition set
  PetscErrorCode RealTimeGenerateFileNames( const int , const int , 
                                         std::vector < std::vector < std::string > > &);
  // CSI Vector values images
  CSIImageType::Pointer GetCSIImage(CSIFilterType::Pointer, const int , 
                                                            const int );
  // write a paraview file from the DA data structures
  PetscErrorCode WriteDAImage(std::ostringstream &);
  int necho, ntime, median, noffset, nslice;    // defaults
  PetscScalar alpha,maxdiff;
  // dimensions for all images
  InputImageType::SpacingType sp;
  InputImageType::PointType orgn;
  // filtering parameters
  InputImageType::SizeType zeroFilterRadius, medianFilterRadius;
  // domain decomposition for DA data structures
  InputImageType::RegionType procRegion;
  // dicom readers
  ReaderType::Pointer  reader; 
  ImageIOType::Pointer gdcmIO;
  // image dimension info
  InputImageType::RegionType::SizeType size;
  // structured grid infrastructure
  DA                dac;
  Vec globalVec,imageVec;
  VecScatter     gather;

  // processor bounding box
  PetscInt ProcStart[3], ProcStartGhost[3];
  PetscInt ProcWidth[3], ProcWidthGhost[3];

  // to hold final image
  InputImageType::Pointer   net_Image;
  // command line arguments
  const char *ExamPath, *OutputDir; 
  int DirId;   
  // mpi rank
  PetscInt       rank;
  // error code
  PetscErrorCode ierr;
};
// constructor
#undef __FUNCT__
#define __FUNCT__ "RealTimeThermalImaging"
RealTimeThermalImaging::RealTimeThermalImaging()
{
  // defaults
  ntime=1; median=2; noffset=0; necho=1; nslice=5;   
  alpha = -0.0097; maxdiff=15.0; 

  //  Software Guide : BeginLatex
  //
  //  As a second step, we define the image type to be used in this example.
  //  This is done by explicitly selecting a pixel type and a dimension. Using
  //  the image type we can define the type of the series reader.
  //  
  // We construct one instance of the series reader object. Set the DICOM
  // image IO object to be use with it, and set the list of filenames to
  // read.
  // Software Guide : EndLatex
  
  // Software Guide : BeginCodeSnippet
  reader = ReaderType::New();
  gdcmIO = ImageIOType::New();
  reader->SetImageIO( gdcmIO );

  zeroFilterRadius[0] = 0; // radius along x
  zeroFilterRadius[1] = 0; // radius along y
  zeroFilterRadius[2] = 0; // radius along z
 
}
// get command line values
#undef __FUNCT__
#define __FUNCT__ "RealTimeThermalImaging::GetCommandLine"
PetscErrorCode RealTimeThermalImaging::GetCommandLine(int argc, char** argv)
{
  PetscFunctionBegin;

  // read input data
  ExamPath   =    argv[1] ;
  OutputDir  =    argv[3] ;
  try
  { 
     DirId                = boost::lexical_cast<int>(   argv[2]);
  }
  catch(const std::exception& e) //catch bad lexical cast
  {
     std::cout << "Error getting Command Line Values:\n"
               << "\nExamPath  = " << argv[1]
               << "\nDirId     = " << argv[2]
               << "\nOutputDir = " << argv[3] << std::endl;
     std::cout << e.what() << std::endl;
     return EXIT_FAILURE;
  }

  //error check
  if( argc < 4 )
    {
     std::cerr << "Usage: " << argv[0] <<
                  " ExamPath DirId OutputDir [options] "   << std::endl;
     std::cerr << "Available options: [default value]"     << std::endl ; 
     std::cerr << "-ntime     Number of Time Points: "     << ntime   << std::endl ;
     std::cerr << "-noffset   Time offset to begin with: " << noffset << std::endl ;
     std::cerr << "-necho     Number of echo Times: "      << necho   << std::endl ;
     std::cerr << "-nslice    Number of slices: "          << nslice  << std::endl ;
     std::cerr << "-alpha     Alpha (ppm/degC): "          << alpha   << std::endl ;
     std::cerr << "-maxdiff   filtering parameter: "       << maxdiff << std::endl ;
     std::cerr << "-median    median filtering parameter: "<< median  << std::endl ;
     std::cerr << "EchoTime and ImagingFrequency Optional\n" 
              << " but may be input in case of any errors\n" ;
     return EXIT_FAILURE;
    }

  /* Read options */
  ierr=PetscOptionsGetInt(PETSC_NULL,"-ntime",&ntime,PETSC_NULL);CHKERRQ(ierr);
  ierr=PetscOptionsGetInt(PETSC_NULL,"-necho",&necho,PETSC_NULL);CHKERRQ(ierr);
  ierr=PetscOptionsGetInt(PETSC_NULL,"-noffset",&noffset,PETSC_NULL);CHKERRQ(ierr);
  ierr=PetscOptionsGetInt(PETSC_NULL,"-nslice",&nslice,PETSC_NULL);CHKERRQ(ierr);
  ierr=PetscOptionsGetInt(PETSC_NULL,"-median",&median,PETSC_NULL);CHKERRQ(ierr);
  ierr=PetscOptionsGetScalar(PETSC_NULL,"-alpha",&alpha,PETSC_NULL);CHKERRQ(ierr);
  ierr=PetscOptionsGetScalar(PETSC_NULL,"-maxdiff",&maxdiff,PETSC_NULL);CHKERRQ(ierr);

  //  Software Guide : BeginLatex
  //
  //  The size of the neighborhood is defined along every dimension by
  //  passing a \code{SizeType} object with the corresponding values. The
  //  value on each dimension is used as the semi-size of a rectangular
  //  box. For example, in $2D$ a size of \(1,2\) will result in a $3 \times
  //  5$ neighborhood.
  //
  //  \index{itk::MedianImageFilter!Radius}
  //  \index{itk::MedianImageFilter!Neighborhood}
  //
  //  Software Guide : EndLatex 
  medianFilterRadius[0] = median; // radius along x
  medianFilterRadius[1] = median; // radius along y
  medianFilterRadius[2] = 0; // radius along z

  // make sure the output directory exist, using the cross
  // platform tools: itksys::SystemTools. In this case we select to create
  // the directory if it does not exist yet.
  //
  // \index{itksys!SystemTools}
  // \index{itksys!MakeDirectory}
  // \index{SystemTools}
  // \index{SystemTools!MakeDirectory}
  // \index{MakeDirectory!SystemTools}
  // \index{MakeDirectory!itksys}
  itksys::SystemTools::MakeDirectory( OutputDir );

  // get rank
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}
// turn on debuggin
#undef __FUNCT__
#define __FUNCT__ "RealTimeThermalImaging::DebugOn"
void RealTimeThermalImaging::DebugOn()
{
  PetscFunctionBegin;
  net_Image->DebugOn();
  //sp[0]= 1.00;
  //sp[1]= 1.00;
  //sp[2]= 1.00;
  //orgn[0] = 0.0;
  //orgn[1] = 0.0;
  //orgn[2] = 0.0;
  PetscFunctionReturnVoid();
}
//get dicom header data
#undef __FUNCT__
#define __FUNCT__ "RealTimeThermalImaging::GetHeaderData"
PetscErrorCode RealTimeThermalImaging::GetHeaderData()
{
  PetscFunctionBegin;
  // setup data structures to contain a list of file names for a single time 
  // instance
  // file names
  std::vector< std::vector<std::string> > filenames( 2 * necho , 
                             std::vector< std::string >::vector(nslice,"") );

  // generate first set of images filenames
  RealTimeGenerateFileNames(0,0,filenames);

  // only need to read in the first set to get the header info
  reader->SetFileNames(filenames[0] );
  InputImageType::Pointer tmpImage; // image pointer for header info
  tmpImage= GetImageData(reader, rank, filenames[0]);
  // scale image to meters
  sp = 0.001 * tmpImage->GetSpacing();
  if(!rank) std::cout << "Spacing = "
                      << sp[0] << ", " << sp[1] << ", " << sp[2] << std::endl;
  
  // scale image to meters
  const InputImageType::PointType& orgn_mm = tmpImage->GetOrigin();
  orgn[0] = 0.001 * orgn_mm[0];
  orgn[1] = 0.001 * orgn_mm[1];
  orgn[2] = 0.001 * orgn_mm[2];
  if(!rank) std::cout << "Origin = " << orgn[0] << ", " 
                                     << orgn[1] << ", " << orgn[2] << std::endl;

  // get size information
  size=tmpImage->GetRequestedRegion().GetSize();

  if(!rank) std::cout << "   ****Dimension Data****   " <<std::endl<<std::endl;
  PetscFunctionReturn(0);
}
//setup Petsc DA data structures
#undef __FUNCT__
#define __FUNCT__ "RealTimeThermalImaging::SetupDA"
PetscErrorCode RealTimeThermalImaging::SetupDA()
{
  PetscFunctionBegin;

  // domain decomposition for DA data structures
  InputImageType::RegionType::IndexType proc_start, proc_start_ghost;
  InputImageType::RegionType::SizeType  proc_width, proc_width_ghost;

  /* Create distributed array and get vectors:
       use -da_view to print out info about the DA */
  DAPeriodicType ptype = DA_NONPERIODIC;
  DAStencilType  stype = DA_STENCIL_BOX;
  ierr = DACreate3d(PETSC_COMM_WORLD,ptype,stype, size[0], size[1], size[2],
                    PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,nvarplot,1,
                    PETSC_NULL,PETSC_NULL,PETSC_NULL,&dac);CHKERRQ(ierr);
  ierr = DAGetCorners(dac,&ProcStart[0],&ProcStart[1],&ProcStart[2], 
                          &ProcWidth[0],&ProcWidth[1],&ProcWidth[2]); 
  CHKERRQ(ierr);
  proc_start[0]=ProcStart[0]; proc_width[0]=ProcWidth[0];
  proc_start[1]=ProcStart[1]; proc_width[1]=ProcWidth[1];
  proc_start[2]=ProcStart[2]; proc_width[2]=ProcWidth[2];
  
  ierr = DAGetGhostCorners(dac, 
          &ProcStartGhost[0],&ProcStartGhost[1],&ProcStartGhost[2], 
          &ProcWidthGhost[0],&ProcWidthGhost[1],&ProcWidthGhost[2] ); 
  CHKERRQ(ierr);
  proc_start_ghost[0]=ProcStartGhost[0]; proc_width_ghost[0]=ProcWidthGhost[0];
  proc_start_ghost[1]=ProcStartGhost[1]; proc_width_ghost[1]=ProcWidthGhost[1];
  proc_start_ghost[2]=ProcStartGhost[2]; proc_width_ghost[2]=ProcWidthGhost[2];
  
  ierr = DASetUniformCoordinates(dac,orgn[0],orgn[0]+sp[0]*size[0],
                                     orgn[1],orgn[1]+sp[1]*size[1],
                                     orgn[2],orgn[2]+sp[2]*size[2] );CHKERRQ(ierr);

  // create DA vectors
  DACreateGlobalVector(dac,&globalVec);
  VecScatterCreateToZero(globalVec,&gather,&imageVec);

  //  Software Guide : BeginLatex
  //  
  //  an \doxygen{ImageRegion} object is created and initialized with
  //  the start and size we just prepared using the slice information.
  //
  //  Software Guide : EndLatex 

  // Software Guide : BeginCodeSnippet
  procRegion.SetSize(  proc_width  );
  procRegion.SetIndex( proc_start );

  net_Image = InputImageType::New();
  net_Image->SetRegions( procRegion );
  net_Image->Allocate();
  net_Image->SetSpacing( sp );
  net_Image->SetOrigin( orgn);
  // Software Guide : EndCodeSnippet
  PetscFunctionReturn(0);
}
// Free Petsc Data structures
#undef __FUNCT__
#define __FUNCT__ "RealTimeThermalImaging::FinalizeDA"
PetscErrorCode RealTimeThermalImaging::FinalizeDA()
{
  PetscFunctionBegin;
  // destroy scatter context, vector, when no longer needed
  ierr = VecScatterDestroy(gather);CHKERRQ(ierr);
  ierr = VecDestroy(imageVec);CHKERRQ(ierr);
  //ierr = VecDestroy(globalVec);CHKERRQ(ierr); DA will destroy this
  ierr = DADestroy(dac);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
// subroutine to write the image to disk
#undef __FUNCT__
#define __FUNCT__ "RealTimeThermalImaging::WriteDAImage"
PetscErrorCode RealTimeThermalImaging::WriteDAImage(std::ostringstream &filename)
{
  PetscFunctionBegin;

  // allocate memory for the output image
  VecOutputImageType::RegionType region;
  region.SetSize(  size  );
  VecOutputImageType::Pointer outputImage = VecOutputImageType::New();
  outputImage->SetRegions( region );
  outputImage->Allocate();
  // set image dimensions
  outputImage->SetSpacing(sp);
  outputImage->SetOrigin(orgn);
    
  // initialize ITK iterators
  VecOutIterType    outputIt( outputImage, outputImage->GetRequestedRegion() );
  outputIt.SetFirstDirection(  0 ); outputIt.SetSecondDirection( 1 );
  outputIt.GoToBegin();

  // get pointer to petsc data structures
  PetscScalar   ****MapPixel;
  VecOutputImageType::PixelType   pixelValue;
  ierr = VecGetArray4d(imageVec,size[2],size[1],size[0],nvarplot,
                                        0,0,0,0,&MapPixel); CHKERRQ(ierr);

  /*
     loop through parallel data structures
  */
  for (PetscInt k=0; k<size[2]; k++) 
  {
    for (PetscInt j=0; j<size[1]; j++) 
    {
      for (PetscInt i=0; i<size[0]; i++) 
      {
        //outputIt.Set( &MapPixel[k][j][i][0] ) ; 
        for(PetscInt ivar = 0; ivar < nvarplot ; ivar++) 
                    pixelValue[ivar] = MapPixel[k][j][i][ivar];
        outputIt.Set( pixelValue ) ; 
        ++outputIt; // update iterators
      }
      outputIt.NextLine(); // get next line
    }
    outputIt.NextSlice(); // get next slice
  }
  ierr = VecRestoreArray4d(imageVec,size[2],size[1],size[0],nvarplot,
                                           0,0,0,0,&MapPixel); CHKERRQ(ierr);

  // setup writer
  VecWriterType::Pointer writer = VecWriterType::New();
  writer->SetFileName( filename.str() );
  std::cout << "writing " << filename.str() << std::endl;
  writer->SetInput( outputImage );
  try
    {
    writer->Update();
    }
  catch (itk::ExceptionObject &ex)
    {
    std::cout << ex; abort();
    }
  PetscFunctionReturn(0);
}
// generate proton resonance frequency thermal images
#undef __FUNCT__
#define __FUNCT__ "RealTimeThermalImaging::GeneratePRFTmap"
PetscErrorCode RealTimeThermalImaging::GeneratePRFTmap()
{
  PetscFunctionBegin;

  std::vector< std::vector<std::string> > filenames( 2 , 
                             std::vector< std::string >::vector(nslice,"") );
  RealTimeGenerateFileNames(0,0,filenames);

  // containers for phase data
  PhaseFilterType::Pointer  basephaseFilter=PhaseFilterType::New(),// n-1 phase data
                         currentphaseFilter=PhaseFilterType::New();// current phase data

  // containers for image pointers
  // get base image first and prepare temperature image
  InputImageType::Pointer baseImage, currentImage;

  {
    // get images
    baseImage = GetPhaseImage(basephaseFilter,filenames[0],filenames[1] );
    // output base phase image as a place holder
    std::ostringstream basephase_filename;
    // output file
    basephase_filename << OutputDir <<"/tmap."<< rank << ".";
    OSSRealzeroright(basephase_filename,4,0,0);
    basephase_filename << ".vtk" ;

    WriteImage(baseImage,basephase_filename,zeroFilterRadius);

  }

  bool EnhanceSNR=true;
  if(EnhanceSNR)
   {
     for( int iii = 1 ; iii <= nslice ; iii++)
      {
       ierr= PetscMatlabEnginePutArray(PETSC_MATLAB_ENGINE_WORLD,
                                 size[0],size[1],baseImage->GetBufferPointer(),
                                                  "PixelBuffer");CHKERRQ(ierr);
       ierr= PetscMatlabEngineEvaluate(PETSC_MATLAB_ENGINE_WORLD,
                             "BufferOut=EnhanceSNR(PixelBuffer)");
       CHKERRQ(ierr);
       ierr= PetscMatlabEngineGetArray(PETSC_MATLAB_ENGINE_WORLD,
                                 size[0],size[1],baseImage->GetBufferPointer(),
                                                  "BufferOut");CHKERRQ(ierr);
      }
   }

  // loop over time instances
  for( int iii = 1 ; iii <= ntime ; iii++)
   {
    // generate list of file names
    RealTimeGenerateFileNames(iii,0,filenames);

    // get images and header info
    double tmap_factor;
    currentImage = GetPhaseImage(currentphaseFilter,filenames[0], filenames[1] );
    // Software Guide : BeginLatex
    // 
    // We can trigger the reading process by calling the \code{Update()}
    // method on the series reader. It is wise to put this invocation inside
    // a \code{try/catch} block since the process may eventually throw
    // exceptions.
    //
    // Software Guide : EndLatex
    
    // scratch storage for header value
    std::string value; 
    // Get Echo Time 
    std::string echotimekey = "0018|0081";
    double echotime;
    try
    {  
       gdcmIO->GetValueFromTag(echotimekey, value) ;
       echotime=boost::lexical_cast<double>(value);
    }
    catch(const std::exception& e) //catch bad lexical cast
    {
       std::cout<<"Error getting echo time " << " (" << echotimekey << ")\n";
       std::cout<<"value returned is "<<value<< "\n";
       std::cout<<e.what() << std::endl;
       ierr = PetscOptionsGetScalar(PETSC_NULL,"-echotime",&echotime,PETSC_NULL);
       CHKERRQ(ierr);
       std::cout << "using value from command line "<< echotime << "\n";
    }
    std::string imagfreqkey = "0018|0084";
    double imagfreq;
    try
    {  
       gdcmIO->GetValueFromTag(imagfreqkey, value) ;
       // trailing space on string causing cast error
       value.erase(value.find(' ')) ; 
       imagfreq=boost::lexical_cast<double>(value);
    }
    catch(const std::exception& e) //catch bad lexical cast
    {
       std::cout << "Error getting Imaging Freq"<<" ("<<imagfreqkey << ")\n";
       std::cout << "value returned is "<<value << "\n";
       std::cout << e.what() << std::endl;
       ierr = PetscOptionsGetScalar(PETSC_NULL,"-imagfreq",&imagfreq,PETSC_NULL);
       CHKERRQ(ierr);
       std::cout << "using value from command line "<< imagfreq << "\n";
    }
    // misc info
    // echo data
    std::string studyIdkey  = "0020|0010" ;
    std::string seriesIdkey = "0020|0011" ;
    gdcmIO->GetValueFromTag(studyIdkey  , value) ;
    std::cout << "Study Id " << value  ;
    gdcmIO->GetValueFromTag(seriesIdkey , value) ;
    std::cout << " Series Number " << value   << "\n" ;
    std::cout << "alpha " << alpha  << " (ppm/degC) "
              << "maxdiff " << maxdiff << " (degC) "
              << "echo time " << echotime << " (ms) "
              << "imaging freq " << imagfreq << " (MHz) \n";
    // faster to multiply
    tmap_factor =  1.0/(2.0*pi*imagfreq*alpha*echotime*1.e-3) ;

   
    // Software Guide : BeginLatex
    // setup real, imaginary, base phase, and temperature map iterators
    // The const slice iterator walks the 3D input image, and the non-const
    // linear iterator walks the 2D output image. The iterators are initialized
    // to walk the same linear path through a slice.  Remember that the
    // \emph{second} direction of the slice iterator defines the direction that
    // linear iteration walks within a slice
    InputIteratorType    currIt(currentImage   ,currentImage->GetRequestedRegion());
    BufferIteratorType   baseIt(   baseImage   ,   baseImage->GetRequestedRegion());
    BufferIteratorType   net_It(   net_Image   ,   net_Image->GetRequestedRegion());
   
    currIt.SetFirstDirection(  0 );   currIt.SetSecondDirection( 1 );
    baseIt.SetFirstDirection(  0 );   baseIt.SetSecondDirection( 1 );
    net_It.SetFirstDirection(  0 );   net_It.SetSecondDirection( 1 );
   
    // Software Guide : EndCodeSnippet 
    // set iterators to the beginning
    currIt.GoToBegin();
    baseIt.GoToBegin();
    net_It.GoToBegin();
    // compute net phase difference
    while( !net_It.IsAtEnd() )
      {
      while ( !net_It.IsAtEndOfSlice() )
        {
        while ( !net_It.IsAtEndOfLine() )
          {
          net_It.Set(
             net_It.Get() +  ( currIt.Get() - baseIt.Get() ) * tmap_factor
                    );
          // update the base to contain the n-1 image for next time
          baseIt.Set( currIt.Get() );
          ++currIt;
          ++baseIt;
          ++net_It;
          }
          currIt.NextLine();
          baseIt.NextLine();
          net_It.NextLine();
        }
      // get next slice
      currIt.NextSlice(); 
      net_It.NextSlice();
      baseIt.NextSlice();
      }
   
    // output unfiltered tmap
    std::ostringstream unfiltered_filename;
    unfiltered_filename << OutputDir <<"/unfiltered."<< rank << ".";
    OSSRealzeroright(unfiltered_filename,4,0,iii);
    unfiltered_filename << ".vtk" ;

    WriteImage(net_Image , unfiltered_filename , zeroFilterRadius);

    // output filtered tmap
    std::ostringstream filtered_filename;
    filtered_filename << OutputDir <<"/filtered."<< rank << ".";
    OSSRealzeroright(filtered_filename,4,0,iii);
    filtered_filename << ".vtk" ;

    WriteImage(net_Image , filtered_filename , medianFilterRadius);

   } // end loop over time instances
  PetscFunctionReturn(0);
}
// generate CSI thermal images
#undef __FUNCT__
#define __FUNCT__ "RealTimeThermalImaging::GenerateCSITmap"
PetscErrorCode RealTimeThermalImaging::GenerateCSITmap() 
{
  PetscFunctionBegin;

  // file names
  std::vector< std::vector<std::string> > filenames( 2 * necho , 
                             std::vector< std::string >::vector(nslice,"") );

  // get header info
  std::vector< PetscScalar > echotime(necho,0.0), imagfreq(necho,0.0);
  for (int jjj = 0 ; jjj < necho ; jjj ++ )
    {
     // generate list of file names
     RealTimeGenerateFileNames(0,jjj,filenames);

     // only need to read in the first set to get the header info
     reader->SetFileNames(filenames[0] );
     InputImageType::Pointer tmpImage; // image pointer for header info
     tmpImage= GetImageData(reader, rank, filenames[0]);

     // scratch storage for header value
     std::string value; 
     // Get Echo Time 
     std::string echotimekey = "0018|0081";
     try
     {  
        gdcmIO->GetValueFromTag(echotimekey, value) ;
        echotime[jjj]=boost::lexical_cast<double>(value);
     }
     catch(const std::exception& e) //catch bad lexical cast
     {
        std::cout<<"Error getting echo time " << " (" << echotimekey << ")\n";
        std::cout<<"value returned is "<<value<< "\n";
        std::cout<<e.what() << std::endl; abort();
     }
     std::string imagfreqkey = "0018|0084";
     try
     {  
        gdcmIO->GetValueFromTag(imagfreqkey, value) ;
        // trailing space on string causing cast error
        value.erase(value.find(' ')) ; 
        imagfreq[jjj]=boost::lexical_cast<double>(value);
     }
     catch(const std::exception& e) //catch bad lexical cast
     {
        std::cout << "Error getting Imaging Freq"<<" ("<<imagfreqkey << ")\n";
        std::cout << "value returned is "<<value << "\n";
        std::cout << e.what() << std::endl;
        ierr = PetscOptionsGetScalar(PETSC_NULL,"-imagfreq",&imagfreq[jjj],PETSC_NULL);
        CHKERRQ(ierr);
        std::cout << "using value from command line "<< imagfreq[jjj] << "\n";
     }
     // echo data
     std::string studyIdkey  = "0020|0010" ;
     std::string seriesIdkey = "0020|0011" ;
     gdcmIO->GetValueFromTag(studyIdkey  , value) ;
     if(!rank) std::cout << "Study Id " << value  ;
     gdcmIO->GetValueFromTag(seriesIdkey , value) ;
     if(!rank) std::cout << " Series Number " << value 
                         << " echo time " << echotime[jjj]<< " (ms)"<<std::endl;
    }

  //Define gamma*B0 and the echo spacing
  PetscScalar      gB0=imagfreq[0], esp=echotime[1]-echotime[0];
  if(!rank) std::cout << "echo spacing " << esp   << " (ms) "
                      << "imaging freq " << gB0   << " (MHz) "
                      << "ntime  "       << ntime << std::endl;

  if(!rank) std::cout << "   ****Header Info****   " << std::endl << std::endl ;
  // containers for real and imaginary data
  CSIFilterType::Pointer   realImageFilter=CSIFilterType::New(),// real images
                           imagImageFilter=CSIFilterType::New();// imag images

  // pointers to real and imaginary vector images for CSI computation
  CSIImageType::Pointer realImages, imagImages;

  // loop over time instances
  for( int iii = 0 ; iii <= ntime ; iii++)
   {

    // get image data
    realImages = GetCSIImage(realImageFilter,iii,0); //real images
    imagImages = GetCSIImage(imagImageFilter,iii,1); //imag images

    // Software Guide : BeginLatex
    //
    // \index{Iterators!construction of} \index{Iterators!and image regions} The
    // necessary images and region definitions are now in place.  All that is
    // left to do is to create the iterators and perform the copy.  Note that
    // image iterators are not accessed via smart pointers so they are
    // light-weight objects that are instantiated on the stack.  Also notice how
    // the input and output iterators are defined over the \emph{same
    // corresponding region}.  Though the images are different sizes, they both
    // contain the same target subregion.
    //
    // Software Guide : EndLatex
    
    // Software Guide : BeginCodeSnippet

    // initialize ITK iterators
    CSIIteratorType realIt( realImages, realImages->GetRequestedRegion() );
    CSIIteratorType imagIt( imagImages, imagImages->GetRequestedRegion() );
    realIt.SetFirstDirection(  0 );   realIt.SetSecondDirection( 1 );
    imagIt.SetFirstDirection(  0 );   imagIt.SetSecondDirection( 1 );

    realIt.GoToBegin(); imagIt.GoToBegin();

    // get pointer to petsc data structures
    PetscScalar   ****MapPixel;
    DAGetGlobalVector(dac,&globalVec);
    DAVecGetArrayDOF(dac,globalVec,&MapPixel);
    //ierr = PetscObjectSetName((PetscObject)MapPixel,"MapPixel");CHKERRQ(ierr);

    PetscInt PixelCounter = 0 ; 
    /*
       loop through parallel data structures
    */
    for (PetscInt k=ProcStart[2]; k<ProcStart[2]+ProcWidth[2]; k++) 
    {
      for (PetscInt j=ProcStart[1]; j<ProcStart[1]+ProcWidth[1]; j++) 
      {
        for (PetscInt i=ProcStart[0]; i<ProcStart[0]+ProcWidth[0]; i++) 
        {
          //std::cout << "pixel # " << PixelCounter++ << std::endl;
          //for(PetscInt ivar = 0; ivar < nvarplot ; ivar++) 
          //          MapPixel[k][j][i][ivar] = realIt.Get()[ivar] ;
          ierr= PetscMatlabEnginePutArray(PETSC_MATLAB_ENGINE_WORLD,
                                          necho,1,&realIt.Get()[0],
                                                     "RealPixel");CHKERRQ(ierr);
          ierr= PetscMatlabEnginePutArray(PETSC_MATLAB_ENGINE_WORLD,
                                          necho,1,&imagIt.Get()[0],
                                                     "ImagPixel");CHKERRQ(ierr);
          ierr= PetscMatlabEnginePutArray(PETSC_MATLAB_ENGINE_WORLD,
                                          nvarplot,1,MapPixel[k][j][i],
                                                     "MapPixel");CHKERRQ(ierr);
          ierr= PetscMatlabEngineEvaluate(PETSC_MATLAB_ENGINE_WORLD,
                "MapPixel=MapCSI(RealPixel,ImagPixel,%18.16e,%18.16e)",gB0,esp);
          CHKERRQ(ierr);
          ierr= PetscMatlabEngineGetArray(PETSC_MATLAB_ENGINE_WORLD,
                                          nvarplot,1,MapPixel[k][j][i],
                                                     "MapPixel");CHKERRQ(ierr);
          ++realIt; ++imagIt; // update iterators
          CHKMEMQ; // check for memory corruption use -malloc_debug to enable
        }
        // get next line
        realIt.NextLine(); imagIt.NextLine();
      }
      // get next slice
      realIt.NextSlice(); imagIt.NextSlice();
    }

    /*
       Restore vector
    */
    DAVecRestoreArrayDOF(dac,globalVec,&MapPixel);
    DARestoreGlobalVector(dac,&globalVec);
   
    // gather the image buffer on processor zero
    VecScatterBegin(gather,globalVec,imageVec,INSERT_VALUES,SCATTER_FORWARD);
    VecScatterEnd(gather,globalVec,imageVec,INSERT_VALUES,SCATTER_FORWARD);
  
    // output various maps
    if(!rank) // write image from root process
     {
       std::ostringstream unfiltered_filename;
       unfiltered_filename << OutputDir <<"/unfiltered."<< rank << ".";
       OSSRealzeroright(unfiltered_filename,4,0,iii);
       unfiltered_filename << ".vtk" ;
       ierr= WriteDAImage(unfiltered_filename); CHKERRQ(ierr);
     }

   } // end loop over time instances

  PetscFunctionReturn(0);
}
// subroutine to generate dicom filenames
#undef __FUNCT__
#define __FUNCT__ "RealTimeThermalImaging::RealTimeGenerateFileNames"
PetscErrorCode 
RealTimeThermalImaging::RealTimeGenerateFileNames(const int istep,
                                                  const int iecho,
                         std::vector < std::vector < std::string > > &filenames)
{
   PetscFunctionBegin;
   // filenames[0][jjj] ==> real images 
   // filenames[1][jjj] ==> imaginary images
   for( int iii = 0 ; iii < 2 ; iii++)
     for( int jjj = 0 ; jjj < nslice ; jjj++)
      {
         std::ostringstream file_name;
         // hard code for now
         file_name << ExamPath << "/s" << DirId << "/i" 
                   <<   DirId     + 2*nslice*(necho*istep+iecho+noffset) + (2*jjj+1) + iii 
                   << ".MRDC."<<    2*nslice*(necho*istep+iecho+noffset) + (2*jjj+1) + iii ;
         filenames[iii][jjj] = file_name.str();
      }
  PetscFunctionReturn(0);
}
// return a phase image from the real/imaginary images
#undef __FUNCT__
#define __FUNCT__ "RealTimeThermalImaging::GetPhaseImage"
InputImageType::Pointer RealTimeThermalImaging::GetPhaseImage(
                                      PhaseFilterType::Pointer phaseFilter,
                                      std::vector<std::string>  &realFilenames,
                                      std::vector<std::string>  &imagFilenames )
{
  PetscFunctionBegin;
  //  The ExtractImageFilter type is instantiated using the input and
  //  output image types. A filter object is created with the New()
  //  method and assigned to a SmartPointer.
  ProcFilterType::Pointer realfilter=ProcFilterType::New();
  ProcFilterType::Pointer imagfilter=ProcFilterType::New();

  //  Software Guide : BeginLatex
  //  Then the region is passed to the filter using the
  //  SetExtractionRegion() method.
  //
  //  \index{itk::ExtractImageFilter!SetExtractionRegion()}
  //  Software Guide : EndLatex 

  // Software Guide : BeginCodeSnippet
  
  // get real images
  reader->SetFileNames(realFilenames );
  realfilter->SetInput( reader->GetOutput() );
  realfilter->SetExtractionRegion( procRegion );
  phaseFilter->SetInput1( GetImageData(realfilter, rank, realFilenames) );
  // get imaginary images
  reader->SetFileNames(imagFilenames );
  imagfilter->SetInput( reader->GetOutput() );
  imagfilter->SetExtractionRegion( procRegion );
  phaseFilter->SetInput2( GetImageData(imagfilter, rank, imagFilenames) );
  // get imaginary images
  // compute phase

  try
    {
    phaseFilter->Update();
    }
  catch (itk::ExceptionObject &excp)
    {
    std::cout << "Exception thrown" << excp << std::endl; abort();
    }

  // Software Guide : EndCodeSnippet
  PetscFunctionReturn(phaseFilter->GetOutput());
}
// return a phase image from the real/imaginary images
#undef __FUNCT__
#define __FUNCT__ "RealTimeThermalImaging::GetCSIImage"
CSIImageType::Pointer RealTimeThermalImaging::GetCSIImage(
                               CSIFilterType::Pointer CSIFilter,
                               const int istep, const int idData)
{
  PetscFunctionBegin;
  //  The ExtractImageFilter type is instantiated using the input and
  //  output image types. A filter object is created with the New()
  //  method and assigned to a SmartPointer.
  std::vector< ProcFilterType::Pointer > 
                                 procfilter( necho , ProcFilterType::New() );

  std::vector< std::vector<std::string> > filenames( 2 * necho , 
                             std::vector< std::string >::vector(nslice,"") );

  for (int jjj = 0 ; jjj < necho ; jjj ++ )
    {
     // generate list of file names
     RealTimeGenerateFileNames(istep,jjj,filenames);

     // get images at each echo time
     reader->SetFileNames( filenames[idData] );
     procfilter[jjj]->SetInput( reader->GetOutput() );
     procfilter[jjj]->SetExtractionRegion( procRegion );
     CSIFilter->PushBackInput( GetImageData( procfilter[jjj], rank, filenames[idData] ) );

    }
  // push back zeros for necho < maxCSIecho
  for (int jjj=necho;jjj<maxCSIecho;jjj++) CSIFilter->PushBackInput(net_Image);

  //  Software Guide : BeginLatex
  //  Then the region is passed to the filter using the
  //  SetExtractionRegion() method.
  //
  //  \index{itk::ExtractImageFilter!SetExtractionRegion()}
  //  Software Guide : EndLatex 

  // Software Guide : BeginCodeSnippet
  
  // get imaginary images
  // compute phase

  try
    {
    CSIFilter->Update();
    }
  catch (itk::ExceptionObject &excp)
    {
    std::cout << "Exception thrown" << excp << std::endl; abort();
    }

  // Software Guide : EndCodeSnippet
  PetscFunctionReturn(CSIFilter->GetOutput());
}
// write an Ini File containg dicom header info 
#undef __FUNCT__
#define __FUNCT__ "RealTimeThermalImaging::WriteIni"
PetscErrorCode RealTimeThermalImaging::WriteIni()
{
  PetscFunctionBegin;
  std::ofstream Inifile;
  std::ostringstream fileID ; // filename

  // output file
  fileID << OutputDir <<"/mrti.ini" ;
  Inifile.open(fileID.str().c_str(), std::ios::out);

  // dicom header parameters...
  Inifile <<"[mrti]" << std::endl;

  /* dimensions of MRTI data */

  Inifile <<"xpixel=" << size[0] << std::endl ;
  Inifile <<"ypixel=" << size[1] << std::endl ;
  Inifile <<"nslice=" << size[2] << std::endl ;

  /* physical space dimensions */

  // origin
  Inifile <<"x0=" << orgn[ 0] << std::endl ;
  Inifile <<"y0=" << orgn[ 1] << std::endl ;
  Inifile <<"z0=" << orgn[ 2] << std::endl ;

  // spacing
  Inifile <<"dx=" << sp[0] << std::endl ;
  Inifile <<"dy=" << sp[1] << std::endl ;
  Inifile <<"dz=" << sp[2] << std::endl ;

  // close file 
  Inifile.close();

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "main"
// main driver routine
int main( int argc, char* argv[] )
{
  PetscErrorCode ierr;

  /* Initialize Petsc */
  ierr = PetscInitialize(&argc,&argv,(char*)0,help);CHKERRQ(ierr); 
  PetscFunctionBegin;

  /* Initialize Main Class */
  RealTimeThermalImaging MRTI;

  //get commandline options
  ierr = MRTI.GetCommandLine(argc, argv);CHKERRQ(ierr); 

  //get dicom header data and write it out
  ierr = MRTI.GetHeaderData();CHKERRQ(ierr);
  ierr = MRTI.WriteIni();CHKERRQ(ierr);

  //debugging...
  PetscTruth  flg,Debug=PETSC_FALSE;
  ierr = PetscOptionsGetTruth(PETSC_NULL,"-debug",&Debug,&flg);
  if(Debug) MRTI.DebugOn();

  //setup DA arrays for parallel processing
  ierr = MRTI.SetupDA();CHKERRQ(ierr);

  //error check echo times
  if(MRTI.get_necho()>maxCSIecho)
   {
     std::cerr << "# of echo times is " << MRTI.get_necho() << std::endl ;
     std::cerr << "Max  echo times is " << maxCSIecho       << std::endl ;
     return EXIT_FAILURE;
   }

  // Software Guide : EndCodeSnippet
  switch(MRTI.get_necho())
   {
    case 1: // 1 echo --> multiplane tmap
      ierr = MRTI.GeneratePRFTmap();CHKERRQ(ierr);
      break;
    default: // 2 <= necho <= 16 --> CSI
      ierr = MRTI.GenerateCSITmap();CHKERRQ(ierr);
      break;
   }
     
  /* Free memory */
  ierr = MRTI.FinalizeDA();CHKERRQ(ierr);
  ierr = PetscFinalize();CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

