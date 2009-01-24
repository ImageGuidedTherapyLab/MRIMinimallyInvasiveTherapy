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
// useful typedefs
typedef double                                                 InputPixelType;
typedef float                                                 OutputPixelType;
typedef itk::Image<  InputPixelType, Dimension >               InputImageType;
typedef itk::Image< OutputPixelType, Dimension >              OutputImageType;
typedef itk::ExtractImageFilter< InputImageType, InputImageType>   FilterType;
typedef itk::ImageSeriesReader<      InputImageType >              ReaderType;
typedef itk::ImageFileWriter<       OutputImageType >              WriterType;
typedef itk::GDCMImageIO                                          ImageIOType;
typedef itk::ImageSliceConstIteratorWithIndex< InputImageType >  
                                                            InputIteratorType;
typedef itk::ImageSliceIteratorWithIndex< InputImageType > BufferIteratorType;
typedef itk::ImageSliceIteratorWithIndex<OutputImageType > OutputIteratorType;

// subroutine to generate dicom filenames
void RealTimeGenerateFileNames(const char * ExamPath, const int DirId, 
                               const int istep, const int nslice,
                               const int necho, const int noffset,
                         std::vector < std::vector < std::string > > &filenames)
{
   // filenames[0][jjj] ==> real images 
   // filenames[1][jjj] ==> imaginary images
   for( int iii = 0 ; iii < 2 ; iii++)
     for( int jjj = 0 ; jjj < nslice ; jjj++)
      {
         std::ostringstream file_name;
         // hard code for now
         file_name << ExamPath << "/s" << DirId << "/i" 
                   <<   DirId     + 2*nslice*(istep+noffset) + (2*jjj+1) + iii 
                   << ".MRDC."<<    2*nslice*(istep+noffset) + (2*jjj+1) + iii ;
         filenames[iii][jjj] = file_name.str();
      }
   return;
}
// subroutine to read the dicom data in real-time

template<class T>
InputImageType::ConstPointer GetImageData(T ITKPointer, std::vector<std::string>  &filenames )
{
  bool DataNotRead = true;
  while(DataNotRead)
    {
     std::cout << std::endl ; 
     for (unsigned int i = 0; i < filenames.size(); i++)
      {
       if(access(filenames[i].c_str(),R_OK))
         std::cout << filenames[i] << " NOT FOUND\n";
       else 
         std::cout << filenames[i] << " found....\n";
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
  return ITKPointer->GetOutput();
}
// write an Ini File containg dicom header info 
void WriteIni( const char * OutputDir,
               InputImageType::PointType            &origin,
               InputImageType::SpacingType          &spacing,
               InputImageType::RegionType::SizeType &size)
{
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
  Inifile <<"x0=" << origin[ 0] << std::endl ;
  Inifile <<"y0=" << origin[ 1] << std::endl ;
  Inifile <<"z0=" << origin[ 2] << std::endl ;

  // spacing
  Inifile <<"dx=" << spacing[0] << std::endl ;
  Inifile <<"dy=" << spacing[1] << std::endl ;
  Inifile <<"dz=" << spacing[2] << std::endl ;

  // close file 
  Inifile.close();

  return;
}
// subroutine to write the image to disk
void WriteImage(InputImageType::Pointer Image, std::ostringstream &filename,
                                               InputImageType::SizeType &filterRadius)
{
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
  return;
}

#undef __FUNCT__
#define __FUNCT__ "main"
// main driver routine
int main( int argc, char* argv[] )
{
  PetscErrorCode ierr;
  PetscInt       rank;

  /* Initialize Petsc */
  ierr = PetscInitialize(&argc,&argv,(char*)0,help);CHKERRQ(ierr); 
  PetscFunctionBegin;
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank); CHKERRQ(ierr);

  // read input data
  const char * ExamPath   =                             argv[1] ;
  const char * OutputDir  =                             argv[3] ;
  int DirId;   
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

  /* Read options */
  int ntime=1, median=2, noffset=0, necho=1, nslice=5;    // defaults
  PetscScalar alpha = -0.0097,maxdiff=15.0; 
  PetscTruth  flg,Debug=PETSC_FALSE;//debugging...
  ierr = PetscOptionsGetInt(PETSC_NULL,"-ntime",&ntime,PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(PETSC_NULL,"-necho",&necho,PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(PETSC_NULL,"-noffset",&noffset,PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(PETSC_NULL,"-nslice",&nslice,PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(PETSC_NULL,"-median",&median,PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetScalar(PETSC_NULL,"-alpha",&alpha,PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetScalar(PETSC_NULL,"-maxdiff",&maxdiff,PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetTruth(PETSC_NULL,"-debug",&Debug,&flg);
 
  //error check
  if( argc < 4 )
    {
     std::cerr << "Usage: " << argv[0] <<
                  " ExamPath DirId OutputDir [options] " << std::endl;
     std::cerr << "Available options "                   <<              std::endl ; 
     std::cerr << "-ntime     Number of Time Points "    << ntime     << std::endl ;
     std::cerr << "-noffset   Time offset to begin with" << noffset   << std::endl ;
     std::cerr << "-necho     Number of echo Times  "    << necho     << std::endl ;
     std::cerr << "-nslice    Number of slices "         << nslice    << std::endl ;
     std::cerr << "-alpha     Alpha (ppm/degC) "         << alpha     << std::endl ;
     std::cerr << "-maxdiff   filtering parameter  "     << maxdiff   << std::endl ;
     std::cerr << "-median    median filtering parameter"<< median    << std::endl ;
     std::cerr << "EchoTime and ImagingFrequency Optional\n" 
              << " but may be input in case of any errors\n" ;
     return EXIT_FAILURE;
    }


  // setup data structures to contain a list of file names for a single time 
  // instance
  std::vector< std::vector<std::string> > 
       filenames( 2 * necho , std::vector< std::string >::vector(nslice,"") );

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

  //  Software Guide : BeginLatex
  //
  //  As a second step, we define the image type to be used in this example.
  //  This is done by explicitly selecting a pixel type and a dimension. Using
  //  the image type we can define the type of the series reader.
  //  
  //  Software Guide : EndLatex 

  // get base image first and prepare temperature image
  InputImageType::Pointer   baseImage = InputImageType::New();
  InputImageType::Pointer   net_Image = InputImageType::New();

  // image pointers
  std::vector< InputImageType::ConstPointer > realImage; realImage.resize(necho);
  std::vector< InputImageType::ConstPointer > imagImage; imagImage.resize(necho);
  // Software Guide : EndCodeSnippet

  // Software Guide : BeginLatex
  // We construct one instance of the series reader object. Set the DICOM
  // image IO object to be use with it, and set the list of filenames to
  // read.
  // Software Guide : EndLatex
  
  // Software Guide : BeginCodeSnippet
  ReaderType::Pointer  reader = ReaderType::New();
  ImageIOType::Pointer gdcmIO = ImageIOType::New();
  reader->SetImageIO( gdcmIO );

  if(Debug)
    {
      baseImage->DebugOn();
      net_Image->DebugOn();
    }

  // dimensions for all images
  InputImageType::SpacingType sp;
  InputImageType::PointType orgn;
 
  // generate first set of images filenames
  RealTimeGenerateFileNames(ExamPath,DirId,0,nslice,necho,noffset,filenames);

  // only need to read in the first set to get the header info
  reader->SetFileNames(filenames[0] );
  realImage[0]= GetImageData(reader, filenames[0]);

  // scale image to meters
  sp = 0.001 * realImage[0]->GetSpacing();
  sp[0]= 1.00;
  sp[1]= 1.00;
  sp[2]= 1.00;
  std::cout << "Spacing = ";
  std::cout << sp[0] << ", " << sp[1] << ", " << sp[2] << std::endl;
  
  // scale image to meters
  const InputImageType::PointType& orgn_mm = realImage[0]->GetOrigin();
  orgn[0] = 0.001 * orgn_mm[0];
  orgn[1] = 0.001 * orgn_mm[1];
  orgn[2] = 0.001 * orgn_mm[2];
  orgn[0] = 0.0;
  orgn[1] = 0.0;
  orgn[2] = 0.0;
  std::cout << "Origin = ";
  std::cout << orgn[0] << ", " << orgn[1] << ", " << orgn[2] << std::endl;

  // write out header info
  InputImageType::RegionType::SizeType
                              size=realImage[0]->GetRequestedRegion().GetSize();

  WriteIni(OutputDir,orgn,sp,size);

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
  InputImageType::SizeType zeroFilterRadius, medianFilterRadius;
  zeroFilterRadius[0] = 0; // radius along x
  zeroFilterRadius[1] = 0; // radius along y
  zeroFilterRadius[2] = 0; // radius along z
  medianFilterRadius[0] = median; // radius along x
  medianFilterRadius[1] = median; // radius along y
  medianFilterRadius[2] = 0; // radius along z


  // processor bounding box
  PetscInt ProcStart[3], ProcStartGhost[3];
  PetscInt ProcWidth[3], ProcWidthGhost[3];
  InputImageType::RegionType::IndexType proc_start, proc_start_ghost;
  InputImageType::RegionType::SizeType  proc_width, proc_width_ghost;

  /* Create distributed array and get vectors:
       use -da_view to print out info about the DA */
  DA             dac;
  DAPeriodicType ptype = DA_NONPERIODIC;
  DAStencilType  stype = DA_STENCIL_BOX;
  ierr = DACreate3d(PETSC_COMM_WORLD,ptype,stype, size[0], size[1], size[2],
                    PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,1,1,
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

  //  Software Guide : BeginLatex
  //  
  //  an \doxygen{ImageRegion} object is created and initialized with
  //  the start and size we just prepared using the slice information.
  //
  //  Software Guide : EndLatex 

  // Software Guide : BeginCodeSnippet
  InputImageType::RegionType procRegion;
  procRegion.SetSize(  proc_width  );
  procRegion.SetIndex( proc_start );
  // Software Guide : EndCodeSnippet

  //  The ExtractImageFilter type is instantiated using the input and
  //  output image types. A filter object is created with the New()
  //  method and assigned to a SmartPointer.
  std::vector< FilterType::Pointer > realfilter(necho,FilterType::New());
  std::vector< FilterType::Pointer > imagfilter(necho,FilterType::New());

  //  Software Guide : BeginLatex
  //  Then the region is passed to the filter using the
  //  SetExtractionRegion() method.
  //
  //  \index{itk::ExtractImageFilter!SetExtractionRegion()}
  //  Software Guide : EndLatex 

  // Software Guide : BeginCodeSnippet
  for (int jjj = 0 ; jjj < necho ; jjj ++ )
    {
      realfilter[jjj]->SetExtractionRegion( procRegion );
      imagfilter[jjj]->SetExtractionRegion( procRegion );
      realfilter[jjj]->SetInput( reader->GetOutput() );
      imagfilter[jjj]->SetInput( reader->GetOutput() );
    }
  // Software Guide : EndCodeSnippet

  // allocate base phase image and output image buffer from first set
  baseImage->SetRegions( procRegion );
  net_Image->SetRegions( procRegion );
  baseImage->Allocate();
  net_Image->Allocate();

  {
    // get images
    for (int jjj = 0 ; jjj < necho ; jjj ++ )
      {
       // get real images
       reader->SetFileNames(filenames[2*jjj  ] );
       realImage[jjj] = GetImageData(realfilter[jjj], filenames[2*jjj  ] );
       // get imaginary images
       reader->SetFileNames(filenames[2*jjj+1] );
       imagImage[jjj] = GetImageData(imagfilter[jjj], filenames[2*jjj+1] );
      }
    // setup real, imaginary, base phase, and temperature map iterators
    InputIteratorType    realIt(   realImage[0],  realImage[0]->GetRequestedRegion());
    InputIteratorType    imagIt(   imagImage[0],  imagImage[0]->GetRequestedRegion());
    BufferIteratorType   baseIt(   baseImage,  baseImage->GetRequestedRegion());
    BufferIteratorType   net_It(   net_Image,  net_Image->GetRequestedRegion());
   
    realIt.SetFirstDirection(  0 );   realIt.SetSecondDirection( 1 );
    imagIt.SetFirstDirection(  0 );   imagIt.SetSecondDirection( 1 );
    baseIt.SetFirstDirection(  0 );   baseIt.SetSecondDirection( 1 );
    net_It.SetFirstDirection(  0 );   net_It.SetSecondDirection( 1 );
   
    // Software Guide : EndCodeSnippet 
    // set iterators to the beginning
    realIt.GoToBegin();
    imagIt.GoToBegin(); 
    baseIt.GoToBegin();
    net_It.GoToBegin();
    
    // create base phase image and initialize output buffer
    while( !baseIt.IsAtEnd() )
      {
      while ( !baseIt.IsAtEndOfSlice() )
        {
        while ( !baseIt.IsAtEndOfLine() )
          {
          baseIt.Set(  std::atan2(  imagIt.Get() , realIt.Get() ) );
          net_It.Set( 0.0 );
          ++realIt;
          ++imagIt;
          ++baseIt;
          ++net_It;
          }
          realIt.NextLine();
          imagIt.NextLine();
          baseIt.NextLine();
          net_It.NextLine();
        }
        // get next slice
        realIt.NextSlice(); 
        imagIt.NextSlice(); 
        baseIt.NextSlice();
        net_It.NextSlice();
      }
    // output base phase image as a place holder
    std::ostringstream basephase_filename;
    // output file
    basephase_filename << OutputDir <<"/tmap";
    OSSRealzeroright(basephase_filename,4,0,0);
    basephase_filename << "."<< rank <<".mha" ;

    WriteImage(baseImage,basephase_filename,zeroFilterRadius);

  }

  // read input data
  // loop over time instances
  for( int iii = 1 ; iii <= ntime ; iii++)
   {
    // generate list of file names
    RealTimeGenerateFileNames(ExamPath,DirId,iii,nslice,necho,noffset,filenames);

    // get images
    for (int jjj = 0 ; jjj < necho ; jjj ++ )
      {
       // get real images
       reader->SetFileNames(filenames[2*jjj  ] );
       realImage[jjj] = GetImageData(realfilter[jjj], filenames[2*jjj  ] );
       // get imaginary images
       reader->SetFileNames(filenames[2*jjj+1] );
       imagImage[jjj] = GetImageData(imagfilter[jjj], filenames[2*jjj+1] );
      }

    // Software Guide : BeginLatex
    // 
    // We can trigger the reading process by calling the \code{Update()} method
    // on the series reader. It is wise to put this invocation inside a
    // \code{try/catch} block since the process may eventually throw
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
       std::cout << "Error getting echo time " << " (" << echotimekey << ")\n";
       std::cout << "value returned is "<<value << "\n";
       std::cout << e.what() << std::endl;
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
    double tmap_factor =  1.0/(2.0 * pi * imagfreq * alpha * echotime * 1.e-3) ;
    // Software Guide : EndCodeSnippet

   
    // Software Guide : BeginLatex
    // setup real, imaginary, base phase, and temperature map iterators
    // The const slice iterator walks the 3D input image, and the non-const
    // linear iterator walks the 2D output image. The iterators are initialized
    // to walk the same linear path through a slice.  Remember that the
    // \emph{second} direction of the slice iterator defines the direction that
    // linear iteration walks within a slice
    InputIteratorType    realIt(   realImage[0],  realImage[0]->GetRequestedRegion());
    InputIteratorType    imagIt(   imagImage[0],  imagImage[0]->GetRequestedRegion());
    BufferIteratorType   baseIt(   baseImage   ,  baseImage->GetRequestedRegion());
    BufferIteratorType   net_It(   net_Image   ,  net_Image->GetRequestedRegion());
   
    realIt.SetFirstDirection(  0 );   realIt.SetSecondDirection( 1 );
    imagIt.SetFirstDirection(  0 );   imagIt.SetSecondDirection( 1 );
    baseIt.SetFirstDirection(  0 );   baseIt.SetSecondDirection( 1 );
    net_It.SetFirstDirection(  0 );   net_It.SetSecondDirection( 1 );
   
    // Software Guide : EndCodeSnippet 
    // set iterators to the beginning
    realIt.GoToBegin();
    imagIt.GoToBegin(); 
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
             net_It.Get() 
                      +  
             ( std::atan2( imagIt.Get(), realIt.Get() ) - baseIt.Get() ) * tmap_factor
                    );
          // update the base to contain the n-1 image for next time
          baseIt.Set( std::atan2( imagIt.Get(), realIt.Get() ) );
          ++realIt;
          ++imagIt;
          ++baseIt;
          ++net_It;
          }
          realIt.NextLine();
          imagIt.NextLine();
          baseIt.NextLine();
          net_It.NextLine();
        }
      // get next slice
      realIt.NextSlice(); 
      imagIt.NextSlice(); 
      net_It.NextSlice();
      baseIt.NextSlice();
      }
   
    // output unfiltered tmap
    std::ostringstream unfiltered_filename;
    unfiltered_filename << OutputDir <<"/unfiltered";
    OSSRealzeroright(unfiltered_filename,4,0,iii);
    unfiltered_filename << "."<< rank <<".mha" ;

    WriteImage(net_Image , unfiltered_filename , zeroFilterRadius);

    // output filtered tmap
    std::ostringstream filtered_filename;
    filtered_filename << OutputDir <<"/filtered";
    OSSRealzeroright(filtered_filename,4,0,iii);
    filtered_filename << "."<< rank <<".mha" ;

    WriteImage(net_Image , filtered_filename , medianFilterRadius);

   } // end loop over time instances
     
  /* Free memory */
  ierr = DADestroy(dac);CHKERRQ(ierr);
  ierr = PetscFinalize();CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
