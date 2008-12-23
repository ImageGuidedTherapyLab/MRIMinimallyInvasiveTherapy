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


// globals 
const unsigned int          Dimension = 3;
const double                       pi = 3.1415926535897931;
// useful typedefs
typedef double                                                 InputPixelType;
typedef unsigned char                                          OutputPixelType;
typedef itk::Image<  InputPixelType, Dimension >               InputImageType;
typedef itk::Image< OutputPixelType, Dimension >              OutputImageType;
typedef itk::ImageSeriesReader<      InputImageType >              ReaderType;
typedef itk::ImageFileWriter<       OutputImageType >              WriterType;
typedef itk::GDCMImageIO                                          ImageIOType;
typedef itk::ImageSliceConstIteratorWithIndex< InputImageType >  
                                                            InputIteratorType;
typedef itk::ImageSliceIteratorWithIndex< InputImageType > BufferIteratorType;
typedef itk::ImageSliceIteratorWithIndex<OutputImageType > OutputIteratorType;

// subroutine to generate dicom filenames
std::vector<std::string>   RealTimeGenerateFileNames(const char * ExamPath,
                                                     int DirId, 
                                                     int istep, int nslice)
{
   std::vector<std::string>   filenames(2*nslice);
   for( int jjj = 0 ; jjj < 2*nslice ; jjj++)
    {
       std::ostringstream file_name;
       // hard code for now
       file_name << ExamPath << "/s" << DirId << "/i" 
                             <<   DirId     + 2*nslice*istep + (jjj+1)
                             << ".MRDC."<<    2*nslice*istep + (jjj+1) ;
       filenames[jjj] = file_name.str();
    }
   return filenames;
}
// subroutine to read the dicom data in real-time
InputImageType::ConstPointer GetImageData(ReaderType::Pointer  reader ,
                                     std::vector<std::string> &filenames )
{
  bool DataNotRead = true;
  InputImageType::ConstPointer Image;
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
       reader->Update();
       DataNotRead = false;
       Image = reader->GetOutput();
       }
     catch (itk::ExceptionObject &excp)
       {
       std::cout << "Exception thrown" << excp << std::endl;
       sleep(1);
       }
    }
  return Image; 
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
void WriteImage(InputImageType::Pointer Image, BufferIteratorType imageIt,
                                             std::ostringstream &filename,
                                          InputImageType::SpacingType &sp,
                                          InputImageType::PointType &orgn)
{
  // create output image buffer
  OutputImageType::Pointer outputImage = OutputImageType::New();
  outputImage->SetRegions( Image->GetRequestedRegion());
  outputImage->Allocate();
  outputImage->SetSpacing( sp );  
  outputImage->SetOrigin( orgn);
  OutputIteratorType  outputIt( outputImage, outputImage->GetRequestedRegion());
  outputIt.SetFirstDirection(  0 );   outputIt.SetSecondDirection( 1 );
  // copy image buffer into an output buffer of a possibly different data type
   imageIt.GoToBegin();
  outputIt.GoToBegin();
  while( !imageIt.IsAtEnd() )
    {
    while ( !imageIt.IsAtEndOfSlice() )
      {
      while ( !imageIt.IsAtEndOfLine() )
        {
        outputIt.Set( (OutputPixelType)  imageIt.Get()  );
        ++imageIt;
        ++outputIt;
        }
        imageIt.NextLine();
        outputIt.NextLine();
      }
      imageIt.NextSlice();
      outputIt.NextSlice();
    }
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( filename.str() );
  std::cout << "writing " << filename.str() << "\n"; 
  writer->SetInput( outputImage );
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
// main driver routine
int main( int argc, char* argv[] )
{
  if( argc < 8 )
    {
    std::cerr << "Usage: " << argv[0] << 
      " ExamPath DirId OutputDir NumberTimePoints NumberSlice Alpha (ppm/degC) Maxdiff (degC) [EchoTime (ms)] [ImagingFreqency (MHz)]" << std::endl;
    std::cerr << "EchoTime and ImagingFrequency Optional\n" 
              << " but may be input in case of any errors\n" ;
    return EXIT_FAILURE;
    }

  // read input data
  const char * ExamPath     = argv[1];
  const char * OutputDir    = argv[3];
  int DirId;    
  int ntime;    
  int nslice;   
  double alpha,maxdiff; 
  try
  {  
     DirId   = boost::lexical_cast<int>(   argv[2]);
     ntime   = boost::lexical_cast<int>(   argv[4]);
     nslice  = boost::lexical_cast<int>(   argv[5]);
     alpha   = boost::lexical_cast<double>(argv[6]);
     maxdiff = boost::lexical_cast<double>(argv[7]);
  }
  catch(const std::exception& e) //catch bad lexical cast
  {
     std::cout << "Error getting Command Line Values:\n"
               << "\nDirId  = " << argv[2]
               << "\nntime  = " << argv[4]
               << "\nnslice = " << argv[5]
               << "\nalpha  = " << argv[6]
               << "\nmaxdiff= " << argv[7] << std::endl;
     std::cout << e.what() << std::endl;
     return EXIT_FAILURE;
  }

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
  ReaderType::Pointer       base_reader    = ReaderType::New();
  ImageIOType::Pointer      base_gdcmIO    = ImageIOType::New();
  InputImageType::Pointer   baseImage = InputImageType::New();
  InputImageType::Pointer   net_Image = InputImageType::New();
  InputImageType::Pointer   filtImage = InputImageType::New();
  // dimensions for all images
  InputImageType::SpacingType sp;
  InputImageType::PointType orgn;
  {
    // generate list of file names and get image data
    std::vector<std::string> filenames =
                             RealTimeGenerateFileNames(ExamPath,DirId,0,nslice);
    base_reader->SetImageIO(  base_gdcmIO );
    base_reader->SetFileNames(  filenames );
    InputImageType::ConstPointer inputImage=GetImageData(base_reader,filenames);

    // allocate base phase image and output image buffer
    InputImageType::RegionType region;
    InputImageType::RegionType::SizeType size;
    InputImageType::RegionType::IndexType index;
    InputImageType::RegionType requestedRegion=inputImage->GetRequestedRegion();
    index = requestedRegion.GetIndex();
    size  = requestedRegion.GetSize();
    size[2]  = nslice ; 
    region.SetSize( size );
    region.SetIndex( index );
    baseImage->SetRegions( region );
    net_Image->SetRegions( region );
    filtImage->SetRegions( region );
    baseImage->Allocate();
    net_Image->Allocate();
    filtImage->Allocate();

    // scale image to meters
    sp = 0.001 * inputImage->GetSpacing();
    std::cout << "Spacing = ";
    std::cout << sp[0] << ", " << sp[1] << ", " << sp[2] << std::endl;
  
    // scale image to meters
    const InputImageType::PointType& orgn_mm = inputImage->GetOrigin();
    orgn[0] = 0.001 * orgn_mm[0];
    orgn[1] = 0.001 * orgn_mm[1];
    orgn[2] = 0.001 * orgn_mm[2];
    std::cout << "Origin = ";
    std::cout << orgn[0] << ", " << orgn[1] << ", " << orgn[2] << std::endl;

    // write out header info
    WriteIni(OutputDir,orgn,sp,size);

    // setup real, imaginary, base phase, and temperature map iterators
    InputIteratorType    realIt(  inputImage, inputImage->GetRequestedRegion());
    InputIteratorType    imagIt(  inputImage, inputImage->GetRequestedRegion());
    BufferIteratorType   baseIt(   baseImage,  baseImage->GetRequestedRegion());
    BufferIteratorType   net_It(   net_Image,  net_Image->GetRequestedRegion());
    BufferIteratorType   filtIt(   filtImage,  filtImage->GetRequestedRegion());
   
    realIt.SetFirstDirection(  0 );   realIt.SetSecondDirection( 1 );
    imagIt.SetFirstDirection(  0 );   imagIt.SetSecondDirection( 1 );
    baseIt.SetFirstDirection(  0 );   baseIt.SetSecondDirection( 1 );
    net_It.SetFirstDirection(  0 );   net_It.SetSecondDirection( 1 );
    filtIt.SetFirstDirection(  0 );   filtIt.SetSecondDirection( 1 );
   
    // Software Guide : EndCodeSnippet 
    realIt.GoToBegin();
    // imaginary image is after the real image
    imagIt.GoToBegin(); imagIt.NextSlice(); // n+1
    // set base and output image to the beginnning
    baseIt.GoToBegin();
    net_It.GoToBegin();
    filtIt.GoToBegin();
    
    // create base phase image and initialize output buffer
    while( !baseIt.IsAtEnd() )
      {
      while ( !baseIt.IsAtEndOfSlice() )
        {
        while ( !baseIt.IsAtEndOfLine() )
          {
          baseIt.Set(  std::atan2(  imagIt.Get() , realIt.Get() ) );
          net_It.Set( 0.0 );
          filtIt.Set( 0.0 );
          ++realIt;
          ++imagIt;
          ++baseIt;
          ++net_It;
          ++filtIt;
          }
          realIt.NextLine();
          imagIt.NextLine();
          baseIt.NextLine();
          net_It.NextLine();
          filtIt.NextLine();
        }
      // skip two on the input image
      realIt.NextSlice(); realIt.NextSlice();
      imagIt.NextSlice(); imagIt.NextSlice();
      // skip one on the base and output image
        baseIt.NextSlice();
        net_It.NextSlice();
        filtIt.NextSlice();
      }
    // output base phase image as a place holder
    std::ostringstream basephase_filename;
    // output file
    basephase_filename << OutputDir <<"/tmap."<< 0 <<".mha" ;
    WriteImage(baseImage,baseIt,basephase_filename,sp,orgn);

  }
  // loop over time instances
  for( int iii = 1 ; iii <= ntime ; iii++)
   {
    std::vector<std::string> filenames =
                          RealTimeGenerateFileNames(ExamPath,DirId,iii,nslice);
   
    // Software Guide : BeginLatex
    // We construct one instance of the series reader object. Set the DICOM
    // image IO object to be use with it, and set the list of filenames to
    // read.
    // Software Guide : EndLatex
    
    // Software Guide : BeginCodeSnippet
    ReaderType::Pointer  reader = ReaderType::New();
    ImageIOType::Pointer gdcmIO = ImageIOType::New();
    reader->SetImageIO( gdcmIO );
    reader->SetFileNames( filenames );
    // Software Guide : EndCodeSnippet
   
    // Software Guide : BeginLatex
    // 
    // We can trigger the reading process by calling the \code{Update()} method
    // on the series reader. It is wise to put this invocation inside a
    // \code{try/catch} block since the process may eventually throw
    // exceptions.
    //
    // Software Guide : EndLatex
   
    InputImageType::ConstPointer inputImage = GetImageData(reader,filenames);
 
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
       std::cout << "using value from command line "<< argv[8] << "\n";
       echotime = boost::lexical_cast<double>(argv[8]);
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
       std::cout << "using value from command line "<< argv[9] << "\n";
       imagfreq = boost::lexical_cast<double>(argv[9]);
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
    InputIteratorType    realIt(  inputImage, inputImage->GetRequestedRegion());
    InputIteratorType    imagIt(  inputImage, inputImage->GetRequestedRegion());
    BufferIteratorType   baseIt(   baseImage,  baseImage->GetRequestedRegion());
    BufferIteratorType   net_It(   net_Image,  net_Image->GetRequestedRegion());
   
    realIt.SetFirstDirection(  0 );   realIt.SetSecondDirection( 1 );
    imagIt.SetFirstDirection(  0 );   imagIt.SetSecondDirection( 1 );
    baseIt.SetFirstDirection(  0 );   baseIt.SetSecondDirection( 1 );
    net_It.SetFirstDirection(  0 );   net_It.SetSecondDirection( 1 );
   
    // Software Guide : EndCodeSnippet 
    realIt.GoToBegin();
    // imaginary image is after the real image
    imagIt.GoToBegin(); imagIt.NextSlice(); // n+1
    // set base and output image to the beginnning
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
             //net_It.Get() 
             //         +  
             (std::atan2( imagIt.Get(), realIt.Get() ) - baseIt.Get()) * tmap_factor
                    );
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
      // skip two on the input image
      realIt.NextSlice(); realIt.NextSlice();
      imagIt.NextSlice(); imagIt.NextSlice();
      // skip one on the base and output image
      net_It.NextSlice();
      baseIt.NextSlice();
      }
   
    // output unfiltered tmap
    std::ostringstream unfiltered_filename;
    unfiltered_filename << OutputDir <<"/unfiltered"<< iii <<".mha" ;
    WriteImage(net_Image , net_It, unfiltered_filename,sp,orgn);

//    // write filtered imaged
//    BufferIteratorType   filtIt(   filtImage,  filtImage->GetRequestedRegion());
//    filtIt.SetFirstDirection(  0 );   filtIt.SetSecondDirection( 1 );
//    filtIt.GoToBegin();
//    net_It.GoToBegin();
//    
//    while( !net_It.IsAtEnd() )
//      {
//      while ( !net_It.IsAtEndOfSlice() )
//        {
//        while ( !net_It.IsAtEndOfLine() )
//          {
//          // only set the temperature if it has changed within a physically
//          // meaningful value. eventually, maxdiff should be predicted from the
//          // kalman filter
//          if ( std::abs( filtIt.Get() - net_It.Get() ) < maxdiff ) filtIt.Set( net_It.Get() );
//          ++net_It;
//          ++filtIt;
//          }
//          filtIt.NextLine();
//          net_It.NextLine();
//        }
//      net_It.NextSlice();
//      filtIt.NextSlice();
//      }
//    // output filtered tmap
//    std::ostringstream filtered_filename;
//    filtered_filename << OutputDir <<"/filtered"<< iii <<".mha" ;
//    WriteImage(filtImage , filtIt, filtered_filename,sp,orgn);
//
   } // end loop over time instances
     
  return EXIT_SUCCESS;
}
