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

#ifdef __BORLANDC__
#define ITK_LEAN_AND_MEAN
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

//system includes
#include <vector>
#include <boost/lexical_cast.hpp>

// Software Guide : BeginCodeSnippet
#include <iostream>
#include <sstream>
#include <string>
#include "itkGDCMImageIO.h"
#include "itkGDCMSeriesFileNames.h"
#include "itkImageSeriesReader.h"
#include "itkImageSeriesWriter.h"
#include "itkImageSliceConstIteratorWithIndex.h"
#include "itkImageSliceIteratorWithIndex.h"
#include "itkScalarToArrayCastImageFilter.h"
#include "itkMagnitudeAndPhaseToComplexImageFilter.h"
#include "itkRealAndImaginaryToComplexImageFilter.h"
#include "itkComplexToModulusImageFilter.h"
#include "itkComplexToPhaseImageFilter.h"
#include "itkComplexToRealImageFilter.h"
#include "itkComplexToImaginaryImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include <itksys/SystemTools.hxx>
#include "gdcmGlobal.h"

#if GDCM_MAJOR_VERSION < 2
#include "gdcmDictSet.h"
#endif
// Software Guide : EndCodeSnippet


//local includes
#include <GetPot>
#include <itkVTKImageVariableNameIO.h>
//  Typedef the image types
typedef double              PixelType;
const unsigned int      Dimension = 3;

typedef itk::Image<              PixelType  , Dimension >        ImageType;
typedef itk::Image< std::complex<PixelType> , Dimension > ComplexImageType;
typedef itk::ImageSeriesReader< ImageType >                     ReaderType;
typedef itk::GDCMImageIO                                       ImageIOType;
//typedef itk::Image< itk::Vector< PixelType , 2 >,
//                                           Dimension > ComplexImageType;
//typedef itk::ScalarToArrayCastImageFilter< ImageType, ComplexImageType >
//                                                      ComplexFilterType;
typedef itk::RealAndImaginaryToComplexImageFilter< 
                   PixelType, PixelType, PixelType> RealImagComplexFilterType;
typedef itk::MagnitudeAndPhaseToComplexImageFilter< 
                   PixelType, PixelType, PixelType> MagnPhasComplexFilterType;
typedef itk::ComplexToModulusImageFilter< ComplexImageType , ImageType >  
                                                        ComplexMagnFilterType;
typedef itk::ComplexToPhaseImageFilter<  ComplexImageType , ImageType >  
                                                        ComplexPhasFilterType;

// write the real imaginary phase magnitude data
// write header data from gdcm file
void WriteRawDataAndHeader( GetPot  &cmdLine,
                            ComplexImageType::Pointer complexImage,
                            ComplexImageType::Pointer phasorImage ,
                            ImageIOType::Pointer          m_gdcmIO )
{
  // setup filters for components
  typedef itk::ComplexToRealImageFilter< ComplexImageType,ImageType > 
                                                    RealFilterType;
  typedef itk::ComplexToImaginaryImageFilter< ComplexImageType,ImageType >
                                                    ImagFilterType;
  RealFilterType::Pointer   realFilter = RealFilterType::New(),
                            magnFilter = RealFilterType::New();
  ImagFilterType::Pointer   imagFilter = ImagFilterType::New(),
                            phasFilter = ImagFilterType::New();

  realFilter->SetInput( complexImage );
  imagFilter->SetInput( complexImage );
  magnFilter->SetInput(  phasorImage );
  phasFilter->SetInput(  phasorImage );
  try
    {
     realFilter->Update( );
     imagFilter->Update( );
     magnFilter->Update( );
     phasFilter->Update( );
    }
  catch (itk::ExceptionObject &excp)
    {
    std::cout << "Exception thrown" << excp << std::endl; abort();
    }
  // setup writer
  typedef itk::ImageFileWriter<         ImageType >     WriterType;
  WriterType::Pointer   realWriter = WriterType::New(),
                        imagWriter = WriterType::New(),
                        magnWriter = WriterType::New(),
                        phasWriter = WriterType::New();

  // custom class for proper variable name
  itk::VTKImageVariableNameIO::Pointer 
                         ioRealPointer = itk::VTKImageVariableNameIO::New(),
                         ioImagPointer = itk::VTKImageVariableNameIO::New(),
                         ioMagnPointer = itk::VTKImageVariableNameIO::New(),
                         ioPhasPointer = itk::VTKImageVariableNameIO::New();
  // set variable name
  ioRealPointer->SetVariableName( "real" );
  ioImagPointer->SetVariableName( "imag" );
  ioMagnPointer->SetVariableName( "magn" );
  ioPhasPointer->SetVariableName( "phas" );
  // custom class for proper data name
  realWriter->SetImageIO( ioRealPointer );
  imagWriter->SetImageIO( ioImagPointer );
  magnWriter->SetImageIO( ioMagnPointer );
  phasWriter->SetImageIO( ioPhasPointer );

  // setup and write
  const int timeid = cmdLine.follow(0, "-timeid");  
  std::string fileBase = cmdLine.follow("./image","--output");
  realWriter->SetInput( realFilter->GetOutput() );
  imagWriter->SetInput( imagFilter->GetOutput() );
  magnWriter->SetInput( magnFilter->GetOutput() );
  phasWriter->SetInput( phasFilter->GetOutput() );

# define OSSRealzeroright(o,v,p,d)  (o).width(v);  (o).precision(p); (o).fill('0'); (o) << std::right << (d)
  // output file name
  std::ostringstream realFile, 
                     imagFile, 
                     magnFile, 
                     phasFile ; // filename
  realFile << fileBase << "Real." ;
  imagFile << fileBase << "Imag.";
  magnFile << fileBase << "Magn.";
  phasFile << fileBase << "Phas.";
  OSSRealzeroright(realFile,4,0, timeid);
  OSSRealzeroright(imagFile,4,0, timeid);
  OSSRealzeroright(magnFile,4,0, timeid);
  OSSRealzeroright(phasFile,4,0, timeid);
  realFile << ".vtk" ;
  imagFile << ".vtk" ;
  magnFile << ".vtk" ;
  phasFile << ".vtk";
  realWriter->SetFileName( realFile.str() );
  imagWriter->SetFileName( imagFile.str() );
  magnWriter->SetFileName( magnFile.str() );
  phasWriter->SetFileName( phasFile.str() );
  try
    {
    realWriter->Update();
    imagWriter->Update();
    magnWriter->Update();
    phasWriter->Update();
    }
  catch (itk::ExceptionObject &excp)
    {
    std::cout << "Exception thrown" << excp << std::endl; abort();
    }

  // output dicom info as an ini file
  std::ostringstream outputFile ;
  std::ofstream Inifile;

  // output file
  outputFile << fileBase << "Header."  ;
  OSSRealzeroright(outputFile,4,0, timeid);
  outputFile << ".ini";
  Inifile.open(outputFile.str().c_str(), std::ios::out);

  // dicom header parameters...
  Inifile <<"[rawdata]" << std::endl;

  // scratch storage for header value
  std::string value; 
  // Get Echo Time 
  std::string echotimekey = "0018|0081";
  double echotime;
  try
  {  
     m_gdcmIO->GetValueFromTag(echotimekey, value) ;
     echotime=boost::lexical_cast<double>(value);
  }
  catch(const boost::bad_lexical_cast& e) //catch bad lexical cast
  {
    try
    {
      // trailing space on string causing cast error
      m_gdcmIO->GetValueFromTag(echotimekey, value) ;
      value.erase(value.find(' ')) ; 
      echotime=boost::lexical_cast<double>(value);
    }
    catch(const std::exception& e) //catch bad lexical cast
    {
      std::cout<< std::endl
               <<"#########################Error getting echo time " 
               << " (" << echotimekey << ")" << std::endl;
      std::cout<<"value returned is "<<value<< std::endl;
      std::cout<< e.what() << std::endl;

      if(!cmdLine.search(1, "-echotime")) // echo time found on the command line
         {
         echotime  = cmdLine.follow(5., "-echotime");   // [rad]
         std::cout << "using value from command line "<< echotime << std::endl;
         }
      else
         {
         std::cout << "input value from command line -echotime ..."<< std::endl;
         abort();
         }
    }
  }
  /* image frequency ie gyromagnetic ratio */
  Inifile <<"echotime = " << echotime << std::endl ;

  std::string imagfreqkey = "0018|0084";
  double imagfreq;
  try
  {  
     m_gdcmIO->GetValueFromTag(imagfreqkey, value) ;
     imagfreq=boost::lexical_cast<double>(value);
  }
  catch(const boost::bad_lexical_cast& e) //catch bad lexical cast
  {
     try
     {
       // trailing space on string causing cast error
       m_gdcmIO->GetValueFromTag(imagfreqkey, value) ;
       value.erase(value.find(' ')) ; 
       imagfreq=boost::lexical_cast<double>(value);
     }
     catch(const std::exception& e) //catch standard exception
     {
       std::cout << std::endl 
                 << "#################Error getting Imaging Freq"
                 <<" ("<<imagfreqkey << ")" << std::endl;
       std::cout << "value returned is "<<value << std::endl;
       std::cout << e.what() << std::endl;

       if(!cmdLine.search(1, "-imagfreq")) // echo time found on the command line
         {
         echotime  = cmdLine.follow(42.98, "-imagfreq");   // [rad]
         std::cout << "using value from command line "<< imagfreq << std::endl;
         }
       else
          {
          std::cout << "input value from command line -imagfreq..."<< std::endl;
          abort();
          }
     }
  }
  /* image frequency ie gyromagnetic ratio */
  Inifile <<"imagfreq = " << imagfreq << std::endl ;

  // get image acquisition time
  double acquisitionTime;
  try
  {  
     std::string acqDurationIdkey = "0019|105a";
     std::string fastPhaseIdkey   = "0019|10f2" ;
     m_gdcmIO->GetValueFromTag(acqDurationIdkey , value) ;
     // FIXME: ITK returns ascii representation of a float
     //        need to cast as a pointer to float then dereference
     float acqDuration = *(reinterpret_cast<const float*>(value.c_str())) ;
     m_gdcmIO->GetValueFromTag(fastPhaseIdkey   , value) ;
     int fastPhase = boost::lexical_cast<int>(value);
     // convert from \mu second to second
     acquisitionTime = acqDuration / fastPhase * 1.e-6;
  }
  catch(const std::exception& e) //catch bad lexical cast
  {
    std::cout<< std::endl
             <<"#########################Error getting acquisition time" 
             << std::endl;
    std::cout<<"value returned is "<<value<< std::endl;
    std::cout<<e.what() << std::endl;
    if( !cmdLine.search(1, "-deltat")) abort();
    acquisitionTime = cmdLine.follow(5.0, "-deltat");   // [s]
  }
  /* acquisition time */
  Inifile <<"deltat = " << acquisitionTime << std::endl ;

  // get number of slices
  int nslice = 0 ; 
  try
  {  
     std::string nsliceIdkey      = "0021|104f" ;
     m_gdcmIO->GetValueFromTag(nsliceIdkey , value) ;
     nslice = boost::lexical_cast<int>(value);
  }
  catch(const std::exception& e) //catch bad lexical cast
  {
    std::cout<< std::endl
             <<"#########################Error getting nslice" 
             << std::endl;
    std::cout<<"value returned is "<<value<< std::endl;
    std::cout<<e.what() << std::endl; 
    if(!cmdLine.search(1, "-nslice")) abort();
    nslice  = cmdLine.follow(5, "-nslice");  
  }
  /* out of plane dimension */
  Inifile << "nslice = " << nslice << std::endl ;

  // get number of echos
  int necho = 1; 
  try
  {  
     std::string nechoIdkey       = "0019|107e" ;
     m_gdcmIO->GetValueFromTag(nechoIdkey  , value) ;
     necho = boost::lexical_cast<int>(value);
  }
  catch(const std::exception& e) //catch bad lexical cast
  {
    std::cout<< std::endl
             <<"#########################Error getting necho" 
             << std::endl;
    std::cout<<"value returned is "<<value<< std::endl;
    std::cout<<e.what() << std::endl; 
    if(!cmdLine.search(1, "-necho")) abort();
    necho  = cmdLine.follow(1, "-necho");   // [rad]
  }
  /* number of echoes */
  Inifile <<"necho = " << necho << std::endl ;

  // close file 
  Inifile.close();
  return;
}

//convert from real and imaginary data
void ConvertFromComplexValue( GetPot  &cmdLine,
                              ReaderType::FileNamesContainer &realFilenames, 
                              ReaderType::FileNamesContainer &imagFilenames )
{
// Software Guide : BeginLatex
//
//  We also declare types for the \doxygen{GDCMImageIO} object that will
//  actually read and write the DICOM images, and the
//  \doxygen{GDCMSeriesFileNames} object that will generate and order all the
//  filenames for the slices composing the volume dataset. Once we have the
//  types, we proceed to create instances of both objects.
//
// Software Guide : EndLatex

// Software Guide : BeginCodeSnippet
  ImageIOType::Pointer realGdcmIO = ImageIOType::New(),
                       imagGdcmIO = ImageIOType::New();
// Software Guide : BeginLatex
// 
// We construct the series reader object. Set the DICOM image
// IO object to be use with it, and set the list of filenames to read.
//
// Software Guide : EndLatex

  ReaderType::Pointer realReader = ReaderType::New(),
                      imagReader = ReaderType::New();
  realReader->SetImageIO(   realGdcmIO );
  imagReader->SetImageIO(   imagGdcmIO );
  realReader->SetFileNames( realFilenames );
  imagReader->SetFileNames( imagFilenames );


  ImageType::ConstPointer realImage,
                          imagImage;
  try
    {
        realReader->Update();
        realImage = realReader->GetOutput();
        imagReader->Update();
        imagImage = imagReader->GetOutput();
    }
  catch (itk::ExceptionObject &excp)
    {
    std::cerr << "Exception thrown while reading the image" << std::endl;
    std::cerr << excp << std::endl;
    abort();
    }

  // data structures to hold complex images
  RealImagComplexFilterType::Pointer complexFilter
                                     = RealImagComplexFilterType::New(), 
                                     phasorFilter
                                     = RealImagComplexFilterType::New(); 
  ComplexMagnFilterType::Pointer magnitudefilter = ComplexMagnFilterType::New();
  ComplexPhasFilterType::Pointer phasefilter = ComplexPhasFilterType::New();

  // just push back complex filter
  complexFilter->PushBackInput( realImage );
  complexFilter->PushBackInput( imagImage );

  // get the magnitude and phase of the complex data
  magnitudefilter->SetInput( complexFilter->GetOutput() );
  phasefilter->SetInput(     complexFilter->GetOutput() );

  // push back the magnitude and phase
  phasorFilter->PushBackInput( magnitudefilter->GetOutput() );
  phasorFilter->PushBackInput(     phasefilter->GetOutput() );

  // get image data 
  try
    {
    complexFilter->Update();
    phasorFilter->Update();
    }
  catch (itk::ExceptionObject &excp)
    {
    std::cout << "Exception thrown" << excp << std::endl; abort();
    }

  //  write data and header
  WriteRawDataAndHeader(cmdLine,
                        complexFilter->GetOutput(), 
                        phasorFilter->GetOutput(), realGdcmIO );
  return;
}

//convert from real and imaginary data
void ConvertFromPhasorValue(  GetPot  &cmdLine,
                              ReaderType::FileNamesContainer &magnFilenames, 
                              ReaderType::FileNamesContainer &phasFilenames )
{
// Software Guide : BeginLatex
//
//  We also declare types for the \doxygen{GDCMImageIO} object that will
//  actually read and write the DICOM images, and the
//  \doxygen{GDCMSeriesFileNames} object that will generate and order all the
//  filenames for the slices composing the volume dataset. Once we have the
//  types, we proceed to create instances of both objects.
//
// Software Guide : EndLatex

// Software Guide : BeginCodeSnippet
  ImageIOType::Pointer phasGdcmIO = ImageIOType::New(),
                       magnGdcmIO = ImageIOType::New();
// Software Guide : EndCodeSnippet

  ReaderType::Pointer phasReader = ReaderType::New(),
                      magnReader = ReaderType::New();
  phasReader->SetImageIO(   phasGdcmIO );
  magnReader->SetImageIO(   magnGdcmIO );
  phasReader->SetFileNames( phasFilenames );
  magnReader->SetFileNames( magnFilenames );
// Software Guide : EndCodeSnippet

  // rescale phase to -pi and pi
  typedef itk::RescaleIntensityImageFilter< 
               ImageType, ImageType > RescaleFilterType;

  RescaleFilterType::Pointer rescaleFilter = RescaleFilterType::New();

  const double     pi = 3.1415926535897931;
  rescaleFilter->SetOutputMinimum(-1.0*pi );
  rescaleFilter->SetOutputMaximum( 1.0*pi );
  rescaleFilter->SetInput( phasReader->GetOutput() );

// Software Guide : BeginLatex
// 
// We can trigger the reading process by calling the \code{Update()} method on
// the series reader. It is wise to put this invocation inside a
// \code{try/catch} block since the process may eventually throw exceptions.
//
// Software Guide : EndLatex

  ImageType::ConstPointer phasImage,
                          magnImage;
  try
    {
// Software Guide : BeginCodeSnippet
        magnReader->Update();
        magnImage = magnReader->GetOutput();
        rescaleFilter->Update();
        phasImage = rescaleFilter->GetOutput();
// Software Guide : EndCodeSnippet
    }
  catch (itk::ExceptionObject &excp)
    {
    std::cerr << "Exception thrown while reading the image" << std::endl;
    std::cerr << excp << std::endl;
    abort();
    }

  // data structures to hold complex images
  MagnPhasComplexFilterType::Pointer complexFilter
                                     = MagnPhasComplexFilterType::New(); 
  RealImagComplexFilterType::Pointer phasorFilter
                                     = RealImagComplexFilterType::New(); 

  // convert magnitude and phase to real imaginary
  complexFilter->PushBackInput( magnImage );
  complexFilter->PushBackInput( phasImage );

  // write out magnitude and phase 
  phasorFilter->PushBackInput( magnImage );
  phasorFilter->PushBackInput( phasImage );

  //// setup iterators
  //itk::ImageRegionIterator< ComplexImageType > 
  //              complexIt(complexFilter->GetOutput(), 
  //                        complexFilter->GetOutput()->GetRequestedRegion());
  //// Software Guide : EndCodeSnippet
  //// set iterators to the beginning
  //complexIt.GoToBegin();
  //while( !complexIt.IsAtEnd() )
  //  {
  //   itk::Vector<  PixelType , 2 > complexValue; 
  //   // from MagnitudeAndPhaseToComplex
  //   std::complex< PixelType > tmpComplex 
  //        = std::polar( complexIt.Get()[0],complexIt.Get()[1]) ;
  //   complexValue[0] =  tmpComplex.real(); 
  //   complexValue[1] =  tmpComplex.imag();

  //   complexIt.Set( complexValue );
  //   ++complexIt;
  //  }

  // get image data 
  try
    {
    complexFilter->Update();
    phasorFilter->Update();
    }
  catch (itk::ExceptionObject &excp)
    {
    std::cout << "Exception thrown" << excp << std::endl; abort();
    }

  //  write data and header
  WriteRawDataAndHeader(cmdLine,
                        complexFilter->GetOutput(), 
                        phasorFilter->GetOutput(), magnGdcmIO );
  return;
}


// Driver routine
int main( int argc, char* argv[] )
{
  GetPot   cmdLine(argc, argv);
  if( cmdLine.size() == 1 || cmdLine.search(2, "--help", "-h") ) 
   {
    std::cerr << "Usage: " << argv[0] << 
      " DicomDirectory  OutputDicomDirectory" << std::endl;
    std::cerr << "Write none, one, or more of the files as arguments:\n\n";
    std::cerr << "  -Mstring  specify magnitude images as string" << std::endl;
    std::cerr << "  -Rstring  specify real images as string"      << std::endl;
    std::cerr << "  -Istring  specify imaginary images as string" << std::endl;
    std::cerr << "  -Pstring  specify phase images as string"     << std::endl;
    std::cerr << "You can use as many arguments as needed " << std::endl;
    return EXIT_FAILURE;
   }

  // load additional dictionary
  std::string dictionaryFile = cmdLine.follow("","-dictionary");
  if( dictionaryFile.length() )
   { 
#if GDCM_MAJOR_VERSION < 2
     // shared lib is loaded so dictionary should have been initialized
     gdcm::Global::GetDicts()->GetDefaultPubDict()->AddDict(dictionaryFile);
#else
     // Newer API, specify a path where XML dicts can be found (Part 3/4 & 6)
     gdcm::Global::GetInstance().Prepend( itksys::SystemTools::GetFilenamePath(dictionaryFile).c_str() );
     // Load them !
     gdcm::Global::GetInstance().LoadResourcesFiles();
#endif
   }
  // Expecting a python script to generate the list of file names
  // vector of strings for real, imaginary, phase, and magnitude files
  ReaderType::FileNamesContainer realFilenames = cmdLine.string_tails("-R");
  ReaderType::FileNamesContainer imagFilenames = cmdLine.string_tails("-I");
  ReaderType::FileNamesContainer phasFilenames = cmdLine.string_tails("-P");
  ReaderType::FileNamesContainer magnFilenames = cmdLine.string_tails("-M");
  //for(ReaderType::FileNamesContainer::const_iterator it = filenames.begin(); 
  //                                                   it != filenames.end(); it++)
  //    std::cout << "   '" << *it << "'\n";

// Software Guide : EndCodeSnippet

  // error check files entered that is an even number
  if( (realFilenames.size() != imagFilenames.size())
                           ||
      (phasFilenames.size() != magnFilenames.size()) ) 
   {
    std::cout << "# of real and imaginary files should be equal" << std::endl; 
    std::cout << "# of phase and magnitude files should be equal"<< std::endl; 
    return EXIT_FAILURE;
   }
  //for(unsigned int fni = 0; fni<numberOfFilenames; fni++)
  //  {
  //  std::cout << "filename # " << fni << " = ";
  //  std::cout << filenames[fni] << std::endl;
  //  }


  // in either case get the real and image image as well as the phasor image
  if( realFilenames.size() )
   {
    ConvertFromComplexValue(cmdLine,realFilenames,imagFilenames);
   }
  else if( magnFilenames.size() )
   {
    ConvertFromPhasorValue(cmdLine,magnFilenames,phasFilenames);
   }
  else 
   {
    std::cout << "no files entered? " << std::endl << std::flush; abort();
   }

  return EXIT_SUCCESS;
}
