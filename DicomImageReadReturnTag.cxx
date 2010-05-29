/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: DicomImageReadPrintTags.cxx,v $
  Language:  C++
  Date:      $Date: 2009-03-17 20:36:50 $
  Version:   $Revision: 1.20 $

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
//  It is often valuable to be able to query the entries from the header of a
//  DICOM file. This can be used for checking for consistency, or simply for
//  verifying that we have the correct dataset in our hands.  This example
//  illustrates how to read a DICOM file and then print out most of the DICOM
//  header information. The binary fields of the DICOM header are skipped.
//
//  \index{DICOM!Header}
//  \index{DICOM!Tags}
//  \index{DICOM!Printing Tags}
//  \index{DICOM!Dictionary}
//  \index{DICOM!GDCM}
//  \index{GDCM!Dictionary}
//
//  Software Guide : EndLatex 

// Software Guide : BeginLatex
// 
// The headers of the main classes involved in this example are specified
// below. They include the image file reader, the GDCM image IO object, the
// Meta data dictionary and its entry element the Meta data object. 
//
// \index{MetaDataDictionary!header}
// \index{MetaDataObject!header}
// \index{GDCMImageIO!header}
//
// Software Guide : EndLatex

// system includes
#include <algorithm>
#include <cctype>
#include <iostream>
#include <string>
#include <vector>

// Software Guide : BeginCodeSnippet
#include "itkImageFileReader.h"
#include "itkGDCMImageIO.h"
#include "itkImageIOBase.h"
#include "itkMetaDataDictionary.h"
#include "itkMetaDataObject.h"
// Software Guide : EndCodeSnippet

#include "gdcmGlobal.h"

#if GDCM_MAJOR_VERSION < 2
#include "gdcmDictSet.h"
#endif

// Software Guide : BeginLatex
// Software Guide : EndLatex

//Return the Dicom TagKey value from InputFile using the DicomDictionary
std::string GetDicomTag(const std::string &InputFile, 
                        const std::string &TagKey, 
                        const std::string &DicomDictionary)
{
  // Software Guide : BeginLatex
  // 
  //  We instantiate the type to be used for storing the image once it is read
  //  into memory.
  //
  // Software Guide : EndLatex

  // Software Guide : BeginCodeSnippet
  typedef signed short       PixelType;
  const unsigned int         Dimension = 2;
  
  typedef itk::Image< PixelType, Dimension >      ImageType;
  // Software Guide : EndCodeSnippet


  // Software Guide : BeginLatex
  // 
  // For a full description of the DICOM dictionary please look at the file.
  //
  // \code{Insight/Utilities/gdcm/Dicts/dicomV3.dic}
  //
  // Software Guide : EndLatex
  if( DicomDictionary.size() )
    {
#if GDCM_MAJOR_VERSION < 2
    // shared lib is loaded so dictionary should have been initialized
    gdcm::Global::GetDicts()->GetDefaultPubDict()->AddDict( DicomDictionary );
#else
    // Newer API, specify a path where XML dicts can be found (Part 3/4 & 6)
    gdcm::Global::GetInstance().Prepend( itksys::SystemTools::GetFilenamePath(
                                                     DicomDictionary).c_str() );
    // Load them !
    gdcm::Global::GetInstance().LoadResourcesFiles();
#endif
    }

  // Software Guide : BeginLatex
  // 
  // Using the image type as template parameter we instantiate the type of the
  // image file reader and construct one instance of it.
  //
  // Software Guide : EndLatex

  // Software Guide : BeginCodeSnippet
  typedef itk::ImageFileReader< ImageType >     ReaderType;

  ReaderType::Pointer reader = ReaderType::New();
  // Software Guide : EndCodeSnippet

  // Software Guide : BeginLatex
  // 
  // The GDCM image IO type is declared and used for constructing one image IO
  // object.
  //
  // Software Guide : EndLatex

  // Software Guide : BeginCodeSnippet
  typedef itk::GDCMImageIO       ImageIOType;
  ImageIOType::Pointer dicomIO = ImageIOType::New();
  // Software Guide : EndCodeSnippet

  // Software Guide : BeginLatex
  // 
  // Here we override the gdcm default value of 0xfff with a value of 0xffff
  // to allow the loading of long binary stream in the DICOM file.
  // This is particularly useful when reading the private tag: 0029,1010
  // from Siemens as it allows to completely specify the imaging parameters
  //
  // Software Guide : EndLatex
  dicomIO->SetMaxSizeLoadEntry(0xffff);

  // Software Guide : BeginLatex
  // 
  // We pass to the reader the filename of the image to be read and connect the
  // ImageIO object to it too.
  //
  // Software Guide : EndLatex

  // Software Guide : BeginCodeSnippet
  reader->SetFileName( InputFile );
  reader->SetImageIO( dicomIO );
  // Software Guide : EndCodeSnippet

  // Software Guide : BeginLatex
  // 
  // The reading process is triggered with a call to the \code{Update()} method.
  // This call should be placed inside a \code{try/catch} block because its
  // execution may result in exceptions being thrown.
  //
  // Software Guide : EndLatex

  // exception handling done in through swig interface file 
  //try
  //  {
  //  // Software Guide : BeginCodeSnippet
  reader->Update();
  //  // Software Guide : EndCodeSnippet
  //  }
  //catch (itk::ExceptionObject &ex)
  //  {
  //  // print error message & throw exception again to catch in python
  //  std::cout << ex << std::endl; 
  //  }

  // Software Guide : BeginLatex
  // 
  // Now that the image has been read, we obtain the Meta data dictionary from
  // the ImageIO object using the \code{GetMetaDataDictionary()} method.
  //
  // \index{MetaDataDictionary}
  // \index{GetMetaDataDictionary()!ImageIOBase}
  // \index{ImageIOBase!GetMetaDataDictionary()}
  //
  // Software Guide : EndLatex

  // Software Guide : BeginCodeSnippet
  typedef itk::MetaDataDictionary   DictionaryType;

  const  DictionaryType & dictionary = dicomIO->GetMetaDataDictionary();
  // Software Guide : EndCodeSnippet

  // Software Guide : BeginLatex
  // 
  // Since we are interested only in the DICOM tags that can be expressed in
  // strings, we declare a MetaDataObject suitable for managing strings.
  //
  // Software Guide : EndLatex

  // Software Guide : BeginCodeSnippet
  typedef itk::MetaDataObject< std::string > MetaDataStringType;
  // Software Guide : EndCodeSnippet

  // get value from TagKey
  // Software Guide : BeginCodeSnippet
  std::string labelId;
  std::string value;
  if( itk::GDCMImageIO::GetLabelFromTag( TagKey, labelId ) )
    {
    std::cout << labelId << " (" << TagKey << ")= ";
    if( dicomIO->GetValueFromTag(TagKey, value) )
      {
      std::cout << value ;
      }
    else
      {
      // print error message & throw exception to catch in python
      throw( std::logic_error("No Value Found in File") );
      }
    std::cout << std::endl;
    }
  else
    {
    throw( std::logic_error("Trying to access inexistant DICOM tag.") );
    }
  // Software Guide : EndCodeSnippet

  return value ;
}
/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: Image2.cxx,v $
  Language:  C++
  Date:      $Date: 2009-03-17 21:11:41 $
  Version:   $Revision: 1.20 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
//  Software Guide : BeginLatex
//
//  The first thing required to read an image from a file is to include
//  the header file of the \doxygen{ImageFileReader} class.
//
//  Software Guide : EndLatex 

std::vector<double>  GetPixelValue(const std::string &InputFile, 
                                   const std::vector<int> inputIndex )
{
  // Software Guide : BeginLatex
  //
  // Then, the image type should be defined by specifying the
  // type used to represent pixels and the dimensions of the image.
  //
  // Software Guide : EndLatex 

  // Software Guide : BeginCodeSnippet
  typedef itk::Vector< double, 3 >             PixelType;
  const unsigned int                       Dimension = 3;

  typedef itk::Image< PixelType, Dimension >   ImageType;
  // Software Guide : EndCodeSnippet


  // Software Guide : BeginLatex
  //
  // Using the image type, it is now possible to instantiate the image reader
  // class. The image type is used as a template parameter to define how the
  // data will be represented once it is loaded into memory. This type does
  // not have to correspond exactly to the type stored in the file. However,
  // a conversion based on C-style type casting is used, so the type chosen
  // to represent the data on disk must be sufficient to characterize it
  // accurately. Readers do not apply any transformation to the pixel data
  // other than casting from the pixel type of the file to the pixel type of
  // the ImageFileReader. The following illustrates a typical
  // instantiation of the ImageFileReader type.
  //
  // \index{itk::ImageFileReader!Instantiation}
  // \index{itk::Image!read}
  //
  // Software Guide : EndLatex 

  // Software Guide : BeginCodeSnippet
  typedef itk::ImageFileReader< ImageType >  ReaderType;
  // Software Guide : EndCodeSnippet


  // Software Guide : BeginLatex
  //
  // The reader type can now be used to create one reader object.  A
  // \doxygen{SmartPointer} (defined by the \code{::Pointer} notation) is used
  // to receive the reference to the newly created reader.  The \code{New()}
  // method is invoked to create an instance of the image reader.
  //
  // \index{itk::ImageFileReader!New()}
  // \index{itk::ImageFileReader!Pointer}
  //
  // Software Guide : EndLatex 

  // Software Guide : BeginCodeSnippet
  ReaderType::Pointer reader = ReaderType::New();
  // Software Guide : EndCodeSnippet


  // Software Guide : BeginLatex
  //
  // The minimum information required by the reader is the filename
  // of the image to be loaded in memory. This is provided through
  // the \code{SetFileName()} method. The file format here is inferred
  // from the filename extension. The user may also explicitly specify the
  // data format explicitly using the \doxygen{ImageIO} (See
  // Chapter~\ref{sec:ImagReadWrite} \pageref{sec:ImagReadWrite} for more
  // information
  //
  // \index{itk::ImageFileReader!SetFileName()}
  //
  // Software Guide : EndLatex 

  // Software Guide : BeginCodeSnippet
  reader->SetFileName( InputFile );
  // Software Guide : EndCodeSnippet


  // Software Guide : BeginLatex
  //
  // Reader objects are referred to as pipeline source objects; they
  // respond to pipeline update requests and initiate the data flow in the
  // pipeline. The pipeline update mechanism ensures that the reader only
  // executes when a data request is made to the reader and the reader has
  // not read any data.  In the current example we explicitly invoke the
  // \code{Update()} method because the output of the reader is not connected
  // to other filters. In normal application the reader's output is connected
  // to the input of an image filter and the update invocation on the filter
  // triggers an update of the reader. The following line illustrates how an
  // explicit update is invoked on the reader.
  //
  // \index{itk::ImageFileReader!Update()}
  //
  // Software Guide : EndLatex 

  // Software Guide : BeginCodeSnippet
  reader->Update();
  // Software Guide : EndCodeSnippet


  // Software Guide : BeginLatex
  //
  // Access to the newly read image can be gained by calling the
  // \code{GetOutput()} method on the reader. This method can also be called
  // before the update request is sent to the reader.  The reference to the
  // image will be valid even though the image will be empty until the reader
  // actually executes.
  //
  // \index{itk::ImageFileReader!GetOutput()}
  //
  // Software Guide : EndLatex 

  // Software Guide : BeginCodeSnippet
  ImageType::Pointer image = reader->GetOutput();
  // Software Guide : EndCodeSnippet

  // The individual position of a pixel inside the image is identified by a
  // unique index. An index is an array of integers that defines the position
  // of the pixel along each coordinate dimension of the image. The IndexType
  // is automatically defined by the image and can be accessed using the
  // scope operator like \doxygen{Index}. The length of the array will match
  // the dimensions of the associated image.
  //
  // The following code illustrates the declaration of an index variable and
  // the assignment of values to each of its components.  Please note that
  // \code{Index} does not use SmartPointers to access it. This is because
  // \code{Index} is a light-weight object that is not intended to be shared
  // between objects. It is more efficient to produce multiple copies of
  // these small objects than to share them using the SmartPointer
  // mechanism.
  // 
  // The following lines declare an instance of the index type and initialize
  // its content in order to associate it with a pixel position in the image.
  //
  // Software Guide : EndLatex 

  // Software Guide : BeginCodeSnippet
  ImageType::IndexType pixelIndex;
 
  pixelIndex[0] = inputIndex[0];   // x position
  pixelIndex[1] = inputIndex[1];   // y position
  pixelIndex[2] = inputIndex[2];   // z position
  // Software Guide : EndCodeSnippet


  // Software Guide : BeginLatex
  //
  // Having defined a pixel position with an index, it is then possible to
  // access the content of the pixel in the image.  The \code{GetPixel()}
  // method allows us to get the value of the pixels.
  //
  // \index{itk::Image!GetPixel()}
  // 
  // Software Guide : EndLatex 

  // Software Guide : BeginCodeSnippet
  ImageType::PixelType   tmpValue = image->GetPixel( pixelIndex );

  // Software Guide : EndCodeSnippet
  std::vector<double>  pixelValue(3,0.0);
  pixelValue[0] = tmpValue[0];
  pixelValue[1] = tmpValue[1];
  pixelValue[2] = tmpValue[2];

  return pixelValue;
}
