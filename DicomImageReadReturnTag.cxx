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
