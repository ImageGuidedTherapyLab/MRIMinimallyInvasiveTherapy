/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: DicomImageReadChangeHeaderWrite.cxx,v $
  Language:  C++
  Date:      $Date: 2005/11/19 16:31:50 $
  Version:   $Revision: 1.8 $

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
//  This example illustrates how to read a single DICOM slice and write it back
//  with some changed header information as another DICOM slice. Header
//  Key/Value pairs can be specified on the command line. The keys are defined
//  in the file
//
//  \code{Insight/Utilities/gdcm/Dicts/dicomV3.dic}
//  
//  Please note that modifying the content of a DICOM header is a very risky
//  operation. The Header contains fundamental information about the patient
//  and therefore its consistency must be protected from any data corruption.
//  Before attempting to modify the DICOM headers of your files, you must make
//  sure that you have a very good reason for doing so, and that you can ensure
//  that this information change will not result in a lower quality of health
//  care to be delivered to the patient.
//
//  \index{DICOM!Changing Headers}
//
//  Software Guide : EndLatex 



// Software Guide : BeginLatex
// 
// We must start by including the relevant header files. Here we include the
// image reader, image writer, the image, the Meta data dictionary and its
// entries the Meta data objects and the GDCMImageIO. The Meta data dictionary
// is the data container that stores all the entries from the DICOM header once
// the DICOM image file is read into an ITK image.
//
// Software Guide : EndLatex

// Software Guide : BeginCodeSnippet
#include "itkImageLinearIteratorWithIndex.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImage.h"
#include "itkMetaDataDictionary.h"
#include "itkMetaDataObject.h"
#include "itkGDCMImageIO.h"
// Software Guide : EndCodeSnippet

typedef itk::MetaDataDictionary   DictionaryType;

void printdictionary(const DictionaryType   & );

#include <list>
#include <fstream>

int main(int argc, char* argv[])
{

  if( argc < 4 )
    {
    std::cerr << "Usage: " << argv[0] << " [Use Header of this File] [Use Data From This File] [Output Image]\n";
    return EXIT_FAILURE;
    }




// Software Guide : BeginLatex
// 
// We declare the image type by selecting a particular pixel type and image
// dimension.
//
// Software Guide : EndLatex

// Software Guide : BeginCodeSnippet
  typedef signed short InputPixelType;
  const unsigned int   Dimension = 2;
  typedef itk::Image< InputPixelType, Dimension > InputImageType;
// Software Guide : EndCodeSnippet



// Software Guide : BeginLatex
// 
// We instantiate the reader type by using the image type as template
// parameter. An instance of the reader is created and the file name to be read
// is taken from the command line arguments.
//
// Software Guide : EndLatex

// Software Guide : BeginCodeSnippet
  typedef itk::ImageFileReader< InputImageType > ReaderType;
  ReaderType::Pointer originread = ReaderType::New();
  ReaderType::Pointer modifyread = ReaderType::New();
  originread->SetFileName( argv[1] );
  modifyread->SetFileName( argv[2] );
// Software Guide : EndCodeSnippet



// Software Guide : BeginLatex
//
// The GDCMImageIO object is created in order to provide the services for
// reading and writing DICOM files. The newly created image IO class is
// connected to the reader.
//
// Software Guide : EndLatex
 
// Software Guide : BeginCodeSnippet
  typedef itk::GDCMImageIO           ImageIOType;
  ImageIOType::Pointer origingdcmImageIO = ImageIOType::New();
  ImageIOType::Pointer modifygdcmImageIO = ImageIOType::New();
  originread->SetImageIO( origingdcmImageIO );
  modifyread->SetImageIO( modifygdcmImageIO );
// Software Guide : EndCodeSnippet




// Software Guide : BeginLatex
// 
// The reading of the image is triggered by invoking \code{Update()} in the
// reader.
//
// Software Guide : EndLatex


  try
    {
// Software Guide : BeginCodeSnippet
    originread->Update();
// Software Guide : EndCodeSnippet
    }
  catch (itk::ExceptionObject & e)
    {
    std::cerr << "exception in file reader " << std::endl;
    std::cerr << e.GetDescription() << std::endl;
    std::cerr << e.GetLocation() << std::endl;
    return EXIT_FAILURE;
    }


  try
    {
// Software Guide : BeginCodeSnippet
    modifyread->Update();
// Software Guide : EndCodeSnippet
    }
  catch (itk::ExceptionObject & e)
    {
    std::cerr << "exception in file reader " << std::endl;
    std::cerr << e.GetDescription() << std::endl;
    std::cerr << e.GetLocation() << std::endl;
    return EXIT_FAILURE;
    }

//// Software Guide : BeginCodeSnippet
//  typedef itk::ImageLinearIteratorWithIndex< InputImageType >    IteratorType;
//  IteratorType originIt( originread->GetOutput(), 
//                         originread->GetOutput()->GetRequestedRegion() );
//  IteratorType modifyIt( modifyread->GetOutput(), 
//                         modifyread->GetOutput()->GetRequestedRegion() );
//
//  originIt.SetDirection(0);
//  modifyIt.SetDirection(0);
//// Software Guide : EndCodeSnippet
//
//// Software Guide: BeginLatex
////
//// Each line in the input is copied to the output.  The input iterator moves
//// forward across columns while the output iterator moves backwards.
////
//// Software Guide : EndLatex
//
//// Software Guide : BeginCodeSnippet
//  for ( originIt.GoToBegin(); ! originIt.IsAtEnd(); originIt.NextLine())
//    {
//    originIt.GoToBeginOfLine();
//    while ( ! originIt.IsAtEndOfLine() )
//      {
//      originIt.Set( modifyIt.Get() );
//      ++originIt;
//      ++modifyIt;
//      }
//    }

// Software Guide : BeginLatex
// 
// Now that the Dictionary has been updated, we proceed to save the image. This
// output image will have the modified data associated to its DICOM header.
//
// Using the image type, we instantiate a writer type and construct a writer.
// A short pipeline between the reader and the writer is connected. The
// filename to write is taken from the command line arguments. The image IO
// object is connected to the writer.
//
// Software Guide : EndLatex

// Software Guide : BeginCodeSnippet
  typedef itk::ImageFileWriter< InputImageType >  Writer1Type;

  Writer1Type::Pointer writer1 = Writer1Type::New();

  writer1->SetFileName( argv[3] );
  writer1->SetInput( modifyread->GetOutput() );
  writer1->SetImageIO( origingdcmImageIO );
  writer1->SetMetaDataDictionary( origingdcmImageIO->GetMetaDataDictionary() );
// Software Guide : EndCodeSnippet


// Software Guide : BeginLatex
// 
// Execution of the writer is triggered by invoking the \code{Update()} method.
//
// Software Guide : EndLatex

  try
    {
// Software Guide : BeginCodeSnippet
    writer1->Update();
// Software Guide : EndCodeSnippet
    }
  catch (itk::ExceptionObject & e)
    {
    std::cerr << "exception in file writer " << std::endl;
    std::cerr << e.GetDescription() << std::endl;
    std::cerr << e.GetLocation() << std::endl;
    return EXIT_FAILURE;
    }



// Software Guide : BeginLatex
// 
// Remember again, that modifying the header entries of a DICOM file involves
// very serious risks for patients and therefore must be done with extreme
// caution.
//
// Software Guide : EndLatex


  return EXIT_SUCCESS;

}


void printdictionary(const DictionaryType   & dictionary){
// Software Guide : BeginLatex
// 
// Since we are interested only in the DICOM tags that can be expressed in
// strings, we declare a MetaDataObject suitable for managing strings.
//
// Software Guide : EndLatex

// Software Guide : BeginCodeSnippet
  typedef itk::MetaDataObject< std::string > MetaDataStringType;
// Software Guide : EndCodeSnippet





// Software Guide : BeginLatex
// 
// We instantiate the iterators that will make possible to walk through all the
// entries of the MetaDataDictionary.
//
// Software Guide : EndLatex

// Software Guide : BeginCodeSnippet
  DictionaryType::ConstIterator itr = dictionary.Begin();
  DictionaryType::ConstIterator end = dictionary.End();
// Software Guide : EndCodeSnippet




//  Software Guide : BeginLatex
//
// For each one of the entries in the dictionary, we check first if its element
// can be converted to a string, a \code{dynamic\_cast} is used for this purpose.
//
//  Software Guide : EndLatex 

// Software Guide : BeginCodeSnippet
  while( itr != end )
    {
    itk::MetaDataObjectBase::Pointer  entry = itr->second;

    MetaDataStringType::Pointer entryvalue = 
      dynamic_cast<MetaDataStringType *>( entry.GetPointer() ) ;
// Software Guide : EndCodeSnippet

    
// Software Guide : BeginLatex
//
// For those entries that can be converted, we take their DICOM tag and pass it
// to the \code{GetLabelFromTag()} method of the GDCMImageIO class. This method
// checks the DICOM dictionary and returns the string label associated to the
// tag that we are providing in the \code{tagkey} variable. If the label is
// found, it is returned in \code{labelId} variable. The method itself return
// false if the tagkey is not found in the dictionary.  For example "0010|0010"
// in \code{tagkey} becomes "Patient's Name" in \code{labelId}.
//
// Software Guide : EndLatex

// Software Guide : BeginCodeSnippet
    if( entryvalue )
      {
      std::string tagkey   = itr->first;
      std::string labelId;
      bool found =  itk::GDCMImageIO::GetLabelFromTag( tagkey, labelId );
// Software Guide : EndCodeSnippet

// Software Guide : BeginLatex
// 
// The actual value of the dictionary entry is obtained as a string with the
// \code{GetMetaDataObjectValue()} method.
// 
// \index{MetaDataObject!GetMetaDataObjectValue()}
// 
// Software Guide : EndLatex
       
// Software Guide : BeginCodeSnippet
      std::string tagvalue = entryvalue->GetMetaDataObjectValue();
// Software Guide : EndCodeSnippet

// Software Guide : BeginLatex
// 
// At this point we can print out an entry by concatenating the DICOM Name or
// label, the numeric tag and its actual value.
//
// Software Guide : EndLatex

// Software Guide : BeginCodeSnippet
      if( found )
        {
        std::cout << "(" << tagkey << ") " << labelId;
        std::cout << " = " << tagvalue.c_str() << std::endl;
        }
// Software Guide : EndCodeSnippet
      else
        {
        std::cout << "(" << tagkey <<  ") " << "Unknown";
        std::cout << " = " << tagvalue.c_str() << std::endl;
        }
      }

// Software Guide : BeginLatex
// 
// Finally we just close the loop that will walk through all the Dictionary
// entries.
//
// Software Guide : EndLatex

// Software Guide : BeginCodeSnippet
    ++itr;
    }
// Software Guide : EndCodeSnippet
  
}


