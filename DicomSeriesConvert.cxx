/*=========================================================================


  this code is a modification of:

     Examples/IO/DicomSeriesReadImageWrite2.cxx

  its use is intended to convert a series of DICOM images to .mha format


  It may also be used to apply a contrast enhancement filter to the imaging
  data and write the contrast enhanced filter to an AVS/RAWIV file. The
  Contrast Enhancement is not applied to the .mha file


 
  this must be linked with the  -L/usr/X11R6/lib -lX11 library



  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: DicomSeriesReadImageWrite2.cxx,v $
  Language:  C++
  Date:      $Date: 2007/07/22 05:40:41 $
  Version:   $Revision: 1.11 $

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
//  Probably the most common representation of datasets in clinical
//  applications is the one that uses sets of DICOM slices in order to compose
//  tridimensional images. This is the case for CT, MRI and PET scanners. It is
//  very common therefore for image analysts to have to process volumetric
//  images that are stored in the form of a set of DICOM files belonging to a
//  common DICOM series. 
//
//  The following example illustrates how to use ITK functionalities in order
//  to read a DICOM series into a volume and then save this volume in another
//  file format.
//
//  The example begins by including the appropriate headers. In particular we
//  will need the \doxygen{GDCMImageIO} object in order to have access to the
//  capabilities of the GDCM library for reading DICOM files, and the
//  \doxygen{GDCMSeriesFileNames} object for generating the lists of filenames
//  identifying the slices of a common volumetric dataset.
//
//  \index{itk::ImageSeriesReader!header}
//  \index{itk::GDCMImageIO!header}
//  \index{itk::GDCMSeriesFileNames!header}
//  \index{itk::ImageFileWriter!header}
//
//  Software Guide : EndLatex 

// Software Guide : BeginCodeSnippet
#include <iostream>
#include <fstream>
#include "CImg.h"
#include "ByteSwapping.h"
#include "itkImage.h"
#include "itkImageSliceConstIteratorWithIndex.h"
#include "itkImageSliceIteratorWithIndex.h"
#include "itkGDCMImageIO.h"
#include "itkGDCMSeriesFileNames.h"
#include "itkImageSeriesReader.h"
#include "itkImageFileWriter.h"
#include "SimpleIni.h" // ini file parser

#define RESISTOR  0.95
// Software Guide : EndCodeSnippet

using namespace cimg_library;
using namespace std;

typedef float PixelType;

// used for error checking
void boundcheck(char *variableid,int variable,int upperbnd){
     if(variable < 0 || variable >= upperbnd){
        std::cout << variableid << " out of images bounds"  <<std::endl;
        std::cout << variable << " \\nin [0,"<<upperbnd<<"]"<<std::endl;
        throw;
     }
     return;
};
// function to write avs file
void WriteAVS(char *,CImg<PixelType> *, std::vector< float > *, bool );
// image enhancement 
void ImageEnhance(PixelType *,int ,int );

int main( int argc, char* argv[] )
{

  const char *dicomdir = cimg_option("-d","DicomImages","DicomDirectory");
  const char *inifile  = cimg_option("-c","none","existing control file");
  const char *seriesname = cimg_option("-s","default","seriesName");
  bool image_enhance = cimg_option("-e",true,"Apply Image Enhancement");
  bool endianess;

//Check endianess of machine at run time
//If the bool variable endianess is set to 1, the machine is little endian
//Otherwise it is big endian

  short int word = 0x0001;
  char *byte = (char *) &word;
  if(byte[0]) endianess = 1;
  else endianess = 0;

  std::cout << std::endl;
  if(endianess){std::cout << "native data type is little endian "<< std::endl;
  } else {      std::cout << "native data type is big endian "   << std::endl; }

// Software Guide : BeginLatex
// 
// We define the pixel type and dimension of the image to be read. In this
// particular case, the dimensionality of the image is 3, and we assume a
// \code{signed short} pixel type that is commonly used for X-Rays CT scanners.
// 
// Software Guide : EndLatex

// Software Guide : BeginCodeSnippet
  const unsigned int      Dimension = 3;

  
  typedef itk::Image< PixelType, Dimension >         ImageType;
  typedef itk::Image< PixelType, Dimension >        OutputType;


// Software Guide : EndCodeSnippet

  // output file name
  char file_o[256];


// Software Guide : BeginLatex
// 
// We use the image type for instantiating the type of the series reader and
// for constructing one object of its type.
//
// Software Guide : EndLatex

// Software Guide : BeginCodeSnippet
  typedef itk::ImageSeriesReader< ImageType >        ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
// Software Guide : EndCodeSnippet



// Software Guide : BeginLatex
// 
// A GDCMImageIO object is created and connected to the reader. This object is
// the one that is aware of the internal intricacies of the DICOM format. 
//
// Software Guide : EndLatex

// Software Guide : BeginCodeSnippet
  typedef itk::GDCMImageIO       ImageIOType;
  ImageIOType::Pointer dicomIO = ImageIOType::New();
  
  reader->SetImageIO( dicomIO );
// Software Guide : EndCodeSnippet



// Software Guide : BeginLatex
//
// Now we face one of the main challenges of the process of reading a DICOM
// series. That is, to identify from a given directory the set of filenames
// that belong together to the same volumetric image. Fortunately for us, GDCM
// offers functionalities for solving this problem and we just need to invoke
// those functionalities through an ITK class that encapsulates a communication
// with GDCM classes. This ITK object is the GDCMSeriesFileNames. Conveniently
// for us, we only need to pass to this class the name of the directory where
// the DICOM slices are stored. This is done with the \code{SetDirectory()}
// method. The GDCMSeriesFileNames object will explore the directory and will
// generate a sequence of filenames for DICOM files for one study/series. 
// In this example, we also call the \code{SetUseSeriesDetails(true)} function
// that tells the GDCMSereiesFileNames object to use additional DICOM 
// information to distinguish unique volumes within the directory.  This is
// useful, for example, if a DICOM device assigns the same SeriesID to 
// a scout scan and its 3D volume; by using additional DICOM information
// the scout scan will not be included as part of the 3D volume.  Note that
// \code{SetUseSeriesDetails(true)} must be called prior to calling
// \code{SetDirectory()}. By default \code{SetUseSeriesDetails(true)} will use
// the following DICOM tags to sub-refine a set of files into multiple series:
// * 0020 0011 Series Number
// * 0018 0024 Sequence Name
// * 0018 0050 Slice Thickness
// * 0028 0010 Rows
// * 0028 0011 Columns
// If this is not enough for your specific case you can always add some more
// restrictions using the \code{AddSeriesRestriction()} method. In this example we will use
// the DICOM Tag: 0008 0021 DA 1 Series Date, to sub-refine each series. The format
// for passing the argument is a string containing first the group then the element
// of the DICOM tag, separed by a pipe (|) sign.
//
//
// \index{itk::GDCMSeriesFileNames!SetDirectory()}
//
// Software Guide : EndLatex

// Software Guide : BeginCodeSnippet
  typedef itk::GDCMSeriesFileNames NamesGeneratorType;
  NamesGeneratorType::Pointer nameGenerator = NamesGeneratorType::New();

  nameGenerator->SetUseSeriesDetails( true );
  nameGenerator->AddSeriesRestriction("0008|0021" );

  nameGenerator->SetDirectory( dicomdir );
// Software Guide : EndCodeSnippet
  

  try
    {
    std::cout << std::endl << "The directory: " << std::endl;
    std::cout << std::endl << dicomdir << std::endl << std::endl;
    std::cout << "Contains the following DICOM Series: ";
    std::cout << std::endl << std::endl;


    
// Software Guide : BeginLatex
// 
// The GDCMSeriesFileNames object first identifies the list of DICOM series
// that are present in the given directory. We receive that list in a reference
// to a container of strings and then we can do things like printing out all
// the series identifiers that the generator had found. Since the process of
// finding the series identifiers can potentially throw exceptions, it is
// wise to put this code inside a try/catch block.
//
// Software Guide : EndLatex

// Software Guide : BeginCodeSnippet
    typedef std::vector< std::string >    SeriesIdContainer;
    
    const SeriesIdContainer & seriesUID = nameGenerator->GetSeriesUIDs();
    
    SeriesIdContainer::const_iterator seriesItr = seriesUID.begin();
    SeriesIdContainer::const_iterator seriesEnd = seriesUID.end();
    while( seriesItr != seriesEnd )
      {
      std::cout << seriesItr->c_str() << std::endl;
      seriesItr++;
      }
// Software Guide : EndCodeSnippet
  


// Software Guide : BeginLatex
// 
// Given that it is common to find multiple DICOM series in the same directory,
// we must tell the GDCM classes what specific series do we want to read. In
// this example we do this by checking first if the user has provided a series
// identifier in the command line arguments. If no series identifier has been
// passed, then we simply use the first series found during the exploration of
// the directory.
//
// Software Guide : EndLatex

// Software Guide : BeginCodeSnippet
    std::string seriesIdentifier;

    if(!strcasecmp(seriesname,"default")){
      seriesIdentifier = seriesUID.begin()->c_str();
    } else {
      seriesIdentifier = seriesname;
    }
    
// Software Guide : EndCodeSnippet


    std::cout << std::endl << std::endl;
    std::cout << "Now reading series: " << std::endl << std::endl;
    std::cout << seriesIdentifier << std::endl;
    std::cout << std::endl << std::endl;




// Software Guide : BeginLatex
// 
// We pass the series identifier to the name generator and ask for all the
// filenames associated to that series. This list is returned in a container of
// strings by the \code{GetFileNames()} method. 
//
// \index{itk::GDCMSeriesFileNames!GetFileNames()}
// 
// Software Guide : EndLatex

// Software Guide : BeginCodeSnippet
    typedef std::vector< std::string >   FileNamesContainer;
    FileNamesContainer fileNames;

    fileNames = nameGenerator->GetFileNames( seriesIdentifier );
// Software Guide : EndCodeSnippet




// Software Guide : BeginLatex
// 
//
// The list of filenames can now be passed to the \doxygen{ImageSeriesReader}
// using the \code{SetFileNames()} method.
//  
//  \index{itk::ImageSeriesReader!SetFileNames()}
//
// Software Guide : EndLatex
  
// Software Guide : BeginCodeSnippet
    reader->SetFileNames( fileNames );
// Software Guide : EndCodeSnippet



// Software Guide : BeginLatex
// 
// Finally we can trigger the reading process by invoking the \code{Update()}
// method in the reader. This call as usual is placed inside a \code{try/catch}
// block.
//
// Software Guide : EndLatex

  ImageType::Pointer inputImage;

// Software Guide : BeginCodeSnippet
    try
      {
      reader->Update();
      inputImage = reader->GetOutput();
      }
    catch (itk::ExceptionObject &ex)
      {
      std::cout << ex << std::endl;
      return EXIT_FAILURE;
      }
// Software Guide : EndCodeSnippet


// Software Guide : BeginLatex
// 
// At this point, we have a volumetric image in memory that we can access by
// invoking the \code{GetOutput()} method of the reader.
//
// Software Guide : EndLatex



// Software Guide : BeginLatex
//
//  A slice iterator type is defined to walk the input image.
//
// Software Guide : EndLatex
  
  
  
  // get sizes of the image
  ImageType::RegionType requestedRegion = inputImage->GetRequestedRegion();

  int dim1 =  requestedRegion.GetSize()[0];
  int dim2 =  requestedRegion.GetSize()[1];
  int dim3 =  requestedRegion.GetSize()[2];

  CImg<PixelType> Filter(dim1,dim2,dim3);

  // scale image to meters
  ImageType::SpacingType sp = 0.001 * inputImage->GetSpacing();
  inputImage->SetSpacing( sp );  
  std::cout << "Spacing = ";
  std::cout << sp[0] << ", " << sp[1] << ", " << sp[2] << std::endl;

  // scale image to meters
  const ImageType::PointType& orgn_mm = inputImage->GetOrigin();
  ImageType::PointType orgn;
  orgn[0] = 0.001 * orgn_mm[0];
  orgn[1] = 0.001 * orgn_mm[1];
  orgn[2] = 0.001 * orgn_mm[2];
  inputImage->SetOrigin(orgn);
  std::cout << "Origin = ";
  std::cout << orgn[0] << ", " << orgn[1] << ", " << orgn[2] << std::endl;


  printf("xpixel = %d  \nypixel = %d \nnslice = %d \n",dim1,dim2,dim3);
  printf("X0 = %f  \nDX = %f \n",orgn[0],sp[0]);
  printf("Y0 = %f  \nDY = %f \n",orgn[1],sp[1]);
  printf("Z0 = %f  \nDZ = %f \n",orgn[2],sp[2]);
  // if control file name given, save dimension info to EXISTING control file
  if(strcasecmp("none",inifile)){ // true if inifile != "none"
     CSimpleIni  Ini; // instantiate ini file
     bool bMulti;
     char imagedata[64]; // character buffer to write out solution
     SI_Error rc = Ini.LoadFile(inifile);
     if (rc < 0) { printf("Failed to open %s\n",inifile); abort(); }

     // write out pixel info of planning images
     //  x,y pixel info is the same for the mrti images
     sprintf(imagedata,"%d",dim1); 
     Ini.SetValue("planning_images","xpixel",imagedata);
     sprintf(imagedata,"%d",dim2); 
     Ini.SetValue("planning_images","ypixel",imagedata);
     sprintf(imagedata,"%d",dim3); 
     Ini.SetValue("planning_images","nslice",imagedata);

     // write out dimensions of planning imaging data in METERS!!!!!!!!!!!!!!!!!
     //  in planes images for planning and mrti images are the same

     sprintf(imagedata,"%e",orgn[0]);
     Ini.SetValue("planning_images","X0",imagedata);
     Ini.SetValue(     "mrti"      ,"X0",imagedata);

     sprintf(imagedata,"%e",orgn[1]);
     Ini.SetValue("planning_images","Y0",imagedata);
     Ini.SetValue(     "mrti"      ,"Y0",imagedata);

     // account for possibly different spacing of planning and thermal images 
     sprintf(imagedata,"%e",sp[0])  ;
     Ini.SetValue("planning_images","DX",imagedata);
     sprintf(imagedata,"%e",sp[0]*double(dim1) / 
                            atof(Ini.GetValue("mrti","xpixel","256",&bMulti)) );
     Ini.SetValue(     "mrti"      ,"DX",imagedata);

     sprintf(imagedata,"%e",sp[1])  ;
     Ini.SetValue("planning_images","DY",imagedata);
     sprintf(imagedata,"%e",sp[1]*double(dim2) / 
                            atof(Ini.GetValue("mrti","ypixel","256",&bMulti)) );
     Ini.SetValue(     "mrti"      ,"DY",imagedata);

     //  out of plane initial position is DIFFERENT 
     sprintf(imagedata,"%e",orgn[2]);
     Ini.SetValue("planning_images","Z0",imagedata);
     int mrti_image_zero = atoi(
             Ini.GetValue("planning_images","mrti_image_zero","-1",&bMulti)  ) ;
     sprintf(imagedata,"%e",orgn[2] + sp[2]*mrti_image_zero);
     Ini.SetValue(     "mrti"      ,"Z0",imagedata);
     boundcheck("thermal image zero plane",mrti_image_zero,dim3); //error check

     //  out of plane spacing is the same for mrti and planning
     sprintf(imagedata,"%e",sp[2])  ;
     Ini.SetValue("planning_images","DZ",imagedata);
     Ini.SetValue(     "mrti"      ,"DZ",imagedata);

     // as a backup plan write out the image bounds to be used 
     // to create a structured mesh
     sprintf(imagedata,"bx=[0,%d];by=[0,%d];bz=[0,%d]",dim1,dim2,dim3);
     Ini.SetValue("planning_images","boundingbox",imagedata);

     // use planning images to make an educated guess of the laser tip location
     int laser_xpixel =
         atoi(Ini.GetValue("planning_images","laser_xpixel"     ,"-1",&bMulti));
     boundcheck("laser xpixel",laser_xpixel,dim1); //error check
     int laser_ypixel =
         atoi(Ini.GetValue("planning_images","laser_ypixel"     ,"-1",&bMulti));
     boundcheck("laser ypixel",laser_ypixel,dim2); //error check
     int laser_slice_plane =
         atoi(Ini.GetValue("planning_images","laser_slice_plane","-1",&bMulti));
     boundcheck("laser slice plane ",laser_slice_plane,dim3); //error check

     sprintf(imagedata,"%e", orgn[0] + sp[0] * laser_xpixel      );
     Ini.SetValue("source_laser","X_0",imagedata);
     sprintf(imagedata,"%e", orgn[1] + sp[1] * laser_ypixel      );
     Ini.SetValue("source_laser","Y_0",imagedata);
     sprintf(imagedata,"%e", orgn[2] + sp[2] * laser_slice_plane );
     Ini.SetValue("source_laser","Z_0",imagedata);

     // write out small ucd file for vis. of laser position
     sprintf(file_o,"%slaser_pos.inp",dicomdir); 
     ofstream avsucd; //output file
     avsucd.open (file_o, ios::out);
     avsucd << "1 0 0 0 0" << endl;
     avsucd << "1 " ;
     avsucd << orgn[0] + sp[0] * laser_xpixel;
     avsucd << " " ;
     avsucd << orgn[1] + sp[1] * laser_ypixel;
     avsucd << " " ;
     avsucd << orgn[2] + sp[2] * laser_slice_plane;
     avsucd << endl;
     avsucd.close();


     // save
     Ini.SaveFile(inifile);
  }

// Software Guide : EndCodeSnippet

// Software Guide : BeginLatex
//
// Next we create the necessary iterators.  The const slice iterator walks
// the 3D input image. 
//
// Software Guide : EndLatex

// Software Guide : BeginCodeSnippet
  typedef itk::ImageSliceConstIteratorWithIndex< ImageType >  SliceConstIteratorType;
  typedef itk::ImageSliceIteratorWithIndex< ImageType >       SliceIteratorType;

  SliceConstIteratorType  inputIt(  inputImage,  inputImage->GetRequestedRegion() );
  SliceIteratorType       outputIt(  inputImage,  inputImage->GetRequestedRegion() );

/*         ________
 begin    |->|->|->|
           --------
          |->|->|->|        this seems to be the default 
   y       --------              direction of iteration of ITK
    ^     |->|->|->|             slice iterator w/
    |      --------              inputIt.SetFirstDirection(  0 );
    |     |->|->|->|             inputIt.SetSecondDirection( 1 );
    |      --------
    |     |->|->|->|   end 
          |--------|              ----------> x
*/
  inputIt.SetFirstDirection(  0 );
  inputIt.SetSecondDirection( 1 );
  outputIt.SetFirstDirection(  0 );
  outputIt.SetSecondDirection( 1 );
  
// At the end of the first slice, the input iterator steps to the first
// row in the next slice.   The process repeats
// until the last slice of the input is processed.
//
// Software Guide : EndLatex

// Software Guide : BeginCodeSnippet

  //  Put Data into CImg Data structure

  
  int iii=0,jjj=0,kkk=0; 
  inputIt.GoToBegin();
  
  while( !inputIt.IsAtEnd() )
    { jjj=0;
    while ( !inputIt.IsAtEndOfSlice() )
      { iii=0;
      while ( !inputIt.IsAtEndOfLine() )
        {
/* massage indices to match avs expected direction of iteration
           ________
          |->|->|->|   end 
           --------
          |->|->|->|        
  y        --------         
   ^      |->|->|->|        
   |       --------         
   |      |->|->|->|        
   |       --------
   | begin|->|->|->|  
          |--------|   ----------> x  
                                                      */
        //Filter(iii,dim2-(jjj+1),kkk) = inputIt.Get() ;
        Filter(iii,jjj,kkk) = inputIt.Get() ;
        iii++;
        ++inputIt;
        }
      inputIt.NextLine(); jjj++;

      }
    inputIt.NextSlice(); kkk++;
    }


  // apply Contrast Enhancement Filter

  const CImgStats stats(Filter);     // Compute basic statistics on the image.
  stats.print("Image statistics");   // Display statistics.


  if(image_enhance){
     // Enhancement algorithm only works with a normed image for some reason?
     Filter.normalize(0,255);
     const CImgStats normstats(Filter); // Compute statistics of normed image.
     normstats.print("Normed Image statistics"); // Display statistics.
     for(int i=0;i<dim3;i++){
        ImageEnhance( Filter.ptr(0,0,i) , Filter.width  , Filter.height );
     }
  }


  
  std::vector< float >    coordata;
  coordata.reserve(dim1+dim2+dim3);
  coordata.resize(dim1+dim2+dim3,0.0);

  for(int i=0;i<dim1;i++) coordata.at(          i) = orgn[0] + float(i)*sp[0];
  for(int i=0;i<dim2;i++) coordata.at(dim1     +i) = orgn[1] + float(i)*sp[1];
  for(int i=0;i<dim3;i++) coordata.at(dim1+dim2+i) = orgn[2] + float(i)*sp[2];

  //for(int i=0;i<dim1;i++) coordata.at(          i) =  float(i)*sp[0];
  //for(int i=0;i<dim2;i++) coordata.at(dim1     +i) =  float(i)*sp[1];
  //for(int i=0;i<dim3;i++) coordata.at(dim1+dim2+i) =  float(i)*sp[2];

  // units given in millimeters convert to meters

  for(int i=0;i<coordata.size();i++) coordata[i] =coordata[i];


  char byteid[6];
  if(endianess){ sprintf(byteid,"litend");  //native data type is little endian 
  } else {       sprintf(byteid,"bigend"); }//native data type is big    endian 
  sprintf(file_o,"%s%s%04d.fld",dicomdir,byteid,0); // for animate filename
  WriteAVS(file_o,&Filter,&coordata,false);
  // write out a byte swapped file as well
  if(endianess){ sprintf(byteid,"bigend");  //native data type is little endian 
  } else {       sprintf(byteid,"litend"); }//native data type is big    endian 
  sprintf(file_o,"%s%s%04d.fld",dicomdir,byteid,0); // for animate filename
  WriteAVS(file_o,&Filter,&coordata,true);

  //  write rawiv file


  float bounds[6],spans[3], rawivorigin[3];
  int  dims[3],numVerts,numCells;
  ofstream rawivfile;
  sprintf(file_o,"%s.rawiv",dicomdir);

  bounds[0]=coordata.at(0)        ; bounds[3]=coordata.at(dim1          -1);
  bounds[1]=coordata.at(dim1)     ; bounds[4]=coordata.at(dim1+dim2     -1);
  bounds[2]=coordata.at(dim1+dim2); bounds[5]=coordata.at(dim1+dim2+dim3-1);
  numVerts = (Filter.width  )*(Filter.height  )*(Filter.depth  );
  numCells = (Filter.width-1)*(Filter.height-1)*(Filter.depth-1);
  dims[0]=Filter.width; dims[1]=Filter.height; dims[2]=Filter.depth;
  rawivorigin[0]=bounds[0];
  rawivorigin[1]=bounds[1];
  rawivorigin[2]=bounds[2];
  spans[0]=(bounds[3]-bounds[0])/(dims[0]-1);
  spans[1]=(bounds[4]-bounds[1])/(dims[1]-1);
  spans[2]=(bounds[5]-bounds[2])/(dims[2]-1);
  if(endianess) swapByteOrder(bounds,(unsigned int)6);
    rawivfile.open(file_o, ios::out|ios::binary);
    rawivfile.write(reinterpret_cast<char *>(bounds), sizeof(float)*6);
    rawivfile.close();
  if(endianess) swapByteOrder(&numVerts,(unsigned int)1);
    rawivfile.open(file_o, ios::out|ios::binary|ios::app);
    rawivfile.write(reinterpret_cast<char *>(&numVerts), sizeof(int)*1);
    rawivfile.close();
  if(endianess) swapByteOrder(&numCells,(unsigned int)1);
    rawivfile.open(file_o, ios::out|ios::binary|ios::app);
    rawivfile.write(reinterpret_cast<char *>(&numCells), sizeof(int)*1);
    rawivfile.close();
  if(endianess) swapByteOrder(dims,(unsigned int)3);
    rawivfile.open(file_o, ios::out|ios::binary|ios::app);
    rawivfile.write(reinterpret_cast<char *>(dims), sizeof(int)*3);
    rawivfile.close();
  if(endianess) swapByteOrder(rawivorigin,(unsigned int)3);
    rawivfile.open(file_o, ios::out|ios::binary|ios::app);
    rawivfile.write(reinterpret_cast<char *>(rawivorigin), sizeof(float)*3);
    rawivfile.close();
  if(endianess) swapByteOrder(spans,(unsigned int)3);
    rawivfile.open(file_o, ios::out|ios::binary|ios::app);
    rawivfile.write(reinterpret_cast<char *>(spans), sizeof(float)*3);
    rawivfile.close();
  if(endianess) swapByteOrder(Filter.data, Filter.size());
    rawivfile.open(file_o, ios::out|ios::binary|ios::app);
    rawivfile.write(reinterpret_cast<char *>(Filter.data), 
                                      sizeof(int)*Filter.size());
    rawivfile.close();
  printf("%s written \n",file_o);
  // swap back
  if(endianess) swapByteOrder(Filter.data, Filter.size());

// Software Guide : BeginLatex
// 
// We proceed now to save the volumetric image in another file, as specified by
// the user in the command line arguments of this program. Thanks to the
// ImageIO factory mechanism, only the filename extension is needed to identify
// the file format in this case.
//
// Software Guide : EndLatex

// Software Guide : BeginCodeSnippet    
    typedef itk::ImageFileWriter< ImageType > WriterType;
    WriterType::Pointer writer = WriterType::New();
    
    sprintf(file_o,"%s.mha",dicomdir);

    writer->SetFileName( file_o );

    writer->SetInput( reader->GetOutput() );
// Software Guide : EndCodeSnippet    

    std::cout  << "Writing the image as " << std::endl << std::endl;
    std::cout  << file_o << std::endl << std::endl;


// Software Guide : BeginLatex
// 
// The process of writing the image is initiated by invoking the
// \code{Update()} method of the writer.
//
// Software Guide : EndLatex

    try
      {
// Software Guide : BeginCodeSnippet    
      writer->Update();
// Software Guide : EndCodeSnippet    
      }
    catch (itk::ExceptionObject &ex)
      {
      std::cout << ex << std::endl;
      return EXIT_FAILURE;
      }
    }
  catch (itk::ExceptionObject &ex)
    {
    std::cout << ex << std::endl;
    return EXIT_FAILURE;
    }



// Software Guide : BeginLatex
// 
// Note that in addition to writing the volumetric image to a file we could
// have used it as the input for any 3D processing pipeline. Keep in mind that
// DICOM is simply a file format and a network protocol. Once the image data
// has been loaded into memory, it behaves as any other volumetric dataset that
// you could have loaded from any other file format. 
//
// Software Guide : EndLatex


  return EXIT_SUCCESS;

}


void swapByteOrder(float* input, unsigned int num)
{
	unsigned int c;
	for (c=0; c<num; c++) {
		swapByteOrder(input[c]);
	}
}

void swapByteOrder(double* input, unsigned int num)
{
	unsigned int c;
	for (c=0; c<num; c++) {
		swapByteOrder(input[c]);
	}
}

void swapByteOrder(unsigned int* input, unsigned int num)
{
	unsigned int c;
	for (c=0; c<num; c++) {
		swapByteOrder(input[c]);
	}
}

void swapByteOrder(int* input, unsigned int num)
{
	unsigned int c;
	for (c=0; c<num; c++) {
		swapByteOrder(input[c]);
	}
}

void swapByteOrder(unsigned short* input, unsigned int num)
{
	unsigned int c;
	for (c=0; c<num; c++) {
		swapByteOrder(input[c]);
	}
}

void WriteAVS(char *Filename,CImg<PixelType> *Filter,
                       std::vector< float > *Coordata, bool ByteSwap)
{
  ofstream avsfile;
  avsfile.open(Filename, ios::out);
  // write out avs file header
  avsfile <<"# AVS field file  " << endl;
  avsfile <<"ndim=3          # number of dimensions in the field" << endl;
  avsfile <<"dim1="<<Filter->width  <<"        # dimension of axis 1 " << endl;
  avsfile <<"dim2="<<Filter->height <<"        # dimension of axis 2 " << endl;
  avsfile <<"dim3="<<Filter->depth  <<"        # dimension of axis 3 " << endl;
  avsfile <<"nspace=3        # number of physical coordinates per point " << endl;
  avsfile <<"veclen=1        # number of components at each point " << endl;
  avsfile <<"data=float      # native format of linux"<<endl;
  avsfile <<"field=rectilinear # fld type(uniform,rectilinear,irregular)"<<endl;
  avsfile <<'\f'<<'\f';
  avsfile.close();
  // write out pixel data
  if(ByteSwap) swapByteOrder(Filter->data, Filter->size()); 
  avsfile.open (Filename, ios::out|ios::binary|ios::app);
  avsfile.write (reinterpret_cast<char *>(Filter->data),
                                         sizeof(float)*Filter->size());
  avsfile.close();
  if(ByteSwap) swapByteOrder(Filter->data, Filter->size()); // swap back
  // write out coordinate data
  if(ByteSwap) swapByteOrder(&(*Coordata)[0], Coordata->size()); 
  avsfile.open (Filename, ios::out|ios::binary|ios::app);
  avsfile.write (reinterpret_cast<char *>(&(*Coordata)[0]), 
                                              sizeof(float)*Coordata->size());
  avsfile.close();
  if(ByteSwap) swapByteOrder(&(*Coordata)[0], Coordata->size()); // swap back
  printf("%s written \n",Filename);
}



void ImageEnhance(PixelType *dataset,int XDIM,int YDIM)
{
  int i,j;
  float *lcmax, *lcmin;
  float *tmpmax, *tmpmin;
  float lmin, lmax, img, avg;
  float tempmax,tempmin,temp;
  float window;
  float a,b,c,alpha;

  float *imgavg;
  


   /* Enhancing the image */
   lcmin = (float*)malloc(sizeof(float)*XDIM*YDIM);
   lcmax = (float*)malloc(sizeof(float)*XDIM*YDIM);
   tmpmin = (float*)malloc(sizeof(float)*XDIM*YDIM);
   tmpmax = (float*)malloc(sizeof(float)*XDIM*YDIM);
   imgavg = (float*)malloc(sizeof(float)*XDIM*YDIM);

   for (j=0; j<YDIM; j++)
     for (i=0; i<XDIM; i++) {
       lcmin[j*XDIM+i] = dataset[j*XDIM+i];
       lcmax[j*XDIM+i] = dataset[j*XDIM+i];
       tmpmin[j*XDIM+i] = dataset[j*XDIM+i];
       tmpmax[j*XDIM+i] = dataset[j*XDIM+i];
       imgavg[j*XDIM+i] = dataset[j*XDIM+i];
     }
   

   
   for (i=1; i<XDIM; i++) {
     imgavg[i] += RESISTOR*(imgavg[i-1]-imgavg[i]);
     if (tmpmin[i-1] < tmpmin[i])
       tmpmin[i] += RESISTOR*(tmpmin[i-1]-tmpmin[i]);
     if (tmpmax[i-1] > tmpmax[i])
       tmpmax[i] += RESISTOR*(tmpmax[i-1]-tmpmax[i]);
   }
   
   for (i=XDIM-2; i>=0; i--) {
     imgavg[i] += RESISTOR*(imgavg[i+1]-imgavg[i]);
     if (tmpmin[i+1] < tmpmin[i])
       tmpmin[i] += RESISTOR*(tmpmin[i+1]-tmpmin[i]);
     if (tmpmax[i+1] > tmpmax[i])
       tmpmax[i] += RESISTOR*(tmpmax[i+1]-tmpmax[i]);
   }


   for (j=1; j<YDIM; j++) {

     for (i=0; i<XDIM; i++) {
       imgavg[j*XDIM+i] += RESISTOR*(imgavg[(j-1)*XDIM+i]-imgavg[j*XDIM+i]);
       if (tmpmin[(j-1)*XDIM+i] < tmpmin[j*XDIM+i])
	 tmpmin[j*XDIM+i] += RESISTOR*(tmpmin[(j-1)*XDIM+i]-tmpmin[j*XDIM+i]);
       if (tmpmax[(j-1)*XDIM+i] > tmpmax[j*XDIM+i])
	 tmpmax[j*XDIM+i] += RESISTOR*(tmpmax[(j-1)*XDIM+i]-tmpmax[j*XDIM+i]);
     }

     for (i=1; i<XDIM; i++) {
       imgavg[j*XDIM+i] += RESISTOR*(imgavg[j*XDIM+i-1]-imgavg[j*XDIM+i]);
       if (tmpmin[j*XDIM+i-1] < tmpmin[j*XDIM+i])
	 tmpmin[j*XDIM+i] += RESISTOR*(tmpmin[j*XDIM+i-1]-tmpmin[j*XDIM+i]);
       if (tmpmax[j*XDIM+i-1] > tmpmax[j*XDIM+i])
	 tmpmax[j*XDIM+i] += RESISTOR*(tmpmax[j*XDIM+i-1]-tmpmax[j*XDIM+i]);
     }
    
     for (i=XDIM-2; i>=0; i--) {
       imgavg[j*XDIM+i] += RESISTOR*(imgavg[j*XDIM+i+1]-imgavg[j*XDIM+i]);
       if (tmpmin[j*XDIM+i+1] < tmpmin[j*XDIM+i])
	 tmpmin[j*XDIM+i] += RESISTOR*(tmpmin[j*XDIM+i+1]-tmpmin[j*XDIM+i]);
       if (tmpmax[j*XDIM+i+1] > tmpmax[j*XDIM+i])
	 tmpmax[j*XDIM+i] += RESISTOR*(tmpmax[j*XDIM+i+1]-tmpmax[j*XDIM+i]);
     }

   }
   

   j=YDIM-1;
   
   for (i=XDIM-2; i>=0; i--) {
     imgavg[j*XDIM+i] += RESISTOR*(imgavg[j*XDIM+i+1]-imgavg[j*XDIM+i]);
     if (lcmin[j*XDIM+i+1] < lcmin[j*XDIM+i])
       lcmin[j*XDIM+i] += RESISTOR*(lcmin[j*XDIM+i+1]-lcmin[j*XDIM+i]);
     if (lcmax[j*XDIM+i+1] > lcmax[j*XDIM+i])
       lcmax[j*XDIM+i] += RESISTOR*(lcmax[j*XDIM+i+1]-lcmax[j*XDIM+i]);
   }

   for (i=1; i<XDIM; i++) {
     imgavg[j*XDIM+i] += RESISTOR*(imgavg[j*XDIM+i-1]-imgavg[j*XDIM+i]);
     if (lcmin[j*XDIM+i-1] < lcmin[j*XDIM+i])
       lcmin[j*XDIM+i] += RESISTOR*(lcmin[j*XDIM+i-1]-lcmin[j*XDIM+i]);
     if (lcmax[j*XDIM+i-1] > lcmax[j*XDIM+i])
       lcmax[j*XDIM+i] += RESISTOR*(lcmax[j*XDIM+i-1]-lcmax[j*XDIM+i]);
   }

    
   for (j=YDIM-2; j>=0; j--) {

     for (i=0; i<XDIM; i++) {
       imgavg[j*XDIM+i] += RESISTOR*(imgavg[(j+1)*XDIM+i]-imgavg[j*XDIM+i]);
       if (lcmin[(j+1)*XDIM+i] < lcmin[j*XDIM+i])
	 lcmin[j*XDIM+i] += RESISTOR*(lcmin[(j+1)*XDIM+i]-lcmin[j*XDIM+i]);
       if (lcmax[(j+1)*XDIM+i] > lcmax[j*XDIM+i])
	 lcmax[j*XDIM+i] += RESISTOR*(lcmax[(j+1)*XDIM+i]-lcmax[j*XDIM+i]);
     }
     
     for (i=XDIM-2; i>=0; i--) {
       imgavg[j*XDIM+i] += RESISTOR*(imgavg[j*XDIM+i+1]-imgavg[j*XDIM+i]);
       if (lcmin[j*XDIM+i+1] < lcmin[j*XDIM+i])
	 lcmin[j*XDIM+i] += RESISTOR*(lcmin[j*XDIM+i+1]-lcmin[j*XDIM+i]);
       if (lcmax[j*XDIM+i+1] > lcmax[j*XDIM+i])
	 lcmax[j*XDIM+i] += RESISTOR*(lcmax[j*XDIM+i+1]-lcmax[j*XDIM+i]);
     }
	
     for (i=1; i<XDIM; i++) {
       imgavg[j*XDIM+i] += RESISTOR*(imgavg[j*XDIM+i-1]-imgavg[j*XDIM+i]);
       if (lcmin[j*XDIM+i-1] < lcmin[j*XDIM+i])
	 lcmin[j*XDIM+i] += RESISTOR*(lcmin[j*XDIM+i-1]-lcmin[j*XDIM+i]);
       if (lcmax[j*XDIM+i-1] > lcmax[j*XDIM+i])
	 lcmax[j*XDIM+i] += RESISTOR*(lcmax[j*XDIM+i-1]-lcmax[j*XDIM+i]);
     }

   }
   
   

   tempmin = 9999;
   tempmax = -9999;
   for (j=0; j<YDIM; j++)
     for (i=0; i<XDIM; i++) {
       
       lmin = min(lcmin[j*XDIM+i], tmpmin[j*XDIM+i]);
       lmax = max(lcmax[j*XDIM+i], tmpmax[j*XDIM+i]);
       img = dataset[j*XDIM+i];
       avg = imgavg[j*XDIM+i];
       
       window = lmax - lmin;
       window = sqrt(window*(510-window));
       if (lmin != lmax) {
	 img = window*(img-lmin)/(lmax-lmin);
	 avg = window*(avg-lmin)/(lmax-lmin);
       }

       alpha = (avg-img)/(181.019 * window);
	                 
       if (alpha != 0) {

	 a = 0.707 * alpha;
	 b = 1.414*alpha*(img - window) - 1;
	 c = 0.707*alpha*img*(img-2*window) + img;
	 
	 lcmin[j*XDIM+i] = lmin+(-b-sqrt(b*b-4*a*c))/(2*a);
	 
       }
       else {

	 lcmin[j*XDIM+i] = lmin+img;

       }
       

       temp = lcmin[j*XDIM+i];
       if (temp > tempmax)
	 tempmax = temp;
       if (temp < tempmin)
	 tempmin = temp;
       
     }
   
   printf("min = %f   max = %f \n",tempmin,tempmax);
   
   /* the following "0.0" could be set to other value (e.g., 0.1)
      in case you want to sharpen the enhanced images */
   temp = tempmax-tempmin;
   tempmin += 0.0*temp;
   tempmax -= 0.0*temp;
   
   for (j=0; j<YDIM; j++)
     for (i=0; i<XDIM; i++) {
       
       if (lcmin[j*XDIM+i] < tempmin)
	 lcmin[j*XDIM+i] = tempmin;
       if (lcmin[j*XDIM+i] > tempmax)
	 lcmin[j*XDIM+i] = tempmax;
       
       dataset[j*XDIM+i] = (PixelType)((lcmin[j*XDIM+i]-tempmin)
			   *255.0/(tempmax-tempmin));
     }
   

   free(lcmin);
   free(tmpmin);
   free(lcmax);
   free(tmpmax);
   free(imgavg);

}
