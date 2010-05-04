import sys
import os
import ConfigParser # file parser
import itkUtilities # custom python module for access to itk library

# use python Set: union, intersection, difference operations
# to deal with real-time files in the file system
#http://docs.python.org/library/sets.html
from sets import Set

class AcquisitionBaseClass:
  def __init__(self,ControlFile):
    #read the main control file
    self.iniFile  = ConfigParser.ConfigParser()
    self.iniFile.readfp( open("%s" % ControlFile , "r") )
    self.instanceTag = "0020|0013"
    #get dicom Dictionary
    try:
      self.dictionary = self.iniFile.get(   "mrti" ,"dictionary"  )
    except ConfigParser.NoOptionError:
      raise IOError("\n\n    No Dicom Dictionary FOUND! " )
    # get exampath info
    self.ExamPath = self.iniFile.get(   "mrti" ,"exampath"    )
    # Set C++ executable to get header data info
    try:
      self.DicomToComplex = self.iniFile.get(   "compexec" ,"Dicomtocomplex"  )
    except ConfigParser.NoOptionError: 
      self.DicomToComplex = "~/DDDAS/trunk/utilities/DicomReadWriteComplexImage"
  def GetHeaderInfo(self):
    # loop over possible slice number tags to find the slice location
    for sliceIDKey in ["0021|104f"]:
      nslice = int( itkUtilities.GetDicomTag( dicomMagFile, sliceIDKey ,
                                                            self.dictionary) )


class SiemensAquisitionGEWrite(AcquisitionBaseClass):
  # try building list of magnitude and phase data
  magnID   = iniFile.getint("mrti" ,"magnitudeid" )
  phasID   = iniFile.getint("mrti" ,"phaseid"     )
  fileListOne = Set(os.listdir( "%s/s%d" % (ExamPath, magnID ) ))
  numFiles = len(fileListOne) 
  fileListTwo = os.listdir( "%s/s%d" % (ExamPath,phasID) ) 
  try:
     assert numFiles == len(fileListTwo)     
  except AssertionError:
     raise IndexError("\n\n    error with number of files in the directory " )
  # set output to magnitude id
  dirID = magnID
  # TODO add realtime abilities as data is written to the ExamPath
  for magnitudeFile,phaseFile in zip(fileListOne,fileListTwo): 
      # copy magnitude first
      dicomMagFile   = "%s/s%d/%s" % (ExamPath, magnID ,magnitudeFile) 
      magTimeID = int( itkUtilities.GetDicomTag( dicomMagFile, instanceTag,
                                                                dictionary) )
      os.system( "cp %s %s/Processed/s%d/i%d.MRDC.%d" % (dicomMagFile,ExamPath,dirID,dirID+2*magTimeID-1,2*magTimeID-1) )
      # now copy phase
      dicomPhaseFile = "%s/s%d/%s" % (ExamPath,phasID,  phaseFile  ) 
      phaseTimeID = int( itkUtilities.GetDicomTag( dicomPhaseFile, instanceTag,
                                                                   dictionary) )
      os.system( "cp %s %s/Processed/s%d/i%d.MRDC.%d" % (dicomPhaseFile,ExamPath,dirID,dirID+2*phaseTimeID,2*phaseTimeID) )
      print  "%s %d %s %d " % (magnitudeFile,magTimeID,phaseFile,phaseTimeID) 
   
   
   # one based
   #for fileid in range(1,numFiles+1):
   #    os.system( "%s -M%s/Processed/s%d/i%d.MRDC.%d -P%s/Processed/s%d/i%d.MRDC.%d --output %s/Processed/s%d/image -timeid %d" % ( 
   #                  DicomToComplex ,
   #                      ExamPath,dirID,dirID+2*fileid-1,2*fileid-1, 
   #                      ExamPath,dirID,dirID+2*fileid  ,2*fileid  , 
   #                      ExamPath,dirID,fileid-1 )  ) 

class GEAquisitionGEWrite(AcquisitionBaseClass):
  dirID    = iniFile.getint("mrti" ,"dirid"     )
  # try building list of real and imaginary data as default
  realID   = iniFile.getint("mrti" ,"realid"      )
  imagID   = iniFile.getint("mrti" ,"imaginaryid" )
  try:
     assert realID == imagID 
  except AssertionError:
     raise IndexError("\n\n    i have never encountered real and imaginary data in different directories... " )

#single file use
if __name__ == "__main__":
   #ini file passed in from the command line
   #controlfile = "/mnt/data/scratch/fuentes/DDDAS/data/biotex/090318_751642_treat/controlBrain.ini"
   controlfile = sys.argv[1] 
   
   #initialize
   try: 
     # try Siemens Acquisition w/ GE Header
     fileProcess = SiemensAquisitionGEWrite(controlfile) 
   except ConfigParser.NoOptionError:
     # try GE Acquisition w/ GE Header
     fileProcess = GEAquisitionGEWrite(controlfile) 
   
   # get some initial header info from the first file...
   fileProcess.GetHeaderInfo()
   
   # make new directory
   os.system('mkdir -p %s/Processed/s%d' % (fileProcess.ExamPath,
                                            fileProcess.dirID) )
   
   # typically assume a high number of images to read in 
   ntime  = fileProcess.iniFile.getint("mrti","ntime")
   timeID = 0 
   while (timeID < ntime): 
      fileProcess.updateDirectoryInfo()
      # parse the new files found, convert to complex format,
      # and update the timeID
      fileProcess.parseAndWriteAsComplexData()
   
