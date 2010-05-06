import sys
import os
import ConfigParser # file parser
import itkUtilities # custom python module for access to itk library

# TODO use python Set: union, intersection, difference operations
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
  def GetHeaderInfo(self,fileName):
    # get file size of "representative file
    self.filesize = os.path.getsize(fileName) 
    # loop over possible slice number tags to find the slice location
    for sliceIDKey in ["0021|104f"]:
      self.nslice = int( itkUtilities.GetDicomTag( fileName, sliceIDKey ,
                                                            self.dictionary) )

class SiemensAquisitionGEWrite(AcquisitionBaseClass):
"""
    This is used for the Siemens Aquisitions data that was transfered to a GE
machine and written in GE DICOM format. This is NOT intended for real-time...
"""
  def __init__(self):
    # try building list of magnitude and phase data
    magnID   = iniFile.getint("mrti" ,"magnitudeid" )
    phasID   = iniFile.getint("mrti" ,"phaseid"     )
    # this should not be done in real time
    # read list of all files
    fileMagnList = os.listdir( "%s/s%d" % (self.ExamPath,magnID) )
    numFiles = len(fileMagnList) 
    filePhasList = os.listdir( "%s/s%d" % (self.ExamPath,phasID) ) 
    try:
       assert numFiles == len(filePhasList)     
    except AssertionError:
       raise IndexError("\n\n    error with number of files in the directory " )
    # preprocess and sort by time instance
    # sorting by number of slices is done in  getFileList
    for magnitudeFile,phaseFile in zip(fileMagnList,filePhasList): 
      # store magnitude time id
      dicomMagnFile   = "%s/s%d/%s" % (self.ExamPath,magnID,magnitudeFile) 
      magnTimeID   = int( itkUtilities.GetDicomTag( dicomMagnFile, instanceTag,
                                                                   dictionary) )
      try: # if key exists then append the file
         self.fileList["%d" % magnTimeID].append( ["-M",dicomMagnFile] )
      except KeyError: 
         self.fileList["%d" % magnTimeID]   =   [ ["-M",dicomMagnFile] ]
      # store phase time id
      dicomPhasFile = "%s/s%d/%s" % (self.ExamPath,phasID,    phaseFile) 
      phasTimeID = int( itkUtilities.GetDicomTag( dicomPhasFile, instanceTag,
                                                                   dictionary) )
      try: 
         self.fileList["%d" % phasTimeID].append( ["-P",dicomPhasFile] )
      except KeyError: 
         self.fileList["%d" % phasTimeID]   =   [ ["-P",dicomPhasFile] ]
    # set output to magnitude id
    self.dirID = magnID
  def getFileList(self,istep):
      nextFileList = ""
      for dataType,fileName in self.fileList["%d" % istep]:
         # TODO sort by slice number, works for single slice for now...
         # sliceID = int( itkUtilities.GetDicomTag( fileName , sliceKey,
         #                                                     dictionary) )
         nextFileList = nextFileList + " %s %s " % (dataType,fileName)
      return nextFileList 
  def GetHeaderInfo(self):
      dataType,fileName  = self.fileList["%d" % 1]
      # call base class to get basic header info
      AcquisitionBaseClass.GetHeaderInfo(self,fileName):
   

class GEAquisitionGEWrite(AcquisitionBaseClass):
"""
    This is used for the GE Aquisition written in GE DICOM format. This IS intended for real-time...
"""
  def __init__(self):
    self.dirID    = iniFile.getint("mrti" ,"dirid"     )
    # try building list of real and imaginary data as default
    self.realID   = iniFile.getint("mrti" ,"realid"      )
    self.imagID   = iniFile.getint("mrti" ,"imaginaryid" )
    try:
       assert realID == imagID 
    except AssertionError:
       raise IndexError("\n\n    i have never encountered real and imaginary data in different directories... " )
  def getFileList(self,istep):
      nextFileList = " "
      iecho = 0  # TODO add CSI  > 1 echo
      for jjj in range(self.nslice):
        realID = 2*nslice*(necho*istep+iecho+noffset) + (2*jjj+1)
        imagID = realID + 1
        realFile = "%s/s%d/i%d.MRDC.%d" % (ExamPath,DirId,DirId+realID,realID) )
        imagFile = "%s/s%d/i%d.MRDC.%d" % (ExamPath,DirId,DirId+imagID,imagID) )
        # check that file is available AND full file is available
        if( os.access(realFile,os.R_OK) and os.access(imagFile,os.R_OK) 
                                        and 
            os.path.getsize(realFile) >= self.filesize 
                                        and
            os.path.getsize(imagFile) >= self.filesize):
          nextFileList = nextFileList + " -R%s/s%d/i%d.MRDC.%d" % (
                                      ExamPath, DirId, DirId + realID, realID) )
          nextFileList = nextFileList + " -I%s/s%d/i%d.MRDC.%d" % (
                                      ExamPath, DirId, DirId + imagID, imagID) )
        else:
          raise IOError
  def GetHeaderInfo(self):
      fileID = 2*nslice*(necho*1+iecho+noffset) + 1
      fileName = "%s/s%d/i%d.MRDC.%d" % (ExamPath,DirId,DirId+fileID,fileID) )
      # call base class to get basic header info
      AcquisitionBaseClass.GetHeaderInfo(self,fileName):

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
   
   # typically assume a high number of images to read in 
   ntime  = fileProcess.iniFile.getint("mrti","ntime")

   # get some initial header info from the first file...
   fileProcess.GetHeaderInfo()
   
   # make new directory
   os.system('mkdir -p Processed/s%d' % (fileProcess.dirID) )
   
   timeID = 0  # intialize and begin reading in files
   while (timeID < ntime): 
      try: # try to update the FileList and convert to complex format
         FileList = fileProcess.getFileList(timeID)
         # write the files
         os.system( "%s %s --output Processed/s%d/image -timeid %d" % ( 
                 DicomToComplex , FileList, fileProcess.dirID, timeID ) ) 
         # update the timeID
         timeID = timeID + 1 
      # IOError should be thrown if files not found
      except IOError: 
         time.sleep(2)  # wait a few seconds and try again
