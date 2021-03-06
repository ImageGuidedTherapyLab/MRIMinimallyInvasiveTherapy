import sys
import os
import ConfigParser # file parser
import itkUtilities # custom python module for access to itk library
# TODO use python Set: union, intersection, difference operations
# to deal with real-time files in the file system
#http://docs.python.org/library/sets.html
from sets import Set

thisScript = sys.argv[0]
class AcquisitionBaseClass:
  """ Base Class for realtime image conversion...  """
  def __init__(self,ControlFile):
    print " base class constructor called \n\n" 
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
      self.DicomToComplex = thisScript[:thisScript.rfind("/")] +"/DicomReadWriteComplexImage"
  def GetValueFromKeyList(self,Section,Entry,fileName,KeyList):
    try : # first check if the user input the value
      Value =  self.iniFile.get( Section , Entry )
      print "Read %s = %s from iniFile " % (Entry,Value) 
      return Value
    # if option not found attempt to get slice info from header
    except ConfigParser.NoOptionError:
      # loop over possible tags to find the key of interest
      KeyList.append(None)
      # KeyList should be of the form ["xxxx|xxxx","xxxx|xxxx",...,None] 
      for potentialIDKey in KeyList:
        if (potentialIDKey == None): #Last Key has been reached
          raise RuntimeError("\n\n    Can't find %s Tag! " % Entry )
        try : 
          return itkUtilities.GetDicomTag(fileName, potentialIDKey ,
                                                    self.dictionary) 
        except ValueError: 
          pass # try the next key
  def GetHeaderInfo(self,fileName):
    # get file size of "representative file
    self.filesize = os.path.getsize(fileName) 
    # number of slice locations
    self.nslice   = int( self.GetValueFromKeyList("mrti","nslice",
                                                     fileName, ["0021|104f"] ) )
    # echo time
    self.echoTime = float( self.GetValueFromKeyList("mrti","echotime",
                                                     fileName, ["0018|0081"] ) )
    # imaging frequency / gyromagnetic ratio
    self.imagFreq = float( self.GetValueFromKeyList("mrti","imagfreq",
                                                     fileName, ["0018|0084"] ) )
    # number of echos 
    self.necho = int(   self.GetValueFromKeyList("mrti","necho",
                                                     fileName, ["0019|107e"] ) )
    # get image acquisition time
    try : # first check if the user input the value
       self.deltat = float(self.iniFile.get( "mrti" , "deltat" ))
    # if option not found attempt to get slice info from header
    except ConfigParser.NoOptionError:
       acqDurationIdkey = "0019|105a"
       fastPhaseIdkey   = "0019|10f2" 
       acqDuration = float( itkUtilities.GetDicomTag(fileName, 
                                         acqDurationIdkey , self.dictionary) )
       fastPhase = int( itkUtilities.GetDicomTag(fileName, 
                                         fastPhaseIdkey   , self.dictionary) )
       # convert from \mu second to second
       self.acquisitionTime = acqDuration / fastPhase * 1.e-6;
    except ValueError: #TODO do i need to nest this
       print "estimating acquistion time as TR*number phase encodes..."
       repetitionTimeIdkey  = "0019|105a" #FIXME need the appropriate key
       numPhaseEncodesIdkey = "0019|10f2" #FIXME need the appropriate key
       repetitionTime = float( itkUtilities.GetDicomTag(fileName, 
                                         repetitionTimeIdkey, self.dictionary) )
       numPhaseEncodes = int( itkUtilities.GetDicomTag(fileName, 
                                      numPhaseEncodesIdkey  , self.dictionary) )
       self.deltat = numPhaseEncodesIdkey  * repetitionTime 
    HeaderIni = ConfigParser.ConfigParser()
    inifile = open("%s/Processed/s%d/imageHeader.ini" % \
                                             (self.ExamPath,self.dirID) , "w")
    HeaderIni.add_section("rawdata")
    HeaderIni.set("rawdata" ,"deltat"  ,"%f" % self.deltat  )
    HeaderIni.set("rawdata" ,"imagfreq","%f" % self.imagFreq)
    HeaderIni.set("rawdata" ,"echotime","%f" % self.echoTime)
    HeaderIni.set("rawdata" ,"necho"   ,"%d" % self.necho   )
    HeaderIni.write(inifile)
    inifile.close; inifile.flush() #ensure entire file written before continuing



class SiemensAquisitionGEWrite(AcquisitionBaseClass):
  """ This is used for the Siemens Aquisitions data that was transfered to a GE
      machine and written in GE DICOM format. 
      This is NOT intended for real-time...
  """
  def __init__(self,ControlFile):
    # call base class constructor to setup control file
    AcquisitionBaseClass.__init__(self,ControlFile)
    # try building list of magnitude and phase data
    magnID   = self.iniFile.getint("mrti" ,"magnitudeid" )
    phasID   = self.iniFile.getint("mrti" ,"phaseid"     )
    # this should not be done in real time
    # read list of all files
    fileMagnList = os.listdir( "%s/s%d" % (self.ExamPath,magnID) )
    numFiles = len(fileMagnList) 
    filePhasList = os.listdir( "%s/s%d" % (self.ExamPath,phasID) ) 
    try:
       assert numFiles == len(filePhasList)     
    except AssertionError:
       raise IndexError("\n\n    error with number of files in the directory " )
    # intialize file list
    self.fileList = {}
    # preprocess and sort by time instance
    # sorting by number of slices is done in  getFileList
    for magnitudeFile,phaseFile in zip(fileMagnList,filePhasList): 
      # store magnitude time id
      dicomMagnFile   = "%s/s%d/%s" % (self.ExamPath,magnID,magnitudeFile) 
      magnTimeID   = int( itkUtilities.GetDicomTag( dicomMagnFile, 
                                        self.instanceTag, self.dictionary) )
      try: # if key exists then append the file
         self.fileList["%d" % magnTimeID].append( ["-M",dicomMagnFile] )
      except KeyError: 
         self.fileList["%d" % magnTimeID]   =   [ ["-M",dicomMagnFile] ]
      # store phase time id
      dicomPhasFile = "%s/s%d/%s" % (self.ExamPath,phasID,    phaseFile) 
      phasTimeID = int( itkUtilities.GetDicomTag( dicomPhasFile, 
                                        self.instanceTag, self.dictionary) )
      try: 
         self.fileList["%d" % phasTimeID].append( ["-P",dicomPhasFile] )
      except KeyError: 
         self.fileList["%d" % phasTimeID]   =   [ ["-P",dicomPhasFile] ]
    # ensure time instance 0 exists
    if "0" not in self.fileList:
      self.fileList["0"] = self.fileList["1"] 
    # set output to magnitude id
    self.dirID = magnID
  def getFileList(self,istep):
      nextFileList = ""
      for dataType,fileName in self.fileList["%d" % istep]:
         # TODO sort by slice number, works for single slice for now...
         # sliceID = int( itkUtilities.GetDicomTag( fileName , sliceKey,
         #                                                     dictionary) )
         nextFileList = nextFileList + " %s%s " % (dataType,fileName)
      return nextFileList 
  def GetHeaderInfo(self):
      dataType,fileName  = (self.fileList["%d" % 1])[0]
      # call base class to get basic header info
      AcquisitionBaseClass.GetHeaderInfo(self,fileName)
   

class GEAquisitionGEWrite(AcquisitionBaseClass):
  """ This is used for the GE Aquisition written in GE DICOM format. 
      This IS intended for real-time...
  """
  def __init__(self,ControlFile):
    # call base class constructor to setup control file
    AcquisitionBaseClass.__init__(self,ControlFile)
    # get directory id...
    self.dirID    = iniFile.getint("mrti" ,"dirid"     )
    try:
       assert realID == imagID 
    except AssertionError:
       raise IndexError("\n\n    i have never encountered real and imaginary data in different directories... " )
    # TODO ensure can handle time instance zero
  def getFileList(self,istep):
      nextFileList = " "
      iecho = 0  # TODO add CSI  > 1 echo
      for jjj in range(self.nslice):
        realID = 2*nslice*(necho*istep+iecho+noffset) + (2*jjj+1)
        imagID = realID + 1
        realFile = "%s/s%d/i%d.MRDC.%d" % (ExamPath,DirId,DirId+realID,realID) 
        imagFile = "%s/s%d/i%d.MRDC.%d" % (ExamPath,DirId,DirId+imagID,imagID) 
        # check that file is available AND full file is available
        if( os.access(realFile,os.R_OK) and os.access(imagFile,os.R_OK) 
                                        and 
            os.path.getsize(realFile) >= self.filesize 
                                        and
            os.path.getsize(imagFile) >= self.filesize):
          nextFileList = nextFileList + " -R%s/s%d/i%d.MRDC.%d" % (
                                      ExamPath, DirId, DirId + realID, realID) 
          nextFileList = nextFileList + " -I%s/s%d/i%d.MRDC.%d" % (
                                      ExamPath, DirId, DirId + imagID, imagID) 
        else:
          raise IOError
  def GetHeaderInfo(self):
      fileID = 2*nslice*(necho*1+iecho+noffset) + 1
      fileName = "%s/s%d/i%d.MRDC.%d" % (ExamPath,DirId,DirId+fileID,fileID) 
      # call base class to get basic header info
      AcquisitionBaseClass.GetHeaderInfo(self,fileName)

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

   # make new directory
   os.system('mkdir -p %s/Processed/s%d' % (fileProcess.ExamPath,fileProcess.dirID) )
   
   # get some initial header info from the first file...
   fileProcess.GetHeaderInfo()
   
   timeID = 0  # codes expect time instances to start from 0
   while (timeID <= ntime): 
      try: # try to update the FileList and convert to complex format
         FileList = fileProcess.getFileList(timeID)
         # write the files
         os.system( "%s %s --output %s/Processed/s%d/image -timeid %d" % ( 
                fileProcess.DicomToComplex , FileList, fileProcess.ExamPath, 
                fileProcess.dirID, timeID ) ) 
         # update the timeID
         timeID = timeID + 1 
      # IOError should be thrown if files not found
      except IOError: 
         time.sleep(2)  # wait a few seconds and try again
