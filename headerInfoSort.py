import sys
import os
import ConfigParser # file parser
# use python Set: union, intersection, difference operations
# to deal with real-time files in the file system
#http://docs.python.org/library/sets.html
from sets import Set

#the jobid is passed in from the command line
controlfile = sys.argv[1]
#controlfile = "/mnt/data/scratch/fuentes/DDDAS/data/biotex/090318_751642_treat/controlBrain.ini"

#read the main control file
iniFile  = ConfigParser.ConfigParser()
iniFile.readfp( open("%s" % controlfile , "r") )

#get dicom Dictionary
try:
  dictionary= iniFile.get(   "mrti" ,"dictionary"  )
except ConfigParser.NoOptionError:
  raise IOError("\n\n    No Dicom Dictionary FOUND! " )
instanceTag = "0020|0013"

# get exampath info
ExamPath   = iniFile.get(   "mrti" ,"exampath"    )

# Set C++ executable to get header data info
# TODO wrap c++ code to return  value directly
printTagExe = "/home/dfuentes/DDDAS/trunk/utilities/DicomImageReadReturnIntTag"
try:
  DicomToComplex = iniFile.get(   "compexec" ,"Dicomtocomplex"  )
except ConfigParser.NoOptionError: 
  DicomToComplex = "~/DDDAS/trunk/utilities/DicomReadWriteComplexImage"

# get some initial header info from the first file...
try: 
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
except ConfigParser.NoOptionError:
  dirID    = iniFile.getint("mrti" ,"dirid"     )
  # try building list of real and imaginary data as default
  realID   = iniFile.getint("mrti" ,"realid"      )
  imagID   = iniFile.getint("mrti" ,"imaginaryid" )
  try:
     assert realID == imagID 
  except AssertionError:
     raise IndexError("\n\n    i have never encountered real and imaginary data in different directories... " )

# make new directory
os.system('mkdir -p %s/Processed/s%d' % (ExamPath,dirID) )




# typically assume a high number of images to read in 
ntime  = iniFile.getint("mrti","ntime")
timeID = 0 
while (timeID < ntime): 


# TODO add realtime abilities as data is written to the ExamPath
for magnitudeFile,phaseFile in zip(fileListOne,fileListTwo): 
    # copy magnitude first
    dicomMagFile   = "%s/s%d/%s" % (ExamPath, magnID ,magnitudeFile) 
    fin,fout = os.popen4( "%s %s '%s' %s 2> /dev/null " % (printTagExe,dicomMagFile,instanceTag,dictionary) )
    result = fout.read()
    magTimeID = int(result.split("=")[1])
    os.system( "cp %s %s/Processed/s%d/i%d.MRDC.%d" % (dicomMagFile,ExamPath,dirID,dirID+2*magTimeID-1,2*magTimeID-1) )
    # now copy phase
    dicomPhaseFile = "%s/s%d/%s" % (ExamPath,phasID,  phaseFile  ) 
    fin,fout = os.popen4( "%s %s '%s' %s 2> /dev/null " % (printTagExe,dicomPhaseFile,instanceTag,dictionary) )
    result = fout.read()
    phaseTimeID = int(result.split("=")[1])
    os.system( "cp %s %s/Processed/s%d/i%d.MRDC.%d" % (dicomPhaseFile,ExamPath,dirID,dirID+2*phaseTimeID,2*phaseTimeID) )
    print  "%s %d %s %d " % (magnitudeFile,magTimeID,phaseFile,phaseTimeID) 


# one based
for fileid in range(1,numFiles+1):
    os.system( "%s -M%s/Processed/s%d/i%d.MRDC.%d -P%s/Processed/s%d/i%d.MRDC.%d --output %s/Processed/s%d/image -timeid %d" % ( 
                  DicomToComplex ,
                      ExamPath,dirID,dirID+2*fileid-1,2*fileid-1, 
                      ExamPath,dirID,dirID+2*fileid  ,2*fileid  , 
                      ExamPath,dirID,fileid-1 )  ) 
