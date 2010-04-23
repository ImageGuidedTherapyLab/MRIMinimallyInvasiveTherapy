import sys
import os
import ConfigParser
#the jobid is passed in from the command line
controlfile = sys.argv[1]
#controlfile = "/mnt/data/scratch/fuentes/DDDAS/data/biotex/090318_751642_treat/controlBrain.ini"

#read the main control file
iniFile  = ConfigParser.ConfigParser()
iniFile.readfp( open("%s" % controlfile , "r") )

# get exampath info
ExamPath   = iniFile.get(   "mrti" ,"exampath"    )
magID      = iniFile.getint("mrti" ,"magnitudeid" )
phaseID    = iniFile.getint("mrti" ,"phaseid"     )
#get dicom Dictionary
try:
  dictionary= iniFile.get(   "mrti" ,"dictionary"  )
except ConfigParser.NoOptionError:
  raise IOError("\n\n    No Dicom Dictionary FOUND! " )
instanceTag = "0020|0013"

# TODO wrap c++ code to return  value directly
printTagExe = "/home/dfuentes/DDDAS/trunk/utilities/DicomImageReadReturnIntTag"
try:
  DicomToComplex = iniFile.get(   "compexec" ,"Dicomtocomplex"  )
except ConfigParser.NoOptionError: 
  DicomToComplex = "~/DDDAS/trunk/utilities/DicomReadWriteComplexImage"

# make new directory
os.system('mkdir -p %s/Processed/s%d' % (ExamPath,magID) )

magnitudeList = os.listdir( "%s/s%d" % (ExamPath, magID ) )
numFiles = len(magnitudeList) 
phaseList     = os.listdir( "%s/s%d" % (ExamPath,phaseID) ) 
try:
   assert numFiles == len(phaseList)     
except AssertionError:
   raise IndexError("\n\n    error with number of files in the directory " )

# set output to magnitude id
dirID = magID

# TODO add realtime abilities as data is written to the ExamPath
for magnitudeFile,phaseFile in zip(magnitudeList,phaseList): 
    # copy magnitude first
    dicomMagFile   = "%s/s%d/%s" % (ExamPath, magID ,magnitudeFile) 
    fin,fout = os.popen4( "%s %s '%s' %s 2> /dev/null " % (printTagExe,dicomMagFile,instanceTag,dictionary) )
    result = fout.read()
    magTimeID = int(result.split("=")[1])
    os.system( "cp %s %s/Processed/s%d/i%d.MRDC.%d" % (dicomMagFile,ExamPath,dirID,dirID+2*magTimeID-1,2*magTimeID-1) )
    # now copy phase
    dicomPhaseFile = "%s/s%d/%s" % (ExamPath,phaseID,  phaseFile  ) 
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
