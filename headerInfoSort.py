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
ExamPath   = iniFile.get(   "mrti" ,"exampath"  )
dirID      = iniFile.getint("mrti" ,"dirid"     )
dictionary = iniFile.get(   "mrti" ,"dictionary")
instanceTag = "0020|0013"
# TODO wrap c++ code to return  value directly
printTagExe = "/home/dfuentes/DDDAS/trunk/utilities/DicomImageReadReturnIntTag"

# make new directory
os.system('mkdir -p %s/Sorted/s%d' % (ExamPath,dirID) )

# TODO add realtime abilities as data is written to the ExamPath
for file in os.listdir( "%s/s%d" % (ExamPath,dirID) ):
    dicomFile = "%s/s%d/%s" % (ExamPath,dirID,file) 
    fin,fout = os.popen4( "%s %s '%s' %s 2> /dev/null " % (printTagExe,dicomFile,instanceTag,dictionary) )
    result = fout.read()
    timeid = int(result.split("=")[1])
    print  "%s %d" % (file,timeid) 
    os.system( "cp %s %s/Sorted/s%d/i%d.MRDC.%d" % (dicomFile,ExamPath,dirID,dirID+timeid,timeid) )
