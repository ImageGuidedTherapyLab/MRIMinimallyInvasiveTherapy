"""simulatemrti.py script

this script is to simulate an mri machine writing files to disk 

the place the files get written to is determined from the control.ini file

"""
print __doc__
# Numerica imports intrinsic math functions(zeros,cos,sin,etc...)
#from Numeric import *
import ConfigParser
import sys
import os
import time

controlfile = sys.argv[1]

config = ConfigParser.ConfigParser()

#read the conf file
fname = open(controlfile,"r")
config.readfp(fname)


# get the local location/file name information for MRTI data
MRTI_filelocation=config.get("mrti","filelocation")
MRTI_nslice=config.getint("mrti","nslice")
MRTI_filebase=config.get("mrti","filebase")
MRTI_filemid=config.get("mrti","filemid")
MRTI_filesuffix=config.get("mrti","filesuffix")
MRTI_ntime=config.getint("mrti","ntime")
# get image aquisition time
acquisitiontime=config.getfloat("timestep","IDEAL_DT") # delay the copy by 
os.system('mkdir -p %s  ' %  MRTI_filelocation )

#assume that files are located in the $CWD
MRTIloc='.'
#transfer files
delay=2
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ MRTI DATA @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
for i in range(MRTI_ntime+1): # range stops at ifileend
   for j in range(1,MRTI_nslice+1): # range stops at MRTI_nslice
      filetmp='%s%d%s%d.%s' % \
         (MRTI_filebase,j,MRTI_filemid,i,MRTI_filesuffix)
      cpcmd='cp %s/%s %s/%s ' % (MRTIloc,filetmp,MRTI_filelocation,filetmp)
      os.system(''.join(['echo ',cpcmd])) 
      t0=time.time() ; status=os.system(cpcmd) ; t1=time.time()
      if(status==0):
        os.system("echo '\t\t\t%s copy complete in %f secs'"% (filetmp,t1-t0))
      else:
        os.system("echo '\t\t\t%s NOT COPIED'" % filetmp);time.sleep(delay)
   print 'image aquisition time = %f seconds' % acquisitiontime
   time.sleep(acquisitiontime) # delay the copy by the image aquisition 
                               # time to simulate mri machine writing to disk
