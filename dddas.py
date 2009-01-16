"""dddas.py script

usage: 
      
   python dddas.py control.ini

   where control.ini exists in the $CWD and contains:
         1) all data needed to gather job execution info 
         2) all data needed to gather file transfer info
         3) relative location of mesh file or directory of mesh files  
         4) relative location of power file or directory of power files  
   
requires that:
   ForwardAgent yes
   ForwardX11 yes
   be set in the file /etc/ssh/ssh_config
          

The script builds the necessary directory hierarchy and pre-processes the 
control files necessary for (batch) job submission. 

It looks at the list of jobs to be run and builds the scripts to transfer
the MRTI/vis data. For Data transfer, does not keep track of restart files. 
ASSUME one data transfer of FEM visualization files per optimization step

*NOTE* to restart a job on the computational host
*NOTE* login to the computational host and use the code execution method
*NOTE* stored in the control files of the job of interest
"""
print __doc__
# Numerica imports intrinsic math functions(zeros,cos,sin,etc...)
import sys
import os
import socket
import utilities
import jobsetup
import ConfigParser
import time

if (len(sys.argv) < 2 ):
    raise __doc__
    

#command to check that file exists
waitcmd='! while [ ! -f %s ]; do echo waiting for %s ; sleep %d; done \n'
#command to check that full file is present on disk (checks the file size)
checkcmd='! while [ $((`du -b %s | cut -f1`)) -lt $((`cat %s`)) ]; do echo file size `du -b %s | cut -f1` expected size `cat %s` ; sleep %d; done \n'
#pause time 
delay=2

#the time will be used to tag each file
profileID = "%d" % int(time.time())
Xferfile = "transfer_times%s.txt" % profileID

#  the following function checks that the file exists and the full data
#  has been written to disk before sending 
#     the variables: location,destination,mainfile
#        are set in the scope before the function is called
def waitandsend(File,Sizefile='NOfile',Tag=False):
    if(Tag):
      fullmain='%s/%s%s'%(location,profileID,mainfile)
    else:
      fullmain='%s/%s'%(location,mainfile)
    fullsize='%s/%s'%(location,Sizefile)
    # checks that file exists
    File.write('! date +%s > TIME0 \n')
    File.write(waitcmd % (fullmain,fullmain,delay) )
    # checks that full file is present on disk
    if(Sizefile!='NOfile'):
       File.write(waitcmd % (fullsize,fullsize,delay) )
       File.write(checkcmd % (fullmain,fullsize,fullmain,fullsize,delay) )
    File.write('! date +%s > TIME1 \n')
    # transfer main file 
    File.write('put %s %s \n' % (fullmain,destination))
    File.write('! date +%s > TIME2 \n')
    # report time
    File.write('! echo %s time wait $((`cat TIME1`-`cat TIME0`))s  transferred in $((`cat TIME2`-`cat TIME1`))s >> %s \n' % (mainfile,Xferfile))
    return

#the jobid is passed in from the command line
controlfile = sys.argv[1]

# initialize config file defaults
dddas = utilities.control_defaults(controlfile) 
#read the main dddas control file
dddas.readfp( open("%s" % controlfile , "r") )

#quick error check to ensure that an initial laser position is given
try:
  x_0    =dddas.getfloat("source_laser","x_0")
  y_0    =dddas.getfloat("source_laser","y_0")
  z_0    =dddas.getfloat("source_laser","z_0")
except ConfigParser.NoOptionError:
  assert False, "\n\n    Thermal Source DATA NOT FOUND !!!!!!!  "

#get the jobid
checkid = dddas.get("compexec","jobid")
jobid   = utilities.checkprevious(checkid)
dddas.set("compexec","jobid",jobid)
dddas.set("compexec","profileID",profileID)

#get visualase parameters
visualase_file    = dddas.get("visualase"    ,"visualase_file")
visualase_dir     = dddas.get("visualase"    ,"visualase_dir" )

# get the local location/file name information for MRTI data
MRTI_filelocation = dddas.get(      "mrti"  ,"filelocation")
os.system('mkdir -p %s' % MRTI_filelocation) # ensure directory exists
MRTI_nslice       = dddas.getint(   "mrti"  ,"nslice")
MRTI_filebase     = dddas.get(      "mrti"  ,"filebase")
MRTI_filemid      = dddas.get(      "mrti"  ,"filemid")
MRTI_filesuffix   = dddas.get(      "mrti"  ,"filesuffix")

# get the file name information for fem data
fem_filename=dddas.get("output","fem_filename")

#get the username
username = dddas.get("compexec","username")

#get the visualization host information
[vishost,locvishost,visendian,viswork] = \
         utilities.hostinformation(dddas,"output","vishost","viswork")

#get the computation host information
[comphost,loccomphost,compendian,workdir] = \
         utilities.hostinformation(dddas,"compexec","comphost","workdir")

# determine if byte swapping needed
if( compendian == visendian):
   byteswap = False
   dddas.set("output","byteswap","false") 
else:
   byteswap = True
   dddas.set("output","byteswap","true") 

# set the MRTI file location on the comphost
if( not loccomphost ): 
  dddas.set("mrti" ,"filein"   ,"%s/mridat/%s%%d%s%%d.%s" % 
        (workdir,MRTI_filebase,MRTI_filemid,MRTI_filesuffix) )
else:
  dddas.set("mrti" ,"filein"   ,"%s/%s/%s%%d%s%%d.%s" % 
        (workdir,MRTI_filelocation,MRTI_filebase,MRTI_filemid,MRTI_filesuffix) )
         

#setup the job:
#   0) create directory heirarchy
#   1) build a list of jobs to run
#   2) copy meshfile(s)/powerfile(s) to correct position in directory hierarchy
JOBS=jobsetup.setupjob(dddas) 

#get information to determine which files are written and/or will be sent
MRTI_transfer = dddas.getboolean("mrti","transfer")
femavs        = dddas.getboolean("output","femavs")
femgmv        = dddas.getboolean("output","femgmv")
femraw        = dddas.getboolean("output","femraw")
MRTIavs       = dddas.getboolean("output","mrtiavs")
MRTIrawiv     = dddas.getboolean("output","mrtirawiv")
 
#get the time instances of the MRTI data 
MRTI_nzero = dddas.getint("mrti","nzero")
MRTI_ntime = dddas.getint("mrti","ntime")
if(MRTI_ntime < MRTI_nzero): 
   raise "\n\n    MRTI_ntime less than MRTI_nzero"
   

#echo data transfer info
if(MRTI_transfer): 
   print "transfering original MRTI data"
   raise "\n\n    haven't added code to transfer the header data in mrti.ini"
   MRTIavs  = True  #should be true transferring
if(femavs): #assume sending 2 avs files per ideal step (1 fem 1 params) 
   print "transfering fem data (avs format)"
if(femgmv): #assume sending 1 gmv file per ideal step
   print "transfering fem data (gmv format)"
if(femraw): #assume sending 1 raw file per ideal step
   print "transfering fem data (raw format)"
if(MRTIavs): #assume sending one MRTI avs file per ideal step
   print "transfering thermal imaging data(avs format)"
if(MRTIrawiv): #assume sending one MRTI rawiv file per ideal step
   print "transfering thermal imaging data(rawiv format)"


#determine code execution method
CODEEXEC=[]
# loop over the jobs names in the list JOBS and run the code for each one
for jobiter in JOBS:
   namejob   = jobiter[0] # the jobs name
   cntrlfile = jobiter[1] # the control file for this job
   #determine the number of processors to run on and execution method
   numproc  = cntrlfile.getint("compexec","numproc")
   run      = cntrlfile.get(   "compexec","run")
   # code execution on lonestar
   if(comphost.split(".")[0] == "lonestar"):
      bsubbase = cntrlfile.get(   "compexec","bsub")
      execcode="cd %s/%s/%s ; %s -J %s -n %d -o out.o%s -e err.o%s %s " %  \
                         (workdir,jobid,namejob,bsubbase,namejob,
                                        numproc,profileID,profileID,run)
   # code execution on shamu
   elif(comphost.split(".")[0] == "shamu"):
      # write a qsub file
      qsubfile=open("%s/%s/%s/%s.qsub" %(workdir,jobid,namejob,namejob) ,"w")
      qsubfile.write("#!/bin/bash              \n"           )
      qsubfile.write("#$ -pe mpich %d          \n" % numproc )
      qsubfile.write("#$ -N %s                 \n" % namejob )
      qsubfile.write("#$ -cwd                  \n"           )
      qsubfile.write("#$ -S /bin/bash          \n"           )
      qsubfile.write("echo 'Got $NSLOTS slots' \n"           )
      qsubfile.write("echo $TMP                \n"           )
      qsubfile.write("export LD_LIBRARY_PATH=%s\n" %  os.getenv('LD_LIBRARY_PATH')         )
      qsubfile.write(run)
      # ensure entire file written before continuing
      qsubfile.close; qsubfile.flush() 
      execcode="cd %s/%s/%s ; qsub %s.qsub  " %  \
                         (workdir,jobid,namejob,namejob)
   else: # default code execution
      execcode="cd %s/%s/%s ; %s " % (workdir,jobid,namejob, run % numproc)
   CODEEXEC.append(execcode)
   # command for vis file transfer
   VISxferlocation=''.join([vishost,":",viswork,'/',jobid])
   vistransfermri="sftp -b %s/%s/%s/VisMriFile.txt %s" % \
                                 (workdir,jobid,namejob,VISxferlocation)
   vistransferfem="sftp -b %s/%s/%s/VisFemFile%%d.txt %s" % \
                                 (workdir,jobid,namejob,VISxferlocation)
   cntrlfile.set("compexec","vistransfermri",vistransfermri)
   cntrlfile.set("compexec","vistransferfem",vistransferfem)
   # store execution method to restart
   cntrlfile.set("compexec","execcode",execcode)
   # tell avs where to look for files
   cntrlfile.set("avs","mrti","%s/%s/%s/mrivis/%s%s%s%%04d.fld" 
             %  (viswork,jobid,namejob,profileID,MRTI_filebase,visendian) )
   num_qoi=cntrlfile.getint("compexec","num_qoi")
   cntrlfile.set("avs","power","%s/%s/%s/mrivis/%s%s_qoi_%d%s%%04d.fld" 
             %  (viswork,jobid,namejob,profileID,"power",num_qoi,visendian) )
   cntrlfile.set("avs","fem","%s/%s/%s/femvis/%s%sqoi_%%dnopt_%%d%s%%04d.inp"  
             % (viswork,jobid,namejob,profileID,fem_filename,visendian ) )
   cntrlfile.set("avs","param","%s/%s/%s/femvis/%s%sqoi_%%dnopt_%%d%s.inp" 
             % (viswork,jobid,namejob,profileID,fem_filename,"param") )
   # set default anatomical mri data location for avs
   interop_dicom_dir=cntrlfile.get("registration","interop_dicom_dir")
   cntrlfile.set("avs","mri","%s/%s%s%%04d.fld"%(viswork,
               filter(len,interop_dicom_dir.split("/")).pop() , visendian))

   # overwrite control file with additional parameters
   inifile=open("%s/%s/files/control.ini" % (jobid,namejob) ,"w")
   cntrlfile.write(inifile)
   inifile.close
   inifile.flush() # ensure the entire file is written before continuing


#determine file transfer method and computational host code execution method
if(loccomphost): #local computations
   xfermethDAT='ls' 
   xferdestDAT=''.join([workdir,'/',jobid]) 
   execMETH= ";".join(CODEEXEC)
else:   #remote computations and local OR remote visualization
   xfermethDAT='sftp -o ControlPath=/tmp/%r@%h:%p -b mriMRTI.txt' ; 
   xferdestDAT=''.join([username,"@",comphost,":",workdir])
   # need to set environment for lonestar
   if(comphost.split(".")[0] == "lonestar"):
      execMETH="ssh -S /tmp/%%r@%%h:%%p %s@%s ' . /opt/modules/init/sh ; . /usr/local/etc/profile; \
                          %s '" %  (username,comphost, ";".join(CODEEXEC))
   else: # default execution method
      execMETH="ssh -S /tmp/%%r@%%h:%%p %s@%s ' %s ' "% (username,comphost,";".join(CODEEXEC))

datxfer='%s %s ' % (xfermethDAT,xferdestDAT)
   


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# write initial visualase file 
vislasfile=open("%s/%s" %(jobid,visualase_file) ,"w")
vislasfile.write("0.0 \n")
vislasfile.close
vislasfile.flush() # ensure the entire file is written before continuing

# write out the filesize file
#MRTIimg_sizefile    = dddas.get("mrti","img_sizefile")
#MRTIfilesize=open("%s/%s" % (MRTI_filelocation,MRTIimg_sizefile) ,"w")
#MRTI_xpixel = dddas.getint(   "mrti"  ,"xpixel")
#MRTI_ypixel = dddas.getint(   "mrti"  ,"ypixel")
#MRTIfilesize.write("%d \n" % (MRTI_xpixel * MRTI_ypixel) ) # assume byte data
#MRTIfilesize.close
#MRTIfilesize.flush() # ensure the entire file is written before continuing

#create DATA transfer script for MRI/MRTI data from local host to comp host
if( not loccomphost ):  # do not transfer if are doing local computations
   mriMRTI=open("mriMRTI.txt","w")
   # initialize file to record data transfer times
   mriMRTI.write('! echo data transfer from %s to %s > %s \n' %
                     ( socket.gethostname().split(".")[0],comphost,Xferfile) )
   # in case file name gets messed up
   mriMRTI.write("""! echo 'to change file names, use:            sed "s/%s\([0-9]\)%s/FILEBASE\\1FILEMID/g" mriMRTI.txt > new.txt' \n""" %
                                                (MRTI_filebase,MRTI_filemid) )
   # remember how to execute new batch sftp script
   mriMRTI.write('! echo sftp -o ControlPath/tmp/%%r@%%h:%%p -b new.txt %s@%s:%s \n' % 
                                                 (username,comphost,workdir) ) 
   #transfer MRTI DATA 
   #   when an image representing t_i is sent the power
   #   that should be output for [t_i,t_{i+1}) is immediately
   #   retrieved by the MRTI server. In order for the immediate
   #   retrieval of the power, the power file for [t_i,t_{i+1}) must
   #   be written before the the image at t_i is sent. Hence, after
   #   reading and writing the image at t_i, the power for
   #   [t_{i+1},t_{i+2}) is written so it can already be there ready
   #   to be retrieved by the MRTI server. NOTE that the power stored
   #   in module constit_data in POW is of the form
   #  
   #                 line 1     line 2 ... line n-1       line n
   #     |----P(t1)----|---P(t2)---|----------|-----P(tn)-----|
   #    t=0           t1          t2   ...   tn-1             tn
   #  
   #    if the image just read in represents time t_i then the
   #    power that should be output is for [t_i,t_{i+1}) which
   #    is stored under index IDfile={i+1}

   if(MRTI_transfer):
     for i in range(MRTI_ntime+1): # range stops at n - 1
       for j in range(1,MRTI_nslice+1): # range stops at n - 1
         location=MRTI_filelocation ; destination='mridat'
         mainfile='%s%d%s%d.%s'%(MRTI_filebase,j,MRTI_filemid,i,MRTI_filesuffix)
         waitandsend(mriMRTI,MRTIimg_sizefile)
         mriMRTI.write('get %s/%s %s/mlat.tmp \n' % 
                              (jobid,visualase_file,visualase_dir)         )
         mriMRTI.write('!  mv %s/mlat.tmp %s/%s \n' % 
                              (visualase_dir,visualase_dir,visualase_file) )
   mriMRTI.close
   mriMRTI.flush() # ensure the entire file is written before continuing

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
#@@@@@@@@ loop through list of control files and create 
#@@@@@@@@ ONE MRTI DATA transfer scripts and 
#@@@@@@@@ num_qoi FEM vis data transfer scripts
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
# get the write directory of the computations
compwritefem        = dddas.get("output","compwritefem")
compwritemri        = dddas.get("output","compwritemri")

# get file names that contain the file size info
MRTIvis_sizefile    = dddas.get("avs","MRTIvis_sizefile")
femvis_sizefile     = dddas.get("avs","femvis_sizefile")
#count the number of vis files that will be sent
numvisfiles=0

# transfer files if are doing remote computations or remote visualization
if( not locvishost or not loccomphost ): 
  # loop over the jobs names in the list JOBS and transfer files for each one
  for jobiter in JOBS:
    namejob   = jobiter[0] # the jobs name
    cntrlfile = jobiter[1] # the control file for this job

    #get the # of qoi's for this job
    num_qoi=cntrlfile.getint("compexec","num_qoi")
    print "%d quantities of interest in %s" % (num_qoi,namejob)

    # open mrti transfer file
    VisMriFile=open("%s/%s/VisMriFile.txt" % (jobid,namejob) ,"w")
    # initialize file to record data transfer times
    VisMriFile.write('! echo data transfer from %s to %s > %s \n' %
                                        (comphost,vishost,Xferfile) )
    # the location and destination are the same for all mri vis files
    if(compwritemri == 'mrivis/'):
       location='%s/%s/%s/mrivis' % (workdir,jobid,namejob) 
    else:
       location=compwritemri
    destination='%s/mrivis' % namejob
    #first must transfer info about the size of the files to be transfered
    #send MRTI size file
    if(MRTIavs):
      mainfile=MRTIvis_sizefile
      waitandsend(VisMriFile)
    # send MRTI data
    for i in range(MRTI_nzero,MRTI_ntime+1): # range stops at n - 1 
       #send power file first
       mainfile='power_qoi_%d%s%04d.fld'%(num_qoi,visendian,i) 
       waitandsend(VisMriFile,"NOfile",True)
       numvisfiles= numvisfiles+1 # update counter
       if(MRTIavs):
          mainfile='%s%s%04d.fld'%(MRTI_filebase,visendian,i) 
          waitandsend(VisMriFile,MRTIvis_sizefile,True)
          numvisfiles= numvisfiles+1 # update counter
       if(MRTIrawiv):
          mainfile='%s%s%04d.rawiv'%(MRTI_filebase,visendian,i)
          waitandsend(VisMriFile,'NOfile',True)
          numvisfiles= numvisfiles+1 # update counter
    # finish by sending power files from remaining groups
    for idgroup in range(num_qoi): # range stops at n - 1
      #get the timestep and optimization step info for each qoi
      noptsteps   = cntrlfile.getint(    "qoi_%d" % idgroup ,"noptsteps" )
      for idopt in range(noptsteps): # range stops at n - 1
          mainfile='power_qoi_%d%s%04d.fld' % (idgroup,visendian,idopt) 
          waitandsend(VisMriFile,'NOfile',True)
    VisMriFile.close # close mrti transfer file
    VisMriFile.flush() # ensure the entire file is written before continuing

    # set the location and destination for the fem vis transfer
    if(compwritefem == 'femvis/'):
       location='%s/%s/%s/femvis' % (workdir,jobid,namejob)  
    else:
       location=compwritefem
    destination='%s/femvis' % namejob

    # loop over groups
    for idgroup in range(num_qoi): # range stops at n - 1
      #get the timestep and optimization step info for each qoi
      noptsteps   = cntrlfile.getint(    "qoi_%d" % idgroup ,"noptsteps" )
      plotoptsolve= cntrlfile.getboolean("qoi_%d" % idgroup ,"plotoptsolve" )
      noffset     = cntrlfile.getint(    "qoi_%d" % idgroup ,"noffset" )
      numideal    = cntrlfile.getint(    "qoi_%d" % idgroup ,"numideal" )
      ideal_nzero = cntrlfile.getint(    "qoi_%d" % idgroup ,"ideal_nzero_init")
      ideal_ntime = cntrlfile.getint(    "qoi_%d" % idgroup ,"ideal_ntime_init")
      fem_ntime   = max(ideal_ntime + numideal * (noptsteps-1),MRTI_ntime)
      #create data transfer file if plotting
      if(plotoptsolve):
         VisFemFile=open("%s/%s/VisFemFile%d.txt" % (jobid,namejob,idgroup),"w")
         mainfile="qoi_%d%s" % (idgroup,femvis_sizefile)
         # first must transfer info about the size of the files to be transfered
         waitandsend(VisFemFile)
         for idopt in range(noptsteps): # range stops at n - 1
            #always begin by sending a FEM vis file to check 
            #laser position, fields, registration , and other parameters
            if(femavs): #sending 3 avs files used to vis the setup  
               mainfile='%sqoi_%dnopt_%d%sinit%04d.inp' % \
                                  (fem_filename,idgroup,idopt,visendian,idopt)
               sizefile="qoi_%d%s" % (idgroup,femvis_sizefile)
               waitandsend(VisFemFile,sizefile,True)
               mainfile='%sqoi_%dnopt_%d%sflds%04d.inp' % \
                                  (fem_filename,idgroup,idopt,visendian,idopt)
               sizefile="qoi_%d%s" % (idgroup,femvis_sizefile)
               waitandsend(VisFemFile,sizefile,True)
               mainfile='%sqoi_%dnopt_%dparam.inp' % \
                                  (fem_filename,idgroup,idopt)
               waitandsend(VisFemFile,'NOfile',True)
               numvisfiles= numvisfiles+3 # update counter
            # always plot out to the end only update the start time
            for idtime in range(ideal_nzero + idopt * noffset, fem_ntime+1):
               # send FEM/PARAMS data in AVS format
               if(femavs): #sending 1 avs file per step 
                  mainfile='%sqoi_%dnopt_%d%s%04d.inp' % \
                                   (fem_filename,idgroup,idopt,visendian,idtime)
                  sizefile="qoi_%d%s" % (idgroup,femvis_sizefile)
                  waitandsend(VisFemFile,sizefile,True)
                  numvisfiles= numvisfiles+1 # update counter
               # send FEM data in GMV format
               if(femgmv): #sending 1 gmv file per step
                  mainfile='%sqoi_%dnopt_%d.%04d' % \
                                            (fem_filename,idgroup,idopt,idtime) 
                  waitandsend(VisFemFile,'NOfile',True)
                  numvisfiles= numvisfiles+1 # update counter
               # send FEM data in RAW format
               if(femraw): #sending 1 raw file per step
                  mainfile='%sqoi_%dnopt_%d.%04d.raw' % \
                                            (fem_filename,idgroup,idopt,idtime) 
                  waitandsend(VisFemFile,'NOfile',True)
                  numvisfiles= numvisfiles+1 # update counter
         VisFemFile.close   # close fem transfer file for the qoi in this job
         VisFemFile.flush() # ensure the entire file written before continuing

print "vishost is %s  " % visendian
print "vishost is %s  " % visendian
print "vishost is %s  " % visendian
print "comphost is %s " % compendian  
print "comphost is %s " % compendian  
print "comphost is %s " % compendian  
if(byteswap): 
   print "BYTESWAPPING..." 
   print "BYTESWAPPING..." 
   print "BYTESWAPPING..." 

if(MRTI_transfer):
   nummrifiles = MRTI_ntime - MRTI_nzero + 1
else:
   nummrifiles = 0
print "transfering %d mri file sets and transfering %d vis file sets \n" % \
                                                       (nummrifiles,numvisfiles)

# echo the commands that will be executed 
print "\nthe following commands will be executed...\n\n"
print "code execution method: %s \n" % execMETH
print "MRTI file transfer method: %s \n" % datxfer

if(dddas.getboolean("output","computefinalobjective")):
   print """
COMPUTING FINAL OBJECTIVE FUNCTION VALUE FOR EACH OPTIMIZATION
STEP THIS WILL NOT WORK IN REAL-TIME!!!!!!!!!!!!!!
         """
try:
   if(dddas.getint("mrti","pixel_size")):
      print "TRANSFERING IMAGE OF SIZE %d BYTES PER PIXEL!!!!" % \
                                        dddas.getint("mrti","pixel_size")
except:
   pass


# create persistent connections if they do not exist...
#    ssh -O exit to kill a connection
if( not locvishost ): 
    if(os.system("ssh -O check -S /tmp/%%r@%%h:%%p %s@%s" % (username,vishost) ) ):
       print "\n\n   Creating Persisent Connection on %s!!!! \n\n " % vishost
       os.system('ssh -MNf -S /tmp/%%r@%%h:%%p %s@%s ' %  (username,vishost) )
    else:
       print "\n\n   Found Persisent Connection on %s!!!! \n\n " % vishost
if( not loccomphost ): 
    if(os.system("ssh -O check -S /tmp/%%r@%%h:%%p %s@%s" % (username,comphost) ) ):
       print "\n\n   Creating Persisent Connection on %s!!!! \n\n " % comphost
       os.system('ssh -MNf -S /tmp/%%r@%%h:%%p %s@%s ' %  (username,comphost) )
    else:
       print "\n\n   Found Persisent Connection on %s!!!! \n\n " % comphost
# sit idle until user inputs ready to continue
utilities.pause_until_ready()

#make directories recursively ignoring if the directory already exists
#  AND transfering all files under directory jobid
#  EVERYTHING needed to run code and should be present in directory jobid
if( not locvishost ): 
    os.system('ssh -S /tmp/%%r@%%h:%%p %s@%s mkdir -p %s  ' %    (username,vishost,viswork) )
    os.system('scp -o ControlPath=/tmp/%%r@%%h:%%p -r %s  %s@%s:%s' %  (jobid,username,vishost,viswork) )

#make directories recursively ignoring if the directory already exists
#  AND transfering all files under directory jobid
#  EVERYTHING needed to run code and should be present in directory jobid
if( not loccomphost ): 
    # create directory to store mri data
    os.system('ssh -S /tmp/%%r@%%h:%%p %s@%s mkdir -p %s/mridat' %  (username,comphost,workdir) ) 
    os.system('scp -o ControlPath=/tmp/%%r@%%h:%%p -r %s  %s@%s:%s' %     (jobid,username,comphost,workdir) )

#transfer MRI/MRTI/FEM files and run the code
os.system(execMETH)
os.system(datxfer)
   
