import sys
import os
import utilities
import jobsetup
import ConfigParser
import time
import subprocess
#the jobid is passed in from the command line
controlfile = sys.argv[1]

#read the main control file
iniFile = utilities.control_defaults(controlfile) 
iniFile.readfp( open("%s" % controlfile , "r") )

#get the jobid
checkid = iniFile.get("compexec","jobid")
jobid   = utilities.checkprevious(checkid)
iniFile.set("compexec","jobid",jobid)
Executable = iniFile.get( "compexec" ,"executable")

print iniFile.get("compexec","comphost")
#get the computation host information
[comphost,loccomphost,compendian,workdir] = \
         utilities.hostinformation(iniFile,"compexec","comphost","workdir")

#constant runtime options
runtime_options = iniFile.get("compexec","runtime_options")

if(Executable == "image"):

  #get input/output parameters
  ExamPath = iniFile.get(   "mrti" ,"exampath")
  listDirId= map(int,iniFile.get("mrti","dirid").split(";"))
  OutputDir= "mrivis"
  #noise estimate
  magIx    = iniFile.getint("kalman" ,"magix")
  magIy    = iniFile.getint("kalman" ,"magiy")

  unsortedJOBS=[]
  for DirId  in listDirId: 
    # runtime options
    BaseOptions = " %s/Processed %d %s %s -magIx %d -magIy %d " % \
                 (ExamPath,DirId,OutputDir,runtime_options,
                  magIx,magIy)

    #build list of jobs to run
    unsortedJOBS = unsortedJOBS + jobsetup.setupkalman(iniFile,BaseOptions,DirId) 
  # hack resort by spread 
  # put faster jobs first
  JOBS= []
  for job in  unsortedJOBS: 
      if( int(job[0][len(job[0])-2:])%3 == 0) : 
        JOBS.insert(0,job) # insert at front
      else:
        JOBS.append(job) # insert at back
  #print JOBS

elif(Executable == "thermalTherapy" or Executable == "dddas"):

  # constant runtime options
  BaseOptions = runtime_options

  #build list of jobs to run
  JOBS=jobsetup.setuplitt(iniFile,BaseOptions) 

else:
  raise ValueError("unknown executable %s " % Executable )

# loop over the job in the list JOBS and run the code for each one
Executable  = "dddas"
CODEEXEC=[] 
for (namejob,numproc,base_options,param_options,cntrlfile,method) in JOBS:
   # code execution on mda cluster
   if(comphost.split(".")[0] == "cn285" 
                    or 
      comphost.split(".")[0] == "cn286" ):
      execcode="cd %s/%s/%s ; bsub -J  %s -n %d -o out.o -e err.o /opt/hpmpi/bin/mpirun -d -prot -srun $WORK/exec/%s_$COMPILER-$MPI_VERSION-cxx-$METHOD  %s %s" % (workdir,jobid,namejob,namejob,numproc,
                       Executable,base_options,param_options)
   # code execution on shamu
   elif(comphost.split(".")[0] == "shamu"):
      # write a qsub file
      qsubfile=open("%s/%s/%s/%s.qsub" %(workdir,jobid,namejob,namejob) ,"w")
      qsubfile.write("#!/bin/bash              \n"           )
      qsubfile.write("#$ -pe mpich %d           \n" % numproc )
      qsubfile.write("#$ -N %s                 \n" % namejob )
      qsubfile.write("#$ -cwd                  \n"           )
      qsubfile.write("#$ -S /bin/bash          \n"           )
      qsubfile.write("#$ -v LD_LIBRARY_PATH,PATH,WORK,COMPILER\n")
      qsubfile.write("#$ -v MPI_VERSION,METHOD \n" )
      qsubfile.write("mpirun_rsh -np $NSLOTS -hostfile $TMP/machines $WORK/exec/%s_$COMPILER-$MPI_VERSION-cxx-$METHOD  %s %s\n" % (Executable,base_options, param_options))
      # ensure entire file written before continuing
      qsubfile.close; qsubfile.flush() 
      execcode="cd %s/%s/%s ; qsub %s.qsub  " %  \
                         (workdir,jobid,namejob,namejob)
   # ranger
   elif(comphost.split(".")[0] == "login3" or 
        comphost.split(".")[0] == "login4" ):
      # write a qsub file
      qsubfile=open("%s/%s/%s/%s.qsub" %(workdir,jobid,namejob,namejob) ,"w")
      qsubfile.write("#!/bin/bash                           \n"           )
      qsubfile.write("# Which account to be charged cpu time\n"           )
      qsubfile.write("#$ -A UTMDACC-DIP                     \n"           )
      qsubfile.write("#  combine stdout stderr              \n"           )
      qsubfile.write("#$ -j y                               \n"           )
      qsubfile.write("#  jobname                            \n"           )
      qsubfile.write("#$ -N %s                              \n" % namejob )
      qsubfile.write("#  inherit submission env             \n"           )
      qsubfile.write("#$ -V                                 \n"           )
      qsubfile.write("# The job is located in the current   \n"           )
      qsubfile.write("# working directory.                  \n"           )
      qsubfile.write("#$ -cwd                             \n\n"           )
      qsubfile.write("#$ -o $JOB_NAME.o$JOB_ID            \n"             )
      qsubfile.write("#$ -q normal                        \n"             )
      qsubfile.write("#$ -pe 16way %d                     \n" % numproc   )
      qsubfile.write("#$ -l h_rt=10:00:00                 \n"             )
      qsubfile.write("set -x                              \n"             )
      qsubfile.write("ibrun ${WORK}/exec/%s_${PETSC_ARCH} %s %s \n" % \
                                              (Executable,base_options,param_options))
      # ensure entire file written before continuing
      qsubfile.close; qsubfile.flush() 
      execcode="cd %s/%s/%s ; qsub %s.qsub  " %  \
                         (workdir,jobid,namejob,namejob)
   elif(comphost.split(".")[0] == "lslogin1"):
      bsubbase = cntrlfile.get(   "compexec","bsub")
      execcode="cd %s/%s/%s ; %s -J %s -n %d -o out.o -e err.o ibrun /work/utexas/iv/fuentes/exec/dddas_em64t-cxx %s " %  \
                         (workdir,jobid,namejob,bsubbase,namejob,
                                        numproc,runtime_options)
   else: # default code execution
      if(numproc > 1): 
        execcode="cd %s/%s/%s; mpirun -n %d $WORK/exec/%s_$COMPILER-$MPI_VERSION-cxx-$METHOD  %s %s > log.txt " % (workdir,jobid,namejob,numproc,
                        Executable,base_options,param_options)
      else:
        execcode="cd %s/%s/%s; $WORK/exec/%s_$COMPILER-$MPI_VERSION-cxx-$METHOD %s %s > log.txt " % (workdir,jobid,namejob,Executable,base_options,param_options)
   # write control file with additional parameters
   inifile=open("%s/%s/files/control.ini" % (jobid,namejob) ,"w")
   cntrlfile.set("compexec","execcode",execcode)
   cntrlfile.write(inifile)
   inifile.close
   inifile.flush() # ensure the entire file is written before continuing
   CODEEXEC.append(execcode)

execMETH= ";".join(CODEEXEC)
print "code execution method: %s \n" % execMETH

# sit idle until user inputs ready to continue
utilities.pause_until_ready()

# create an error log file
errorFile=open("%s/error.log" %(jobid) ,"w")
#run the code
numlocalProc = 12 
process = []
while ( len(CODEEXEC) > 0 or len(process) > 0 ):
    # only run numlocalProc at a time
    if (len(process) > numlocalProc):
      raise RuntimeError("\n\n running too many jobs at a time??")
    elif (len(process) == numlocalProc):
      print len(CODEEXEC), " jobs remaining..."
      time.sleep(30) # pause wait for jobs to finish
    elif( len(CODEEXEC) > 0 ):
      cmd = CODEEXEC.pop(0)
      print "running " , cmd
      process.append( [subprocess.Popen(cmd,shell=True),cmd] )
    if( len(process) > 0 ):
      runningJob = process.pop(0)
      if ( runningJob[0].poll() == None ):
        # job not done put it back in the list
        # not that we pop from the front of the list and pushback at the 
        # end to cycle through
        print " pid ",runningJob[0].pid, " still running"
        process.append( runningJob )
      elif ( runningJob[0].poll() == 0 ):
        pass # job is done 
      else:
        print "job exiting with ", runningJob[0].poll() 
        errorFile.write("error in  %s   \n" % runningJob[1] )
        #raise RuntimeError("\n\n unknown exit code ")
errorFile.close; errorFile.flush() 
