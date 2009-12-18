import sys
import os
import utilities
import jobsetup
import ConfigParser
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
  DirId    = iniFile.getint("mrti" ,"dirid")
  OutputDir= "mrivis"
  #get dicom Dictionary
  try:
    Dictionary= iniFile.get("mrti" ,"dictionary")
  except ConfigParser.NoOptionError:
    raise IOError("\n\n    No Dicom Dictionary FOUND! " )
  #noise estimate
  magIx    = iniFile.getint("kalman" ,"magix")
  magIy    = iniFile.getint("kalman" ,"magiy")

  # constant runtime options
  base_options = " %s %d %s %s %s -magIx %d -magIy %d" % \
               (ExamPath,DirId,Dictionary,OutputDir,runtime_options,magIx,magIy)
  
  #build list of jobs to run
  JOBS=jobsetup.setupkalman(iniFile) 

elif(Executable == "thermalTherapy" or Executable == "dddas"):

  #build list of jobs to run
  JOBS=jobsetup.setuplitt(iniFile) 

  # constant runtime options
  base_options = runtime_options

else:
  raise ValueError("unknown executable %s " % Executable )

# loop over the job in the list JOBS and run the code for each one
CODEEXEC=[] 
for (namejob,numproc,param_options,cntrlfile,method) in JOBS:
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
      qsubfile.write("#$ -pe 4way %d           \n" % numproc )
      qsubfile.write("#$ -N %s                 \n" % namejob )
      qsubfile.write("#$ -cwd                  \n"           )
      qsubfile.write("#$ -S /bin/bash          \n"           )
      qsubfile.write("#$ -v LD_LIBRARY_PATH,PATH,WORK,COMPILER\n")
      qsubfile.write("#$ -v MPI_VERSION,METHOD=%s \n" % method)
      qsubfile.write("mpirun -np $NSLOTS -machinefile $TMP/machines $WORK/exec/%s_$COMPILER-$MPI_VERSION-cxx-$METHOD  %s %s\n" % (Executable,base_options, param_options))
      qsubfile.write("pkill -9 dddas \n")
      # ensure entire file written before continuing
      qsubfile.close; qsubfile.flush() 
      execcode="cd %s/%s/%s ; qsub %s.qsub  " %  \
                         (workdir,jobid,namejob,namejob)
   # ranger
   elif(comphost.split(".")[0] == "login3"):
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
      qsubfile.write("#$ -q normal                   \n"             )
      qsubfile.write("#$ -pe 8way %d                      \n" % numproc   )
      qsubfile.write("#$ -l h_rt=02:00:00                 \n"             )
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
      execcode="cd %s/%s/%s ; mpirun -n %d $WORK/exec/%s_$COMPILER-$MPI_VERSION-cxx-$METHOD  %s %s" % (workdir,jobid,namejob,numproc,
                       Executable,base_options,param_options)
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

#run the code
os.system(execMETH)

