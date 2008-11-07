"""mesh_registration: (requires that ITK installed on reghost)

   This script is for the rigid registration of the FEM mesh 
   from pre-op data with intra-op mri data. 
   The script will transfer the fem mesh, pre-op data, and intra-op data
   to the regwork directory on reghost. Uses ITK to convert the 
   dicom data to mha/rawiv/fld files. 
   The mha files may be subsampled to improve code performance.
   Then the ITK rigid
   registration code is run on the possibly subsampled pre-op/intra-op data.
   ITK will read and write to the same control.ini file that avs uses
   to visualize the registration

   The required input: 
   1) the fem mesh
   2) pre-op dicom directory file from which the fem mesh was created
   3) the directory containing the intra-op mri data
"""
import os  # make os a global variable in this modules namespace
import sys
import utilities
import copy
import ConfigParser 


#the jobid is passed in from the command line
controlfile = sys.argv[1]

# initialize config file defaults
config = utilities.control_defaults(controlfile) 
#read the main dddas control file
config.readfp( open("%s" % controlfile , "r") )

#run some quick error checks to ensure that the thermal images and 
# laser are registered with the planning images
try:
  mrti_image_zero  =config.getint("planning_images" , "mrti_image_zero")        
  laser_xpixel     =config.getint("planning_images" , "laser_xpixel")        
  laser_ypixel     =config.getint("planning_images" , "laser_ypixel")        
  laser_slice_plane=config.getint("planning_images" , "laser_slice_plane") 
except ValueError:
  assert False, "\n\n    INCORRECT DATA INPUT FOR PLANNING IMAGES !!!!!!!  "
except ConfigParser.NoOptionError:
  assert False, "\n\n    PLANNING_IMAGE DATA NOT FOUND !!!!!!!  "

checkid   = config.get( "compexec" , "jobid" ) 
checkid   = "reg_%s" % checkid   
jobid     = utilities.checkprevious(checkid)

#get the username
username = config.get("compexec","username")

#meshdata and powerdata used to point to files
transform_script  = config.get("registration","transform_script" ) 
createbackupmesh= config.get("registration","createbackupmesh" ) 
unregistered_mesh = config.get("compexec","meshdata" ) 
meshfilename = filter(len,unregistered_mesh.split("/")).pop() 
try:
   fileext  = meshfilename.split(".")[1]
except IndexError:
   raise "MESH FILE TYPE IS EXPECTED TO HAVE AN EXTENSION" 
else:
   if(fileext != 'raw' and fileext != 'rawn'):
      raise "%s UNKNOWN MESH FILE TYPE" % fileext
# upon succesful read of mesh file overwrite the mesh file with
# its expected new location
config.set("compexec","meshdata","%s/input_compact" % jobid) 
# data directories
interop_dicom_dir = config.get("registration","interop_dicom_dir") 
preop_dicom_dir   = config.get("registration", "preop_dicom_dir" ) 
# path to suite of codes used for registration
registration_dir  = config.get("registration", "registration_dir" ) 
# dicom convertion executable 
dicomconvert     = config.get("registration", "dicomconvert" ) 
# subsample volume executable 
subsamplevolume  = config.get(   "registration", "subsamplevolume"  ) 
subsamplefactorX = config.getint("registration", "subsamplefactorX" ) 
subsamplefactorY = config.getint("registration", "subsamplefactorY" ) 
subsamplefactorZ = config.getint("registration", "subsamplefactorZ" ) 
# registration executable
registration_exe = config.get("registration", "registration_exe" ) 
# set mesh transform command
config.set("registration","meshXformCmd","python %s/%s %s %%s" % \
                     (registration_dir,transform_script, meshfilename ) )
 
#before transfering all files copy the unregistered mesh file 
# the power file, and dicom data  to the working directory
os.system('mkdir -p %s/reg_data' % jobid)
if(os.system('cp %s %s/reg_data' % (unregistered_mesh ,jobid))):
   raise "\n\n    error with cp %s %s/reg_data"  % (unregistered_mesh,jobid)
# to avoid copying svn meta data
os.system('mkdir -p %s/%s ' % (jobid ,
                         filter(len,interop_dicom_dir.split("/")).pop()  ))
os.system('mkdir -p %s/%s ' % (jobid ,
                         filter(len, preop_dicom_dir.split("/")).pop()   ))
if(os.system('cp %s/* %s/%s ' % (interop_dicom_dir,jobid ,
                       filter(len,interop_dicom_dir.split("/")).pop()  ))):
   raise "\n\n    error with cp %s/i* %s/%s"  %  \
   (interop_dicom_dir,jobid,filter(len,interop_dicom_dir.split("/")).pop() )
if(os.system('cp %s/* %s/%s ' % ( preop_dicom_dir ,jobid ,
                       filter(len, preop_dicom_dir.split("/")).pop()   ))):
   raise "\n\n    error with cp %s/i* %s/%s"  %  \
   (preop_dicom_dir  ,jobid,filter(len, preop_dicom_dir.split("/")).pop()  )
# write control file to the work directory
inifile=open("%s/control.ini" % jobid ,"w")
config.write(inifile)
inifile.close
inifile.flush() # ensure the entire file is written before continuing

# get the registration host information 
[reghost,locreghost,regendian,regwork] = \
         utilities.hostinformation(config,"registration","reghost","regwork")

# code execution strategy - cd to working directory ; 
#                           convert both dicoms to mha/raw/fld ;
#                           subsample the volume data; 
#                           save dimension info of interop data 
#                           to control file ; create a backup mesh ;
#                           run registration code on mha files
REGexec=["cd %s/%s; %s/%s -d %s; %s/%s -d %s -c control.ini;   \
          %s/%s %s.mha %s_%d%d%d.mha %d %d %d; \
          %s/%s %s.mha %s_%d%d%d.mha %d %d %d; \
          cd reg_data ; python %s/%s %s ; \
          python %s/%s ../control.ini;" % \
           (regwork,jobid,
            registration_dir,dicomconvert,
            filter(len,preop_dicom_dir.split("/")).pop()  ,
            registration_dir,dicomconvert,
            filter(len,interop_dicom_dir.split("/")).pop(), 
            registration_dir, subsamplevolume,
            filter(len,preop_dicom_dir.split("/")).pop()  ,
            filter(len,preop_dicom_dir.split("/")).pop()  ,
            subsamplefactorX,subsamplefactorY,subsamplefactorZ,
            subsamplefactorX,subsamplefactorY,subsamplefactorZ, 
            registration_dir, subsamplevolume,
            filter(len,interop_dicom_dir.split("/")).pop()  ,
            filter(len,interop_dicom_dir.split("/")).pop()  ,
            subsamplefactorX,subsamplefactorY,subsamplefactorZ,
            subsamplefactorX,subsamplefactorY,subsamplefactorZ, 
            registration_dir,transform_script, meshfilename ,
            registration_dir,createbackupmesh ),
         "cd %s/%s/reg_data; time \
          %s/%s ../%s_%d%d%d.mha ../%s_%d%d%d.mha ../control.ini > reg.log " % \
           (regwork,jobid, registration_dir, registration_exe, 
                    filter(len,preop_dicom_dir.split("/")).pop()  ,
                    subsamplefactorX,subsamplefactorY,subsamplefactorZ,
                    filter(len,interop_dicom_dir.split("/")).pop(),
                    subsamplefactorX,subsamplefactorY,subsamplefactorZ)  ]

#data transfer script 
if( not locreghost ):  # do not transfer if are doing local computations
   # open the data transfer script that MAY be used registration
   registeredfem=open("registeredfem.txt","w")
   registeredfem.write('get %s/input_compact %s \n' % (jobid,jobid) )
   registeredfem.write('get %s/%s/files/control.ini %s/%s/files \n' % (jobid,jobid,jobid,jobid) )
   # close registration script
   registeredfem.close
   registeredfem.flush() # ensure entire file written before continuing

# write avs command line file
clifile=open("%s/mesh_reg.cli" % jobid,"w")
clifile.write('net_clear \n' )
clifile.write('present "Network Editor"         \n' )
clifile.write('#net_read opens the network file \n' )
clifile.write('net_read mesh_reg.net            \n' )
clifile.write('script -close                    \n' )
clifile.write('parm_set "read ucd.user.9":"read file" %s/%s/reg_data/%s.inp \n'  % \
                             ( regwork,jobid, meshfilename.split(".")[0] ) )
clifile.write('parm_set "read ucd.user.17":"read file" %s/%s/%slaser_pos.inp \n'  % \
                     ( regwork,jobid, interop_dicom_dir.split("/").pop() ) )
clifile.write('parm_set "read field.user.0":"Read Field Browser" %s/%s/%s%s0000.fld \n' % \
  (regwork,jobid,filter(len,interop_dicom_dir.split("/")).pop(),regendian) )
# close registration script
clifile.close
clifile.flush() # ensure entire file written before continuing

# write the registration network
clifile=open("%s/mesh_reg.net" % jobid,"w")
clifile.write(utilities.registration_network())
clifile.close # close registration script
clifile.flush() # ensure entire file written before continuing

#determine file transfer method and registration host code execution method
if(locreghost): #local registration
   xfermethDAT='ls' 
   xferdestDAT=''.join([regwork,'/',jobid,'/input_compact']) 
   execREG= [REGexec[0]  , " %s & " % REGexec[1]  ]
else:   #remote registration
   xfermethDAT='sftp -b registeredfem.txt' ; 
   xferdestDAT=''.join([username,"@",reghost,":",regwork])
   #need to setup environment 
   if(reghost.split(".")[0] == "maverick" or reghost.split(".")[0] == "mav1"):
      execREG=["ssh %s@%s ' . /etc/profile ; %s '    " %  \
                                            (username , reghost , REGexec[0]),
               "ssh %s@%s ' . /etc/profile ; %s ' &  " %  \
                                            (username , reghost , REGexec[1]) ]
   else: # default execution method
      execREG=["ssh %s@%s ' %s '    "% (username , reghost , REGexec[0]),
               "ssh %s@%s ' %s ' &  "% (username , reghost , REGexec[1]) ]
print "registration execution commands: %s \n" % execREG

regxfer='%s %s ' % (xfermethDAT,xferdestDAT)
print "registered mesh transfer method: %s \n" % regxfer
   
#get the visualization host information
[vishost,locvishost,visendian,viswork] = \
         utilities.hostinformation(config,"output","vishost","viswork")
#setup anatomy data on vishost
anatomydata="%s/%s/%s%s0000.fld"%(regwork,jobid,
                     filter(len,interop_dicom_dir.split("/")).pop() , visendian)
if(locvishost): #local visualization
   anatxfermeth='cp %s %s' % (anatomydata,viswork)
else:   #remote visualization
   anatxfermeth='scp %s %s@%s:%s &' % (anatomydata,username,vishost,viswork)
print "anatomy transfer commands: %s \n" % anatxfermeth

# sit idle until user inputs ready to continue
utilities.pause_until_ready()

#run the registration code on reg host then bring data back when ready
if( os.system(execREG[0]) ): # wait for this to finish
    raise "\n\n     error setting up registration data\n"
os.system(anatxfermeth) # do NOT wait for this to finish
os.system(execREG[1])   # do NOT wait for this to finish

status = "i" ; mesh_not_retrieved=1 ; 
while(status != "c" or mesh_not_retrieved ):
   try: 
      if(mesh_not_retrieved):
         print "\ncannot continue until registered mesh has been retrieved"
      print "(g)et the registered mesh (c)ontinue (q)uit"
      status = sys.stdin.read()
      if(  status == "g"):
         print "\ngetting the registered mesh"
         # use error code from system call to determine if can continue
         mesh_not_retrieved = os.system(regxfer) 
      elif(status == "q"):
         print "\n\n    quitting"
         sys.exit(0) 
   except KeyboardInterrupt: 
      pass

print "\n\n    Registration Complete......"
