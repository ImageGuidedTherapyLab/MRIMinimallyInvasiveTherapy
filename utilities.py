"""utilities functions for dddas scripts

      control_defaults - set control file defaults

      waitandsend - checks that the file exists and the full data
                    has been written to disk before sending 
                       the variables: location,destination,mainfile
                          are set in the scope before the function is called
"""
import ConfigParser
import os      # make os a global variable in this modules namespace
import socket  # make socket a global variable in this modules namespace
from math import pi,atan

def control_defaults(controlfile):
   cntrldflt = ConfigParser.ConfigParser()
   #default control file values 
   cntrldflt.add_section("mrti")
   cntrldflt.add_section("compexec")
   cntrldflt.add_section("output")
   cntrldflt.add_section("avs")
   cntrldflt.add_section("source_laser")
   cntrldflt.add_section("registration")
   cntrldflt.add_section("thermal_conductivity")
   cntrldflt.add_section("perfusion")
   cntrldflt.add_section("bc")
   cntrldflt.add_section("field")
   cntrldflt.add_section("qoi_0")
   cntrldflt.add_section("visualase")
   cntrldflt.set("field"       ,  "w_0_field"      ,      "False"      )
   cntrldflt.set("field"       ,  "k_0_field"      ,      "False"      )
   cntrldflt.set("visualase"   ,"visualase_file"   ,"mlat.dat"         )
   cntrldflt.set("visualase"   ,"visualase_dir"    ,    "/tmp"         )
   cntrldflt.set("visualase"   ,  "override"       ,      "False"      )
   cntrldflt.set("visualase"   ,  "override_power" ,       "0.0"       )
   cntrldflt.set("compexec"    ,    "method"       ,      "opt"        )
   cntrldflt.set("compexec"    ,     "comphost"    , socket.gethostname().split(".")[0]) 
   cntrldflt.set("mrti"        ,"img_sizefile"     ,"mrtiimg.size"     )    
   cntrldflt.set("mrti"        ,"nzero"            ,        "0"        )    
   cntrldflt.set("avs"         ,"mrtivis_sizefile" ,"mrtivis.size"     )    
   cntrldflt.set("avs"         ,"femvis_sizefile"  ,"ucdvis.size"      )     
   cntrldflt.set("output"      ,"compwritemri"     ,"mrivis/"          )     
   cntrldflt.set("output"      ,"compwritefem"     ,"femvis/"          )     
   cntrldflt.set("output"      ,"fem_filename"     ,"fem"              )
   cntrldflt.set("output"      ,"computefinalobjective" ,   "false"    )
   cntrldflt.set("output"      ,      "femavs"     ,"False"            )
   cntrldflt.set("output"      ,      "femgmv"     ,"False"            )
   cntrldflt.set("output"      ,      "femraw"     ,"False"            )
   cntrldflt.set("output"      ,      "mrtiavs"    ,"False"            )
   cntrldflt.set("output"      ,      "mrtirawiv"  ,"False"            )
   cntrldflt.set("output"      ,      "mriavs"     ,"False"            )
   cntrldflt.set("output"      ,      "mrirawiv"   ,"False"            )
   cntrldflt.set("qoi_0"       ,   "noptsteps"     ,       "1"         )
   cntrldflt.set("qoi_0"       ,   "optimize_w_0"  ,   "False"         )
   cntrldflt.set("qoi_0"       ,   "optimize_k_0"  ,   "False"         )
   cntrldflt.set("qoi_0"       ,   "optimize_k_1"  ,   "False"         )
   cntrldflt.set("qoi_0"       ,   "optimize_k_2"  ,   "False"         )
   cntrldflt.set("qoi_0"       ,   "optimize_k_3"  ,   "False"         )
   cntrldflt.set("qoi_0"       ,   "optimize_pow"  ,   "False"         )
   cntrldflt.set("qoi_0"       ,   "optimize_mu_a" ,   "False"         )
   cntrldflt.set("qoi_0"       ,   "optimize_mu_s" ,   "False"         )
   cntrldflt.set("thermal_conductivity" , "k_0_lb"   ,   "0.40e0 "       )
   cntrldflt.set("thermal_conductivity" , "k_1_lb"   ,   "0.0e0  "       )
   cntrldflt.set("thermal_conductivity" , "k_2_lb"   ,   "1.0e-2 "       )
   cntrldflt.set("thermal_conductivity" , "k_3_lb"   ,   "305.0e0"       )
   cntrldflt.set("thermal_conductivity" , "k_0_ub"   ,   "0.72e0 "       )
   cntrldflt.set("thermal_conductivity" , "k_1_ub"   ,   "0.33e0 "       )
   cntrldflt.set("thermal_conductivity" , "k_2_ub"   ,   "10.0e0 "       )
   cntrldflt.set("thermal_conductivity" , "k_3_ub"   ,   "325.0e0"       )
   cntrldflt.set("perfusion"            , "w_0_lb"   ,   " 6.0e0 "       )
   cntrldflt.set("perfusion"            , "w_0_ub"   ,   " 50.0e0"       )
   cntrldflt.set("perfusion"            , "w_n_lb"   ,   " 4.0e0 "       )
   cntrldflt.set("perfusion"            , "w_n_ub"   ,   "20.0e0 "       )
   cntrldflt.set("perfusion"            , "w_i_lb"   ,   "20.0e0 "       )
   cntrldflt.set("perfusion"            , "w_i_ub"   ,   " 1.0e2 "       )
   cntrldflt.set("perfusion"            , "w_d_lb"   ,   " 0.0e0 "       )
   cntrldflt.set("perfusion"            , "w_d_ub"   ,   " 4.0e0 "       )
   cntrldflt.set("perfusion"            , "w_2_lb"   ,   " 1.0e-1"       )
   cntrldflt.set("perfusion"            , "w_2_ub"   ,   " 10.0e0"       )
   cntrldflt.set("perfusion"            , "W_NID_LB" ,   "300.0e0"       )
   cntrldflt.set("perfusion"            , "W_NID_MD" ,   "315.0e0"       )
   cntrldflt.set("perfusion"            , "W_NID_UB" ,   "335.0e0"       )
   cntrldflt.set("source_laser",     "anfact"      ,   "0.71e0"        )
   cntrldflt.set("source_laser",     "anfact_lb"   ,   "0.6e0"         )
   cntrldflt.set("source_laser",     "anfact_ub"   ,   "1.0e0"         )
   cntrldflt.set("source_laser",     "POW_LB"      ,   "-0.1e0"        )
   cntrldflt.set("source_laser",     "POW_UB"      ,   "15.0e0"        )
   cntrldflt.set("source_laser",     "mu_a_LB"     ,   "0.30e+2"       )
   cntrldflt.set("source_laser",     "mu_a_UB"     ,   "0.60e+2"       )
   cntrldflt.set("source_laser",     "mu_s_LB"     ,   "2.0e+2"        )
   cntrldflt.set("source_laser",     "mu_s_UB"     ,   "2820.0e+2"     )
   cntrldflt.set(    "bc"      ,  "newton_coeff_LB",     "0.1e0"       )
   cntrldflt.set(    "bc"      ,  "newton_coeff_UB",   "112.0e0"       )
   cntrldflt.set(    "bc"      ,   "u_flux_LB"     ,     "0.1e0"       )
   cntrldflt.set(    "bc"      ,   "u_flux_UB"     ,   "112.0e0"       )
   cntrldflt.set("registration",  "registration"   ,    "False"        )
   cntrldflt.set("registration","interop_dicom_dir",       ""          )
   cntrldflt.set("registration","dicomconvert"     ,"DicomSeriesConvert" )
   cntrldflt.set("registration","registration_exe" ,"RigidRegistration"  )
   cntrldflt.set("registration","transform_script" ,"LBIE_2_hp3d.py 1 1" )
   cntrldflt.set("registration","createbackupmesh" ,"createboxmesh.py"   )
   cntrldflt.set("registration", "subsamplevolume" ,"SubsampleNoFilter") 
   cntrldflt.set("registration", "subsamplefactorX",      "2"          ) 
   cntrldflt.set("registration", "subsamplefactorY",      "2"          ) 
   cntrldflt.set("registration", "subsamplefactorZ",      "1"          ) 

   #set some defaults based on an initial read of the control file
   #scratchini = ConfigParser.ConfigParser()
   #scratchini.readfp( open("%s" % controlfile , "r") )
   #comphost = scratchini.get("compexec","comphost")
   #reghost  = scratchini.get("registration","reghost")
   #if(reghost.split(".")[0]=="mav1" or reghost.split(".")[0]=="maverick"):
   #  registration_dir  = "/home/utexas/iv/fuentes/registration"
   #elif(reghost.split(".")[0]=="DIPWS019"):
   #  registration_dir  = "/home/dfuentes/DDDAS/trunk/registration"
   #else:
   #  registration_dir  = "/org/groups/oden/fuentes/DDDAS/trunk/registration"
   #if(comphost.split(".")[0]=="lonestar"):
   #  comp_rank_begin=4
   #elif(comphost.split(".")[0]=="shamu"):
   #  comp_rank_begin=8
   #else:
   #  comp_rank_begin=1

   #cntrldflt.set("registration","registration_dir" ,   registration_dir    )
   #cntrldflt.set("compexec"    ,"data_rank_begin"  ,         "1"           )
   #cntrldflt.set("compexec"    ,"comp_rank_begin"  ,"%d" % comp_rank_begin )

   return cntrldflt

#check for a previous run of the same jobid to avoid overwrite
def checkprevious(JobID):
   if(os.path.isdir(JobID)):
     print ""
     print "jobid %s ALREADY EXISTS!!!!" % JobID
     print ""
     numjobid = 0
     while(os.path.isdir("%s%d" % (JobID,numjobid)) ):
        numjobid = numjobid + 1
     newID = "%s%d" %(JobID,numjobid)
     print ""
     print "jobid is now %s !!!!" % newID
     print ""
   else: 
     newID = JobID
   os.mkdir(newID)
   return newID 

# sit idle until user inputs ready to continue
def pause_until_ready():
   import sys
   status = "i"
   while(status != "c"):
      try: 
         print "(c)ontinue (q)uit"
         status = sys.stdin.read()
         if(status == "q"):
           print "\n\n    quitting"
           sys.exit(0) 
      except KeyboardInterrupt: 
         pass
   print "\ncontinuing......"
   return

# don't run too many
def verify_job_submission(numJob,Nmax):
   if( numJob  > Nmax ) : 
     print "\n\n    %d jobs > %d job \n\n" % (numJob,Nmax)
     # sit idle until user inputs ready to continue
     status = "i"
     while(status != "c"):
        try: 
           print "(c)ontinue (q)uit"
           status = sys.stdin.read()
           if(status == "q"):
             print "\n\n    quitting"
             sys.exit(0) 
        except KeyboardInterrupt: 
           pass
     print "\ncontinuing......"
   return

# create directory hierarchy to store files
def create_directories(JobID,JobName):
   os.system('mkdir -p %s/%s/files  ' %  (JobID,JobName)  )
   os.system('mkdir -p %s/%s/mrivis ' %  (JobID,JobName)  )
   os.system('mkdir -p %s/%s/femvis ' %  (JobID,JobName)  )
   return

def location(Host,HostID):
   if(Host.split(".")[0] == socket.gethostname().split(".")[0]):
      print "%s: %s is local  \n" % (Host,HostID)
      return True
   else:
      print "%s: %s is remote \n" % (Host,HostID)
      return False

def endian(Host):
   # define a list of big endian machines
   # most machine used are little endian so will default to little endian
   # this will probably cause a bug later though :-(
   table = ["maverick","mav1"]
   if( Host.split(".")[0] in table ):
      return "bigend" 
   else:
      return "litend" 

def hostinformation(IniFile,Section,Hosttype,Workdir):
   print "hello" , Section, Hosttype
   host    = IniFile.get(Section,Hosttype) # host name
   print "host", host
   lochost = location(host,Hosttype) # lochost = false ==> remote host
   hostendian = endian(host)   # endianess of host
   if(lochost): # for local host set work directory to $CWD
      curdir = os.getcwd();  hostwork=curdir;
      IniFile.set(Section,Workdir,hostwork) # to pass to jobsetup
   else:  # get user input work directory on remote host
      hostwork=IniFile.get(Section,Workdir)
   return [host,lochost,hostendian,hostwork]

def variable_range(IniFile,ParamID,mapping="linear",LB="lb",UB="ub"):
   Section = ParamID[0]
   Keyword = ParamID[1]
   beta_lb=IniFile.getfloat(Section,"%s_%s" % (Keyword,LB) )
   beta_ub=IniFile.getfloat(Section,"%s_%s" % (Keyword,UB) )
   try:
     ndiv = IniFile.getint(Section,"%s_nd" % Keyword)
   except:
     ndiv = 2
   #skip endpoints
   if(mapping=="linear"):
     def mappingfunction(i):
         return ( beta_lb*(1.0-float(i)/float(ndiv+1)) + 
                  beta_ub*     float(i)/float(ndiv+1)   )
   elif(mapping=="atan"):
     def mappingfunction(i):
         return (beta_ub+beta_lb)/2.0  + ( beta_ub - beta_lb ) / pi * atan(15.0*(float(i)- float(ndiv+1)/2.0 ))
   else:
     print "unknown mapping function" 
     raise
   return map(mappingfunction,range(1,ndiv+1))

#assume ListString of the form ="[18,234,13];[58];[873,12,321,235];[1,0]"
def ExtractListData(x):
   entries = filter(len,filter(len,x.split("["))[0].split("]"))[0].split(",")
   try:
     dataList = map(int, entries)
   except ValueError:
     dataList = map(float, entries)
   return dataList

# the data structure is setup as follows and ASSUMES equispace time 
#     distances, IDEAL_DT   \forall i
#  
#                                           *NOTE* closed at beginning
# time = 0    ---------                        BUT open   at end
#                 |                                |
#                 |                               \|/
#              Power[0]    power between time = [0,1)  is Power[0]
#                 |            
#                 |
# time = 1    ---------
#                 |
#                 |          
#              Power[1]    power between time = [1,2)  is Power[1]   
#                 |         
#                 |
# time = 2    ---------
#                 |
#                 |        
#              Power[2]    power between time = [2,3)  is Power[2]
#                 |       
#                 |
# time = 3    ---------
#           .
#           .
#           .
#           .
def write_power_file(Maxtime,timePowerList,Filename):
   # error check sorting
   assert timePowerList[0] == sorted(timePowerList[0])
   # extend if necessary
   if ( timePowerList[0][len(timePowerList[0])-1] < Maxtime ):
      timePowerList[0].append(Maxtime)
      timePowerList[1].append(0.0)
   intervalID = 0
   powerFile=open(Filename ,"w")
   for iBound in timePowerList[0]:
     while (intervalID < iBound):
      powerFile.write("%d %f\n"% \
                (intervalID, timePowerList[1][timePowerList[0].index(iBound)]) )
      intervalID = intervalID + 1 
   powerFile.close; powerFile.flush()
   return 
   

#debugging
if __name__ == "__main__":
   x = "[18,234,13];[58];[873,12,321,235];[1,0]"
   y = map(ExtractListData,x.split(";"))
   print y

def singlehistory():
    return """# Gnuplot script file for plotting optimization history
#   usage:
#         sed "s/namejob/<jobid>/g" history.plt | gnuplot
#                          -or-
#         just edit namejob
set   autoscale    # scale axes automatically
unset log          # remove any log-scaling
unset label        # remove any previous labels
set xtic auto      # set xtics automatically
set ytic auto      # set ytics automatically
set xlabel "iteration #"
#set xr [0:30]
#set yr [0:100]
set ylabel "function value [qoi units]"
set key right

set ylabel "function value [qoi units]"
set title "optimization iteration history"
set term x11 1
plot "namejob/iter.dat" using 1:2 title "namejob iteration history" w lp lw 2

set ylabel "function value [qoi units]"
set title "optimization function history"
set term x11 2
plot "namejob/func.dat" using 1:2 title "namejob evaluation history" w lp lw 2

set ylabel "W_0 [kg/s/m^3]"
set title "perfusion history"
set term x11 3
plot "namejob/W_0.dat" using 1:2 title "namejob evaluation history" w lp lw 2

set ylabel "k_0 [J/s/m/K]"
set title "thermal conductivity history"
set term x11 4
plot "namejob/K_0.dat" using 1:2 title "namejob evaluation history" w lp lw 2

set ylabel "k_1 [J/s/m/K]"
set title "thermal conductivity history"
set term x11 5
plot "namejob/K_1.dat" using 1:3 title "namejob evaluation history" w lp lw 2

set ylabel "k_2 [1/K]"
set title "thermal conductivity history"
set term x11 6
plot "namejob/K_2.dat" using 1:3 title "namejob evaluation history" w lp lw 2

set ylabel "k_3 [K]"
set title "thermal conductivity history"
set term x11 7
plot "namejob/K_3.dat" using 1:3 title "namejob evaluation history" w lp lw 2

set ylabel "w_n [kg/s/m^3]"
set title "perfusion history"
set term x11 8
plot "namejob/W_N.dat" using 1:3 title "namejob evaluation history" w lp lw 2

set ylabel "w_i [kg/s/m^3]"
set title "perfusion history"
set term x11 9
plot "namejob/W_I.dat" using 1:3 title "namejob evaluation history" w lp lw 2

set ylabel "w_d [kg/s/m^3]"
set title "perfusion history"
set term x11 10
plot "namejob/W_D.dat" using 1:3 title "namejob evaluation history" w lp lw 2

set ylabel "w_2 [kg/s/m^3]"
set title "perfusion history"
set term x11 11
plot "namejob/W_2.dat" using 1:3 title "namejob evaluation history" w lp lw 2

set ylabel "x_0 [m]"
set title "position history"
set term x11 12
plot "namejob/X_0.dat" using 1:3 title "namejob evaluation history" w lp lw 2

set ylabel "y_0 [m]"
set title "position history"
set term x11 13
plot "namejob/Y_0.dat" using 1:3 title "namejob evaluation history" w lp lw 2

set ylabel "z_0 [m]"
set title "position history"
set term x11 14
plot "namejob/Z_0.dat" using 1:3 title "namejob evaluation history" w lp lw 2

set ylabel "mu_a [m^-1]"
set title "absorbtion history"
set term x11 15
plot "namejob/MU_A.dat" using 1:3 title "namejob evaluation history" w lp lw 2

set ylabel "mu_s [m^-1]"
set title "scattering history"
set term x11 16
plot "namejob/MU_S.dat" using 1:3 title "namejob evaluation history" w lp lw 2

set ylabel "w_ni [1/K]"
set title "perfusion history"
set term x11 17
plot "namejob/W_NID.dat" using 1:3 title "namejob evaluation history" w lp lw 2

set ylabel "w_id [1/K]"
set title "perfusion history"
set term x11 18
plot "namejob/W_NID.dat" using 1:5 title "namejob evaluation history" w lp lw 2


#set term x11 3
#set title "optimization history"
#plot  \\
#"namejob/iter.dat" using 1:2 title "namejob iterations" w lp lw 2 , \\
#"namejob/func.dat" using 1:2 title "namejob function evaluation" w lp lw 2
"""

def plotpower():
   return """# system module
import sys
import os

if (len(sys.argv) < 2 ):
    raise "[usage]: python power.py <powerfile> \\n"

#import data
powerfile=open(sys.argv[1],"r")
data=[(float(line.split()[0]),float(line.split()[1])) for line in powerfile]
powerfile.close()
##################process power file and write a tmp file to view in gnuplot
gnuplotfile=open("power.plt" ,"w")
gnuplotfile.write('# Gnuplot script file for power file plotting \\n' )
gnuplotfile.write('set   autoscale    # scale axes automatically   \\n' )
gnuplotfile.write('unset log          # remove any log-scaling     \\n' )
gnuplotfile.write('unset label        # remove any previous labels \\n' )
gnuplotfile.write('set xtic auto      # set xtics automatically    \\n' )
gnuplotfile.write('set ytic auto      # set ytics automatically    \\n' )
gnuplotfile.write('#set xr [0:30] \\n#set yr [0:100]\\n' )
gnuplotfile.write('set ylabel "Power [Watts]" \\n' )
gnuplotfile.write('set key right\\n' )
gnuplotfile.write('set title "Power History" \\n' )

##################file extension determines file format
fileext=sys.argv[1].split(".")[1]
if(fileext=="dat"): #single data file
 ##################process power file and write a tmp file to view in gnuplot
 t_prev = data[0][0]
 tmppower=open("/tmp/tmp.pow","w")
 for (time,pow) in data[1:]:
     tmppower.write("%f %f %f\\n" % ((time+t_prev)/2.,pow,time-t_prev))
     t_prev = time
 tmppower.close; tmppower.flush() # close file
 #finish gnuplot file
 gnuplotfile.write('set xlabel "time [s]" \\n' )
 gnuplotfile.write('plot "/tmp/tmp.pow" using 1:2:3 title "power" w boxes \\n')
 gnuplotfile.close; gnuplotfile.flush()
 gnuplotfile.close; gnuplotfile.flush()

elif(fileext=="his"): # history data file
 numsteps=int(sys.argv[2])
 #finish gnuplot file
 gnuplotfile.write('set xlabel "idstep" \\n' )
 gnuplotfile.write('plot \\\\\\n')
 for i in range(len(data)/numsteps):
    tmppower=open("/tmp/tmp%d.pow" % i,"w")
    for (idstep,pow) in data[i*numsteps:(i+1)*numsteps-1]:
        tmppower.write("%f %f\\\\\\n" % (idstep,pow))
    tmppower.close; tmppower.flush() # close file
    gnuplotfile.write('"/tmp/tmp%d.pow" using 1:2 title "iter=%03d" w lp'%(i,i))
    if(i < len(data)/numsteps -1 ):
      gnuplotfile.write(', \\\\')
    gnuplotfile.write('\\n')
 gnuplotfile.close; gnuplotfile.flush()
else:
 raise "unknown file"
os.system("gnuplot -persist power.plt - ") """

      
def registration_network():
   return """#!/usr/bin/avs -network
version 5.5 (50.86 i686 ogl)
module "read field.user.0" -xy 182,41 -ex $Path/avs_library/mongo
module "generate colormap.user.1" -xy 0,72
module "color range.user.2" -xy 0,182 -ex $Path/avs_library/mongo
module "orthogonal slicer.user.3" -xy 78,122 -ex $Path/avs_library/mongo
module "field to mesh.user.4" -xy 28,322 -ex $Path/avs_library/mongo
module downsize.user.5 -xy 218,122 -ex $Path/avs_library/mongo
module isosurface.user.6 -xy 198,252 -ex $Path/avs_library/tile
module "set view.user.7" -xy 0,2 -ex $Path/avs_library/mongo
module "print field.user.8" -xy 328,82 -ex $Path/avs_library/mongo
module "read ucd.user.9" -xy 448,32 -ex $Path/avs_library/ucd_multm
module "ucd to geom.user.10" -xy 411,351 -ex $Path/avs_library/ucd_multm
module "ucd legend.user.11" -xy 563,311 -ex $Path/avs_library/ucd_legend
module "generate colormap.user.12" -xy 418,192
module "ucd crop.user.13" -xy 398,152 -ex $Path/avs_library/ucd_multm
module "geometry viewer.user.14" -xy 245,393
module "ucd probe.user.15" -xy 628,142 -ex $Path/avs_library/ucd_multm
module crop.user.16 -xy 45,251 -ex $Path/avs_library/mongo
module "read ucd.user.17" -xy 628,42 -ex $Path/avs_library/ucd_multm
port_connect "read field.user.0":0 "print field.user.8":0
port_connect "read field.user.0":0 "orthogonal slicer.user.3":0
port_connect "read field.user.0":0 downsize.user.5:0
port_connect "generate colormap.user.1":0 isosurface.user.6:2
port_connect "generate colormap.user.1":0 "color range.user.2":1
port_connect "color range.user.2":0 "field to mesh.user.4":1
port_connect "orthogonal slicer.user.3":0 "color range.user.2":0
port_connect "orthogonal slicer.user.3":0 crop.user.16:0
port_connect "field to mesh.user.4":0 "geometry viewer.user.14":0
port_connect downsize.user.5:0 isosurface.user.6:1
port_connect downsize.user.5:0 isosurface.user.6:0
port_connect "read ucd.user.9":0 "ucd probe.user.15":0
port_connect "read ucd.user.9":0 "ucd legend.user.11":0
port_connect "read ucd.user.9":0 "ucd crop.user.13":0
port_connect "ucd to geom.user.10":0 "geometry viewer.user.14":0
port_connect "ucd legend.user.11":0 "ucd to geom.user.10":1
port_connect "generate colormap.user.12":0 "ucd legend.user.11":1
port_connect "ucd crop.user.13":0 "ucd to geom.user.10":0
port_connect "ucd crop.user.13":1 "geometry viewer.user.14":0
port_connect "ucd probe.user.15":0 "geometry viewer.user.14":0
port_connect crop.user.16:0 "field to mesh.user.4":0
parm_set "generate colormap.user.1":colormap "{ \
{0.666666687,0,0.0599999987,0}\
{0.664052308,0,0.0549999997,1.53787023e-05}\
{0.661437929,0,0.0599999987,6.15148092e-05}\
{0.65882355,0,0.0649999976,0.000138408315}\
{0.656209171,0,0.0700000003,0.000246059237}\
{0.653594792,0,0.075000003,0.000384467538}\
{0.650980413,0,0.0799999982,0.000553633261}\
{0.648366034,0,0.0799999982,0.00075355632}\
{0.645751655,0,0.0850000009,0.000984236947}\
{0.643137276,0,0.0900000036,0.00124567491}\
{0.640522897,0,0.0900000036,0.00153787015}\
{0.637908518,0,0.100000001,0.0018608229}\
{0.635294139,0,0.104999997,0.00221453304}\
{0.63267976,0,0.104999997,0.00259900047}\
{0.630065382,0,0.107500002,0.00301422528}\
{0.627451003,0,0.109999999,0.00346020772}\
{0.624836624,0,0.115000002,0.00393694779}\
{0.622222245,0,0.115000002,0.00444444502}\
{0.619607866,0,0.119999997,0.00498269964}\
{0.616993487,0,0.119999997,0.00555171119}\
{0.614379108,0,0.125,0.0061514806}\
{0.611764729,0,0.129999995,0.0067820074}\
{0.60915035,0,0.129999995,0.0074432916}\
{0.606535971,0,0.135000005,0.00813533273}\
{0.603921592,0,0.137500003,0.00885813218}\
{0.601307213,0,0.140000001,0.00961168855}\
{0.598692834,0,0.150000006,0.0103960019}\
{0.596078455,0,0.155000001,0.011211073}\
{0.593464077,0,0.155000001,0.0120569011}\
{0.590849698,0,0.159999996,0.0129334871}\
{0.588235319,0,0.162499994,0.0138408309}\
{0.58562094,0,0.165000007,0.0147789316}\
{0.583006561,0,0.165000007,0.0157477912}\
{0.580392182,0,0.170000002,0.0167474076}\
{0.577777803,0,0.174999997,0.0177777801}\
{0.575163424,0,0.180000007,0.0188389104}\
{0.572549045,0,0.185000002,0.0199307986}\
{0.569934666,0,0.189999998,0.0210534427}\
{0.567320287,0,0.194999993,0.0222068448}\
{0.564705908,0,0.200000003,0.0233910047}\
{0.562091529,0,0.204999998,0.0246059224}\
{0.55947715,0,0.209999993,0.025851598}\
{0.556862772,0,0.215000004,0.0271280296}\
{0.554248393,0,0.215000004,0.0284352191}\
{0.551634014,0,0.219999999,0.0297731664}\
{0.549019635,0,0.224999994,0.0311418697}\
{0.546405256,0,0.230000004,0.0325413309}\
{0.543790877,0,0.234999999,0.0339715518}\
{0.541176498,0,0.239999995,0.0354325287}\
{0.538562119,0,0.25,0.0369242616}\
{0.53594774,0,0.254999995,0.0384467542}\
{0.533333361,0,0.25999999,0.0400000028}\
{0.530718982,0,0.264999986,0.0415840074}\
{0.528104603,0,0.270000011,0.0431987718}\
{0.525490224,0,0.270000011,0.0448442921}\
{0.522875845,0,0.275000006,0.0465205684}\
{0.520261467,0,0.280000001,0.0482276045}\
{0.517647088,0,0.284999996,0.0499654002}\
{0.515032709,0,0.284999996,0.0517339483}\
{0.51241833,0,0.289999992,0.0535332561}\
{0.509803951,0,0.289999992,0.0553633235}\
{0.507189572,0,0.305000007,0.0572241433}\
{0.504575193,0,0.305000007,0.0591157265}\
{0.501960814,0,0.310000002,0.061038062}\
{0.499346405,0,0.314999998,0.0629911646}\
{0.496732026,0,0.314999998,0.0649750158}\
{0.494117647,0,0.319999993,0.0669896305}\
{0.491503268,0,0.319999993,0.0690349936}\
{0.48888889,0,0.324999988,0.0711111203}\
{0.486274511,0,0.330000013,0.0732180029}\
{0.483660132,0,0.330000013,0.0753556415}\
{0.481045753,0,0.335000008,0.0775240362}\
{0.478431374,0,0.340000004,0.0797231942}\
{0.475816995,0,0.340000004,0.0819531009}\
{0.473202616,0,0.349999994,0.0842137709}\
{0.470588237,0,0.354999989,0.086505197}\
{0.467973858,0,0.360000014,0.088827379}\
{0.465359479,0,0.362500012,0.0911803246}\
{0.4627451,0,0.36500001,0.0935640186}\
{0.460130721,0,0.36500001,0.0959784761}\
{0.457516342,0,0.370000005,0.0984236896}\
{0.454901963,0,0.375,0.100899659}\
{0.452287585,0,0.375,0.103406392}\
{0.449673206,0,0.379999995,0.105943874}\
{0.447058827,0,0.38499999,0.108512118}\
{0.444444448,0,0.389999986,0.111111119}\
{0.441830069,0,0.389999986,0.113740876}\
{0.43921569,0,0.394999981,0.116401389}\
{0.436601311,0,0.400000006,0.119092666}\
{0.433986932,0,0.405000001,0.121814691}\
{0.431372553,0,0.405000001,0.124567479}\
{0.428758174,0,0.409999996,0.127351031}\
{0.426143795,0,0.414999992,0.130165324}\
{0.423529416,0,0.414999992,0.133010387}\
{0.420915037,0,0.419999987,0.135886207}\
{0.418300658,0,0.425000012,0.138792783}\
{0.41568628,0,0.430000007,0.141730115}\
{0.413071901,0,0.430000007,0.144698203}\
{0.410457522,0,0.435000002,0.147697046}\
{0.407843143,0,0.439999998,0.150726646}\
{0.405228764,0,0.449999988,0.153787017}\
{0.402614385,0,0.455000013,0.156878129}\
{0.400000006,0,0.455000013,0.160000011}\
{0.397385627,0,0.460000008,0.163152635}\
{0.394771248,0,0.460000008,0.16633603}\
{0.392156869,0,0.469999999,0.16955018}\
{0.38954249,0,0.469999999,0.172795087}\
{0.386928111,0,0.474999994,0.17607075}\
{0.384313732,0,0.479999989,0.179377168}\
{0.381699353,0,0.485000014,0.182714343}\
{0.379084975,0,0.49000001,0.186082274}\
{0.376470596,0,0.5,0.189480975}\
{0.373856217,0,0.502499998,0.192910418}\
{0.371241838,0,0.504999995,0.196370631}\
{0.368627459,0,0.504999995,0.199861601}\
{0.36601308,0,0.514999986,0.203383312}\
{0.363398701,0,0.514999986,0.206935793}\
{0.360784322,0,0.524999976,0.210519031}\
{0.358169943,0,0.529999971,0.214133024}\
{0.355555564,0,0.535000026,0.217777774}\
{0.352941185,0,0.535000026,0.221453294}\
{0.350326806,0,0.540000021,0.225159556}\
{0.347712427,0,0.550000012,0.228896573}\
{0.345098048,0,0.555000007,0.232664362}\
{0.34248367,0,0.560000002,0.236462906}\
{0.339869291,0,0.564999998,0.240292192}\
{0.337254912,0,0.569999993,0.244152248}\
{0.334640533,0,0.574999988,0.24804306}\
{0.332026124,0,0.579999983,0.251964658}\
{0.329411745,0,0.586250007,0.255916983}\
{0.326797366,0,0.592499971,0.259900063}\
{0.324182987,0,0.598749995,0.2639139}\
{0.321568608,0,0.605000019,0.267958522}\
{0.318954229,0,0.610000014,0.27203387}\
{0.31633985,0,0.61500001,0.276139975}\
{0.313725471,0,0.620000005,0.280276835}\
{0.311111093,0,0.625,0.284444481}\
{0.308496714,0,0.629999995,0.288642853}\
{0.305882335,0,0.632499993,0.292872012}\
{0.303267956,0,0.63499999,0.297131896}\
{0.300653577,0,0.63499999,0.301422566}\
{0.298039198,0,0.639999986,0.305743963}\
{0.295424819,0,0.649999976,0.310096145}\
{0.29281044,0,0.654999971,0.314479083}\
{0.290196061,0,0.660000026,0.318892777}\
{0.287581682,0,0.660000026,0.323337197}\
{0.284967303,0,0.670000017,0.327812403}\
{0.282352924,0,0.675000012,0.332318366}\
{0.279738545,0,0.680000007,0.336855084}\
{0.277124166,0,0.682500005,0.341422558}\
{0.274509788,0,0.685000002,0.346020788}\
{0.271895409,0,0.699999988,0.350649774}\
{0.26928103,0,0.702499986,0.355309516}\
{0.266666651,0,0.704999983,0.360000014}\
{0.264052272,0,0.709999979,0.364721298}\
{0.261437893,0,0.714999974,0.369473308}\
{0.258823514,0,0.720000029,0.374256074}\
{0.256209135,0,0.725000024,0.379069626}\
{0.253594756,0,0.730000019,0.383913904}\
{0.250980377,0,0.735000014,0.388788968}\
{0.248365998,0,0.74000001,0.393694758}\
{0.245751619,0,0.746250033,0.398631334}\
{0.24313724,0,0.752499998,0.403598636}\
{0.240522861,0,0.758749962,0.408596724}\
{0.237908483,0,0.764999986,0.413625568}\
{0.235294104,0,0.769999981,0.418685138}\
{0.232679725,0,0.772499979,0.423775494}\
{0.230065346,0,0.774999976,0.428896606}\
{0.227450967,0,0.779999971,0.434048474}\
{0.224836588,0,0.785000026,0.439231098}\
{0.222222209,0,0.790000021,0.444444478}\
{0.21960783,0,0.790000021,0.449688613}\
{0.216993451,0,0.800000012,0.454963505}\
{0.214379072,0,0.805000007,0.460269153}\
{0.211764693,0,0.810000002,0.465605557}\
{0.209150314,0,0.814999998,0.470972717}\
{0.206535935,0,0.819999993,0.476370662}\
{0.203921556,0,0.824999988,0.481799334}\
{0.201307178,0,0.827499986,0.487258762}\
{0.198692799,0,0.829999983,0.492748976}\
{0.19607842,0,0.834999979,0.498269916}\
{0.193464041,0,0.841250002,0.503821611}\
{0.190849662,0,0.847499967,0.509404123}\
{0.188235283,0,0.85374999,0.515017331}\
{0.185620904,0,0.860000014,0.520661294}\
{0.183006525,0,0.86500001,0.526336074}\
{0.180392146,0,0.870000005,0.53204155}\
{0.177777767,0,0.872500002,0.537777781}\
{0.175163388,0,0.875,0.543544829}\
{0.172549009,0,0.879999995,0.549342573}\
{0.16993463,0,0.88499999,0.555171132}\
{0.167320251,0,0.88499999,0.561030388}\
{0.164705873,0,0.889999986,0.566920459}\
{0.162091494,0,0.899999976,0.572841227}\
{0.159477115,0,0.899999976,0.57879281}\
{0.156862736,0,0.902499974,0.58477509}\
{0.154248357,0,0.904999971,0.590788186}\
{0.151633978,0,0.915000021,0.596832037}\
{0.149019599,0,0.917500019,0.602906585}\
{0.14640522,0,0.920000017,0.609011948}\
{0.143790841,0,0.930000007,0.615148067}\
{0.141176462,0,0.930000007,0.621314883}\
{0.138562083,0,0.930000007,0.627512515}\
{0.135947704,0,0.935000002,0.633740902}\
{0.133333325,0,0.939999998,0.640000045}\
{0.130718946,0,0.949999988,0.646289885}\
{0.128104568,0,0.954999983,0.65261054}\
{0.125490189,0,0.959999979,0.658961952}\
{0.12287581,0,0.959999979,0.665344119}\
{0.120261431,0,0.964999974,0.671757042}\
{0.117647052,0,0.970000029,0.678200722}\
{0.115032673,0,0.975000024,0.684675157}\
{0.112418294,0,0.980000019,0.691180348}\
{0.109803915,0,0.985000014,0.697716296}\
{0.107189536,0,0.987500012,0.704282999}\
{0.104575157,0,0.99000001,0.710880458}\
{0.101960778,0,1,0.717508674}\
{0.0993463993,0,1,0.724167645}\
{0.0967320204,0,1,0.730857372}\
{0.0941176414,0,1,0.737577856}\
{0.0915032625,0,1,0.744329095}\
{0.0888888836,0,1,0.75111115}\
{0.0862745047,0,1,0.757923901}\
{0.0836601257,0,1,0.764767408}\
{0.0810457468,0,1,0.771641672}\
{0.0784313679,0,1,0.778546751}\
{0.0758169889,0,1,0.785482526}\
{0.07320261,0,1,0.792449057}\
{0.0705882311,0,1,0.799446404}\
{0.0679738522,0,1,0.806474447}\
{0.0653594732,0,1,0.813533247}\
{0.0627450943,0,1,0.820622861}\
{0.0601307154,0,1,0.827743173}\
{0.0575163364,0,1,0.8348943}\
{0.0549019575,0,1,0.842076123}\
{0.0522875786,0,1,0.849288762}\
{0.0496731997,0,1,0.856532097}\
{0.0470588207,0,1,0.863806248}\
{0.0444444418,0,1,0.871111095}\
{0.0418300629,0,1,0.878446758}\
{0.0392156839,0,1,0.885813177}\
{0.036601305,0,1,0.893210292}\
{0.0339869261,0,1,0.900638223}\
{0.0313725471,0,1,0.90809691}\
{0.0287581682,0,1,0.915586293}\
{0.0261437893,0,1,0.923106492}\
{0.0235294104,0,1,0.930657446}\
{0.0209150314,0,1,0.938239157}\
{0.0183006525,0,1,0.945851624}\
{0.0156862736,0,1,0.953494787}\
{0.0130718946,0,1,0.961168766}\
{0.0104575157,0,1,0.968873501}\
{0.00784313679,0,1,0.976608992}\
{0.00522875786,0,1,0.984375238}\
{0.00261437893,0,1,0.992172241}\
{0,0,1,1}\
}"
parm_set "orthogonal slicer.user.3":"slice plane" 12 -range 0 12
parm_set isosurface.user.6:level 253.8120117 -range 0 253.8120117
parm_set "ucd crop.user.13":choice space
parm_set "ucd crop.user.13":"Do Crop" true
parm_set "ucd probe.user.15":x 0.0531144999
parm_set "ucd probe.user.15":y 0.02439289913
parm_set "ucd probe.user.15":z 0.006012199912
parm_set "ucd probe.user.15":"Probe Type" "Probe Coordinates"
parm_set "ucd probe.user.15":type Crosshair
parm_set crop.user.16:"max x" 255 -range 0 255
parm_set crop.user.16:"max y" 255 -range 0 255
geom_set_scene -scene "geometry viewer.user.14"
geom_set_camera_name "Camera 1"
geom_resize_camera -view "Camera 1" 811 649
#
# State for view: Camera 1
#
geom_set_position -view "Camera 1" 0 0 -12
geom_set_view_modes -depth_cue 0 -view "Camera 1"
geom_set_view_modes -polygonal_spheres 0 -view "Camera 1"
geom_set_view_modes -stereo 0 -view "Camera 1"
geom_set_view_modes -head_tracking 1 -view "Camera 1"
geom_set_view_modes -z_buffer 1 -view "Camera 1"
geom_set_camera_params -view "Camera 1" -front -88 -back 112
geom_set_depth_cue_params "Camera 1" -scale 0.1
#
# Light state
#
geom_set_light -light 1 -type bi-directional -state 1
geom_set_light -type ambient -state 1
#
# State for object: top
#
geom_set_cur_cli_obj top
geom_set_matrix   -mat \
                  29.6377    -0.127345    -1.87304    0 \
                    0.0905008    -29.4974    3.43741    0 \
                    -1.87518    -3.43624    -29.438    0 \
                    -0.628068    0.479982    0.2846    1 
geom_set_position   -1.19209e-07 1.19209e-07 -7.45058e-08
geom_set_obj_window -0.0117481 0.0553091 -0.0146682 0.0446862 -0.0196865 0.0397555
#
# State for object: "field mesh"
#
geom_set_cur_cli_obj -push
geom_set_name_context "field to mesh.user.4"
geom_create_obj "field mesh" -mod "field to mesh.user.4"
geom_set_trans_mode parent
geom_set_name_context
geom_set_cur_cli_obj -pop
#
# State for object: "text for probe"
#
geom_set_cur_cli_obj -push
geom_set_name_context probe.user.7
geom_create_obj "text for probe" -mod probe.user.7
geom_set_matrix   -mat \
                  0.871929    0.031311    -0.488629    0 \
                    -0.109141    0.985274    -0.131619    0 \
                    0.477312    0.168092    0.862507    0 \
                    0    0    0    1 
geom_set_name_context
geom_set_cur_cli_obj -pop
#
# State for object: probe
#
geom_set_cur_cli_obj -push
geom_set_name_context probe.user.7
geom_create_obj probe -mod probe.user.7
geom_set_trans_mode redirect
geom_set_position   -0.00350188 0.0297662 0
geom_set_name_context
geom_set_cur_cli_obj -pop
#
# State for object: "ucd geom"
#
geom_set_cur_cli_obj -push
geom_set_name_context "ucd to geom.user.10"
geom_create_obj "ucd geom" -mod "ucd to geom.user.10"
geom_set_trans_mode parent
geom_set_obj_window -0.0117481 0.0553091 -0.0146682 0.0446862 -0.0196865 0.0397555
geom_set_name_context
geom_set_cur_cli_obj -pop
#
# State for object: "probe geom"
#
geom_set_cur_cli_obj -push
geom_set_name_context "ucd probe.user.15"
geom_create_obj "probe geom" -mod "ucd probe.user.15"
geom_set_trans_mode redirect
geom_set_name_context
geom_set_cur_cli_obj -pop
#
# State for object: "ucd probe labels"
#
geom_set_cur_cli_obj -push
geom_set_name_context "ucd probe.user.15"
geom_create_obj "ucd probe labels" -mod "ucd probe.user.15"
geom_set_trans_mode parent
geom_set_obj_window -0.00616 0.049721 -0.009722 0.03974 -0.014733 0.034802
geom_set_name_context
geom_set_cur_cli_obj -pop
#
# State for object: probe
#
geom_set_cur_cli_obj -push
geom_set_name_context "ucd probe.user.15"
geom_create_obj probe -mod "ucd probe.user.15"
geom_set_trans_mode redirect
geom_set_name_context
geom_set_cur_cli_obj -pop
#
# State for object: "crop probe"
#
geom_set_cur_cli_obj -push
geom_set_name_context "ucd crop.user.13"
geom_create_obj "crop probe" -mod "ucd crop.user.13"
geom_set_trans_mode redirect
geom_set_properties -amb 0.300 -diff 0.700 -spec 0.000 -exp 50.000 \
 -trans 0.500 -spec_col 0.000 0.000 0.000
geom_set_obj_window -0.00616 0.049721 -0.009722 0.03974 -0.014733 0.034802
geom_set_name_context
geom_set_cur_cli_obj -pop
shell "ui" shell
 panel Application -w app_panel -p ui -xy 0,0 -wh 260,1024
  panel "Top Level Stack" -w master_stack -p Application -xy 2,99 -wh 254,880\
   -P columns integer 1
   panel "set view.user.7" -w panel -p "Top Level Stack" -xy 0,64 -wh 194,86
    manipulator "set view.user.7:User" -w oneshot -p "set view.user.7" \
        -xy 7,9 -wh 58,22
    manipulator "set view.user.7:Top" -w oneshot -p "set view.user.7" \
        -xy 65,9 -wh 58,22
    manipulator "set view.user.7:Bottom" -w oneshot -p "set view.user.7" \
        -xy 124,9 -wh 58,22
    manipulator "set view.user.7:Front" -w oneshot -p "set view.user.7" \
        -xy 7,28 -wh 58,22
    manipulator "set view.user.7:Back" -w oneshot -p "set view.user.7" \
        -xy 65,28 -wh 58,22
    manipulator "set view.user.7:Right" -w oneshot -p "set view.user.7" \
        -xy 124,28 -wh 58,22
    manipulator "set view.user.7:Left" -w oneshot -p "set view.user.7" \
        -xy 7,54 -wh 58,22
    manipulator "set view.user.7:Bounds" -w toggle -p "set view.user.7" \
        -xy 65,54 -wh 58,22
    manipulator "set view.user.7:Persp" -w toggle -p "set view.user.7" \
        -xy 124,54 -wh 58,22
   panel stack.0 -w stack -p "Top Level Stack" -xy 0,109 -wh 253,749\
   -P columns integer 1\
   -P title string "MRI Image"
    panel "read field.user.0" -w panel -p stack.0 -xy 0,176 -wh 252,329
     manipulator "read field.user.0:Read Field Browser" -w browser -p "read field.user.0" \
         -xy 7,9 -wh 233,199
     manipulator "read field.user.0:Data Conversion" -w radio_buttons -p "read field.user.0" \
         -xy 7,206 -wh 117,45
     manipulator "read field.user.0:Read Status" -w textblock -p "read field.user.0" \
         -xy 7,252 -wh 235,67
    panel "generate colormap.user.1" -w panel -p stack.0 -xy 0,64 -wh 252,572
     manipulator "generate colormap.user.1:colormap" -w color_editor -p "generate colormap.user.1" \
         -xy 7,9 -wh 235,423
     manipulator "generate colormap.user.1:lo value" -w dial -p "generate colormap.user.1" \
         -xy 7,436 -wh 90,130
     manipulator "generate colormap.user.1:hi value" -w dial -p "generate colormap.user.1" \
         -xy 97,436 -wh 90,130
    panel downsize.user.5 -w panel -p stack.0 -xy 0,64 -wh 110,145
     manipulator downsize.user.5:downsize -w idial -p downsize.user.5 \
         -xy 7,9 -wh 90,130
    panel "orthogonal slicer.user.3" -w panel -p stack.0 -xy 0,64 -wh 225,145
     manipulator "orthogonal slicer.user.3:axis" -w radio_buttons -p "orthogonal slicer.user.3" \
         -xy 97,9 -wh 117,67
    panel isosurface.user.6 -w panel -p stack.0 -xy 0,64 -wh 225,215
     manipulator isosurface.user.6:level -w dial -p isosurface.user.6 \
         -xy 7,9 -wh 90,130
     manipulator "isosurface.user.6:optimize surf" -w toggle -p isosurface.user.6 \
         -xy 97,9 -wh 117,22
     manipulator "isosurface.user.6:optimize wire" -w toggle -p isosurface.user.6 \
         -xy 7,138 -wh 117,22
     manipulator "isosurface.user.6:generate normals" -w toggle -p isosurface.user.6 \
         -xy 7,160 -wh 117,22
     manipulator "isosurface.user.6:flip normals" -w toggle -p isosurface.user.6 \
         -xy 7,183 -wh 117,22
    panel "field to mesh.user.4" -w panel -p stack.0 -xy 0,176 -wh 225,145
     manipulator "field to mesh.user.4:Z scale" -w dial -p "field to mesh.user.4" \
         -xy 7,9 -wh 90,130
     manipulator "field to mesh.user.4:Normalize" -w toggle -p "field to mesh.user.4" \
         -xy 97,9 -wh 117,22
   panel stack.1 -w stack -p "Top Level Stack" -xy 0,109 -wh 254,771\
   -P columns integer 1\
   -P title string "FEM mesh"
    panel "read ucd.user.9" -w panel -p stack.1 -xy 0,199 -wh 252,215
     manipulator "read ucd.user.9:read file" -w browser -p "read ucd.user.9" \
         -xy 7,9 -wh 233,199
    panel "generate colormap.user.12" -w panel -p stack.1 -xy 0,64 -wh 252,572
     manipulator "generate colormap.user.12:colormap" -w color_editor -p "generate colormap.user.12" \
         -xy 7,9 -wh 235,423
     manipulator "generate colormap.user.12:lo value" -w dial -p "generate colormap.user.12" \
         -xy 7,436 -wh 90,130
     manipulator "generate colormap.user.12:hi value" -w dial -p "generate colormap.user.12" \
         -xy 97,436 -wh 90,130
    panel "ucd crop.user.13" -w panel -p stack.1 -xy 0,64 -wh 135,160
     manipulator "ucd crop.user.13:choice" -w radio_buttons -p "ucd crop.user.13" \
         -xy 7,9 -wh 117,45
     manipulator "ucd crop.user.13:side" -w radio_buttons -p "ucd crop.user.13" \
         -xy 7,67 -wh 117,45
     manipulator "ucd crop.user.13:Do Crop" -w toggle -p "ucd crop.user.13" \
         -xy 7,125 -wh 117,22
    panel "ucd to geom.user.10" -w panel -p stack.1 -xy 0,64 -wh 225,413
     manipulator "ucd to geom.user.10:Shrink" -w toggle -p "ucd to geom.user.10" \
         -xy 7,9 -wh 117,22
     manipulator "ucd to geom.user.10:Shrink Factor" -w idial -p "ucd to geom.user.10" \
         -xy 124,9 -wh 90,130\
   -P title string "Shrink Factor"
     manipulator "ucd to geom.user.10:Geometry Display Mode" -w text -p "ucd to geom.user.10" \
         -xy 7,138 -wh 176,22
     manipulator "ucd to geom.user.10:mode" -w radio_buttons -p "ucd to geom.user.10" \
         -xy 7,160 -wh 176,67
     manipulator "ucd to geom.user.10:Explode Materials" -w toggle -p "ucd to geom.user.10" \
         -xy 7,228 -wh 117,22
     manipulator "ucd to geom.user.10:Explode Factor" -w idial -p "ucd to geom.user.10" \
         -xy 124,228 -wh 90,130
     manipulator "ucd to geom.user.10:Save Geometry" -w toggle -p "ucd to geom.user.10" \
         -xy 7,355 -wh 117,22
     manipulator "ucd to geom.user.10:Color Cells" -w toggle -p "ucd to geom.user.10" -hide \
         -xy 7,380 -wh 117,22
    panel "ucd legend.user.11" -w panel -p stack.1 -xy 0,64 -wh 223,323
     manipulator "ucd legend.user.11:Node Data" -w text -p "ucd legend.user.11" \
         -xy 7,9 -wh 117,22
     manipulator "ucd legend.user.11:node data" -w radio_buttons -p "ucd legend.user.11" \
         -xy 7,28 -wh 117,22
     manipulator "ucd legend.user.11:value" -w dial -p "ucd legend.user.11" \
         -xy 124,28 -wh 90,130
     manipulator "ucd legend.user.11:lo value" -w dial -p "ucd legend.user.11" -hide \
         -xy 7,160 -wh 90,130
     manipulator "ucd legend.user.11:hi value" -w dial -p "ucd legend.user.11" -hide \
         -xy 97,160 -wh 90,130
     manipulator "ucd legend.user.11:range" -w toggle -p "ucd legend.user.11" \
         -xy 7,291 -wh 117,22
    panel "ucd probe.user.15" -w panel -p stack.1 -xy 0,64 -wh 215,494
     manipulator "ucd probe.user.15:Pick Geometry" -w toggle -p "ucd probe.user.15" \
         -xy 34,278 -wh 117,22
     manipulator "ucd probe.user.15:Label Options" -w text -p "ucd probe.user.15" \
         -xy 30,142 -wh 117,22
     manipulator "ucd probe.user.15:label nodes" -w toggle -p "ucd probe.user.15" \
         -xy 33,192 -wh 117,22
     manipulator "ucd probe.user.15:label id" -w toggle -p "ucd probe.user.15" \
         -xy 34,212 -wh 117,22
     manipulator "ucd probe.user.15:label value" -w toggle -p "ucd probe.user.15" \
         -xy 34,233 -wh 117,22
     manipulator "ucd probe.user.15:label cell" -w toggle -p "ucd probe.user.15" \
         -xy 34,252 -wh 117,22
     manipulator "ucd probe.user.15:Text Size" -w idial -p "ucd probe.user.15" \
         -xy 14,303 -wh 90,130
     manipulator "ucd probe.user.15:Text Offset" -w dial -p "ucd probe.user.15" \
         -xy 116,309 -wh 90,130
     manipulator "ucd probe.user.15:Node Data" -w text -p "ucd probe.user.15" \
         -xy 40,438 -wh 117,22
     manipulator "ucd probe.user.15:node data" -w radio_buttons -p "ucd probe.user.15" \
         -xy 37,462 -wh 117,22
     manipulator "ucd probe.user.15:type" -w radio_buttons -p "ucd probe.user.15" \
         -xy 32,68 -wh 117,67
    panel "read ucd.user.17" -w panel -p stack.1 -xy 0,199 -wh 254,214
     manipulator "read ucd.user.17:read file" -w browser -p "read ucd.user.17" \
         -xy 10,9 -wh 234,199
 panel "print field.user.8" -w panel -p ui -xy 613,73 -wh 584,796
  manipulator "print field.user.8:Min X" -w typein_integer -p "print field.user.8" \
      -xy 7,9 -wh 117,22
  manipulator "print field.user.8:Max X" -w typein_integer -p "print field.user.8" \
      -xy 7,28 -wh 117,22
  manipulator "print field.user.8:Min Y" -w typein_integer -p "print field.user.8" \
      -xy 7,54 -wh 117,22
  manipulator "print field.user.8:Max Y" -w typein_integer -p "print field.user.8" \
      -xy 7,73 -wh 117,22
  manipulator "print field.user.8:Min Z" -w typein_integer -p "print field.user.8" \
      -xy 7,99 -wh 117,22
  manipulator "print field.user.8:Max Z" -w typein_integer -p "print field.user.8" \
      -xy 7,119 -wh 117,22
  manipulator "print field.user.8:Min W" -w typein_integer -p "print field.user.8" \
      -xy 7,143 -wh 117,22
  manipulator "print field.user.8:Max W" -w typein_integer -p "print field.user.8" \
      -xy 7,167 -wh 117,22
  manipulator "print field.user.8:Max Elements" -w idial -p "print field.user.8" \
      -xy 124,167 -wh 90,130
  manipulator "print field.user.8:Display Header" -w toggle -p "print field.user.8" \
      -xy 7,297 -wh 117,22
  manipulator "print field.user.8:Display Data" -w toggle -p "print field.user.8" \
      -xy 7,316 -wh 117,22
  manipulator "print field.user.8:Output File" -w typein -p "print field.user.8" \
      -xy 7,342 -wh 235,22
  manipulator "print field.user.8:Output Browser" -w text_browser -p "print field.user.8" \
      -xy 7,362 -wh 462,367
 panel "geometry viewer.user.14!display" -w container -p ui -xy 650,3 -wh 943,847\
   -P zoom_coords string "0 0 0 0 0 <$NULL> 0 0 0 0"
  manipulator "orthogonal slicer.user.3:slice plane" -w idial -p "geometry viewer.user.14!display" \
      -xy 26,593 -wh 90,130
  manipulator "ucd probe.user.15:x" -w typein_real -p "geometry viewer.user.14!display" \
      -xy 20,758 -wh 117,22
  manipulator "ucd probe.user.15:y" -w typein_real -p "geometry viewer.user.14!display" \
      -xy 20,777 -wh 117,22
  manipulator "ucd probe.user.15:z" -w typein_real -p "geometry viewer.user.14!display" \
      -xy 20,803 -wh 117,22
  manipulator "ucd probe.user.15:Probe Type" -w text -p "geometry viewer.user.14!display" \
      -xy 20,732 -wh 117,22
  panel crop.user.16 -w panel -p "geometry viewer.user.14!display" \
      -xy 712,520 -wh 197,297
   manipulator "crop.user.16:min x" -w idial -p crop.user.16 \
       -xy 7,9 -wh 90,130
   manipulator "crop.user.16:max x" -w idial -p crop.user.16 \
       -xy 97,9 -wh 90,130
   manipulator "crop.user.16:min y" -w idial -p crop.user.16 \
       -xy 7,137 -wh 90,130
   manipulator "crop.user.16:max y" -w idial -p crop.user.16 \
       -xy 97,137 -wh 90,130
   manipulator "crop.user.16:min z" -w idial -p crop.user.16 -hide \
       -xy 7,266 -wh 90,130
   manipulator "crop.user.16:max z" -w idial -p crop.user.16 -hide \
       -xy 97,266 -wh 90,130
   manipulator "crop.user.16:size to fit" -w toggle -p crop.user.16 \
       -xy 5,271 -wh 117,22
manipulator "geometry viewer.user.14":object -w none
manipulator "geometry viewer.user.14":"Update Always" -w none
manipulator "geometry viewer.user.14":"Update Image" -w none
# End of file
"""

