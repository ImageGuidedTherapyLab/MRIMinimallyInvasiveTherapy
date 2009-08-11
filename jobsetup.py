"""various functions to setup jobs for dddas scripts
 
      distribute_nopts

"""
import os  # make os a global variable in this modules namespace
import sys
import utilities
import copy
import ConfigParser 

Delimiter=";"

def setuplitt(config):
   print """param_study:
     default setup, may vary the intial guess of parameters 
     into several batch jobs 
   """
   joblist=[]
   workdir = config.get( "compexec" , "workdir" ) 
   jobid   = config.get( "compexec" , "jobid" ) 
   #default lists
   listk_0  = map(float,config.get("thermal_conductivity","k_0_healthy").split(Delimiter))
   listk_1  = map(float,config.get("thermal_conductivity","k_1").split(Delimiter))
   listk_2  = map(float,config.get("thermal_conductivity","k_2").split(Delimiter))
   listk_3  = map(float,config.get("thermal_conductivity","k_3").split(Delimiter))
   listw_0  = map(float,config.get("perfusion","w_0_healthy").split(Delimiter))
   listw_n  = map(float,config.get("perfusion","w_n").split(Delimiter))
   listw_i  = map(float,config.get("perfusion","w_i").split(Delimiter))
   listw_d  = map(float,config.get("perfusion","w_d").split(Delimiter))
   listw_2  = map(float,config.get("perfusion","w_2").split(Delimiter))
   listw_ni = map(float,config.get("perfusion","w_ni").split(Delimiter))
   listw_id = map(float,config.get("perfusion","w_id").split(Delimiter))
   listx_0  = map(float,config.get("probe","x_0").split(Delimiter))
   listy_0  = map(float,config.get("probe","y_0").split(Delimiter))
   listz_0  = map(float,config.get("probe","z_0").split(Delimiter))
   listu_flux     = map(float,config.get("bc","u_flux").split(Delimiter)) 
   listnewton_coeff = map(float,config.get("bc","newton_coeff").split(Delimiter)) 
   listoptimize_w_0 = config.get("qoi_0","optimize_w_0").split(Delimiter)
   listoptimize_k_0 = config.get("qoi_0","optimize_k_0").split(Delimiter)
   listoptimize_k_1 = config.get("qoi_0","optimize_k_1").split(Delimiter)
   listoptimize_k_2 = config.get("qoi_0","optimize_k_2").split(Delimiter)
   listoptimize_k_3 = config.get("qoi_0","optimize_k_3").split(Delimiter)
   listoptimize_pow = config.get("qoi_0","optimize_pow").split(Delimiter)
   listoptimize_mu_a= config.get("qoi_0","optimize_mu_a").split(Delimiter)
   listoptimize_mu_s= config.get("qoi_0","optimize_mu_s").split(Delimiter)
   listk_1_ub = map(float,config.get("thermal_conductivity","k_1_ub").split(Delimiter))
   listk_0_ub = map(float,config.get("thermal_conductivity","k_0_ub").split(Delimiter))
   listpde = config.get("method","pde").split(Delimiter)
   listqoi = config.get("method","qoi").split(Delimiter)
   listw_0_field = config.get("field","w_0_field").split(Delimiter)
   listk_0_field = config.get("field","k_0_field").split(Delimiter)
   #assumed of the form time_window_list "[18,234];[13,58];[12,321];[0,1]"
   time_window_list= config.get("qoi_0","init_time_window")
   listtime_window = map(utilities.ExtractListData,time_window_list.split(Delimiter))

   # method list
   listmethod  = config.get("compexec","method").split(Delimiter)

   # laser params default from control file
   listmu_a   = map(float,config.get("optical","mu_a_healthy").split(Delimiter))
   listmu_s   = map(float,config.get("optical","mu_s_healthy").split(Delimiter))
   listanfact = map(float,config.get("optical","anfact").split(Delimiter))

   # get mesh data and # proc info
   # repeat mesh for speedup plot
   listmeshfile =config.get("compexec","meshdata" ).split(Delimiter) 
   listnumproc = map(int,config.get("compexec","numproc").split(Delimiter))

   # error check
   if( len(listnumproc) != len(listmeshfile) ): 
       raise IndexError("\n\n  len(listnumproc) != len(listmeshfile)  Need number of procs for each mesh... repeat mesh for speedup plot")

   # power data
   listpowerfile=config.get("compexec","powerdata").split(Delimiter)  

   # echo params
   print "listk_0       "      , listk_0       
   print "listk_1       "      , listk_1       
   print "listk_2       "      , listk_2       
   print "listk_3       "      , listk_3       
   print "listw_0       "      , listw_0       
   print "listw_n       "      , listw_n       
   print "listw_i       "      , listw_i       
   print "listw_d       "      , listw_d       
   print "listw_2       "      , listw_2       
   print "listw_ni      "      , listw_ni      
   print "listw_id      "      , listw_id      
   print "listx_0       "      , listx_0       
   print "listy_0       "      , listy_0       
   print "listz_0       "      , listz_0       
   print "listmu_a      "      , listmu_a      
   print "listmu_s      "      , listmu_s      
   print "listanfact    "      , listanfact    
   print "listu_flux    "      , listu_flux    
   print "listnewton_coeff"    , listnewton_coeff
   print "listmethod"          , listmethod
   print "listoptimize_w_0"    , listoptimize_w_0
   print "listoptimize_k_0"    , listoptimize_k_0
   print "listoptimize_k_1"    , listoptimize_k_1
   print "listoptimize_k_2"    , listoptimize_k_2
   print "listoptimize_k_3"    , listoptimize_k_3
   print "listoptimize_pow"    , listoptimize_pow
   print "listoptimize_mu_a"   , listoptimize_mu_a
   print "listoptimize_mu_s"   , listoptimize_mu_s
   print "listk_0_ub"          , listk_0_ub
   print "listk_1_ub"          , listk_1_ub
   print "listmeshfile"        , listmeshfile
   print "listpowerfile"       , listpowerfile
   print "listpde"             , listpde
   print "listqoi"             , listqoi
   print "listw_0_field"       , listw_0_field 
   print "listk_0_field"       , listk_0_field 
   print "listtime_window"     , listtime_window 
   #use list comprehension to build entire set of parameter list
   paramlist =[ (meshdata,powerdata,
                 k_0,k_1,k_2,k_3,w_0,w_n,w_i,w_d,w_2,w_ni,w_id,x_0,y_0,z_0,
                 mu_s,mu_a,anfact,u_flux,newton_coeff, method, optimize_w_0,
                 optimize_k_0, optimize_k_1, optimize_k_2, optimize_k_3, 
                 optimize_pow, optimize_mu_a, optimize_mu_s, 
                 w_0_field, k_0_field,k_0_ub,k_1_ub,objective,time_window,pde)
                    for meshdata         in listmeshfile   
                    for powerdata        in listpowerfile   
                    for k_0              in listk_0 
                    for k_1              in listk_1 
                    for k_2              in listk_2 
                    for k_3              in listk_3 
                    for w_0              in listw_0 
                    for w_n              in listw_n 
                    for w_i              in listw_i 
                    for w_d              in listw_d 
                    for w_2              in listw_2 
                    for w_ni             in listw_ni
                    for w_id             in listw_id
                    for x_0              in listx_0 
                    for y_0              in listy_0 
                    for z_0              in listz_0 
                    for mu_s             in listmu_s
                    for mu_a             in listmu_a 
                    for anfact           in listanfact 
                    for u_flux           in listu_flux 
                    for newton_coeff     in listnewton_coeff
                    for method           in listmethod      
                    for optimize_w_0     in listoptimize_w_0  
                    for optimize_k_0     in listoptimize_k_0  
                    for optimize_k_1     in listoptimize_k_1  
                    for optimize_k_2     in listoptimize_k_2  
                    for optimize_k_3     in listoptimize_k_3  
                    for optimize_pow     in listoptimize_pow  
                    for optimize_mu_a    in listoptimize_mu_a  
                    for optimize_mu_s    in listoptimize_mu_s  
                    for w_0_field        in listw_0_field  
                    for k_0_field        in listk_0_field  
                    for k_0_ub           in listk_0_ub        
                    for k_1_ub           in listk_1_ub        
                    for objective        in listqoi        
                    for time_window      in listtime_window 
                    for pde              in listpde       ]

   # don't run too many
   utilities.verify_job_submission(len(paramlist),40)

   #create batch job list
   id = 0 
   #create script to process the optimization iteration function value outputs
   fcnvalfile=open("%s/printvalue.txt" % jobid ,"w")

   for (meshdata,powerdata,
          k_0,k_1,k_2,k_3,w_0,w_n,w_i,w_d,w_2,w_ni,w_id,x_0,y_0,z_0,
          mu_s,mu_a,anfact,u_flux,newton_coeff, method, optimize_w_0,
          optimize_k_0, optimize_k_1, optimize_k_2, optimize_k_3, 
          optimize_pow, optimize_mu_a, optimize_mu_s, w_0_field, 
          k_0_field,k_0_ub,k_1_ub,objective,time_window,pde) in paramlist:
      # create directory hierarchy to store files
      namejob= "%s%02d" % (jobid,id)
      # create directories
      utilities.create_directories(jobid,namejob)
      fcnvalfile.write("echo %s using pde %s \n" % (namejob,pde))
      fcnvalfile.write("grep 'Function value' %s/out.o* | head -1 \n" % namejob)
      fcnvalfile.write("grep  func_eval       %s/out.o* | tail -1 \n" % namejob)
      fcnvalfile.write("grep 'Function value' %s/out.o* | tail -1 \n" % namejob)
      fcnvalfile.write("grep 'Function value' %s/out.o* | cut  -d ' ' -f3,6 | cut -d ',' -f'1 2' --output-delimiter=' ' > %s/iter.dat \n" % (namejob,namejob))
      fcnvalfile.write("grep 'Function evaluation' %s/out.o* | cut  -d ':' -f3 | awk '{print NR $0}' > %s/func.dat \n" % (namejob,namejob))
      print "creating job %s using pde %s " % (namejob,pde)
      id = id + 1 # update counter

      # functions are call by reference need to deepcopy
      cntrlfile = copy.deepcopy(config) 
      cntrlfile.set(     "field"      ,"k_0_field"       ,  k_0_field        )
      cntrlfile.set(     "field"      ,"w_0_field"       ,  w_0_field        )
      cntrlfile.set("thermal_conductivity" ,"k_0_healthy", "%f" % k_0        )
      cntrlfile.set("thermal_conductivity" ,"k_0_ub"     , "%f" % k_0_ub     )
      cntrlfile.set("thermal_conductivity" ,"k_1"        , "%f" % k_1        )
      cntrlfile.set("thermal_conductivity" ,"k_1_ub"     , "%f" % k_1_ub     )
      cntrlfile.set("thermal_conductivity" ,"k_2"        , "%f" % k_2        )
      cntrlfile.set("thermal_conductivity" ,"k_3"        , "%f" % k_3        )
      cntrlfile.set("perfusion"     ,"w_0_healthy"       , "%f" % w_0        )
      cntrlfile.set("perfusion"     ,"w_n"             , "%f" % w_n        )
      cntrlfile.set("perfusion"     ,"w_i"             , "%f" % w_i        )
      cntrlfile.set("perfusion"     ,"w_d"             , "%f" % w_d        )
      cntrlfile.set("perfusion"     ,"w_2"             , "%f" % w_2        )
      cntrlfile.set("perfusion"     ,"w_ni"            , "%f" % w_ni       )
      cntrlfile.set("perfusion"     ,"w_id"            , "%f" % w_id       )
      cntrlfile.set("probe"         ,"x_0"             , "%f" % x_0        ) 
      cntrlfile.set("probe"         ,"y_0"             , "%f" % y_0        )
      cntrlfile.set("probe"         ,"z_0"             , "%f" % z_0        )
      cntrlfile.set("optical"       ,"mu_a_healthy"    , "%f" % mu_a       )
      cntrlfile.set("optical"       ,"mu_s_healthy"    , "%f" % mu_s       )
      cntrlfile.set("optical"       ,"anfact"          , "%f" % anfact     )
      cntrlfile.set("bc"            ,"u_flux"          , "%f" % u_flux     )
      cntrlfile.set("bc"            ,"newton_coeff"    , "%f" % newton_coeff )
      # set mesh file 
      try:
         assert os.path.exists( meshdata )
         if( meshdata[0] == "/" ): # full path given
           cntrlfile.set("compexec","meshdata",meshdata )
         else:
           cntrlfile.set("compexec","meshdata","%s/%s" % (os.getcwd(),meshdata))
      except AssertionError:
         raise IOError("\n\n    error with %s " % meshdata)
      # number of proc in listnumproc shoud coincide with mesh number
      numproclistid = listmeshfile.index(meshdata)
      numproc = listnumproc[numproclistid]
      cntrlfile.set("compexec" ,"numproc", numproc )
      # set power file 
      try:
        assert os.path.exists( powerdata )
        if( powerdata[0] == "/" ): # full path given
          cntrlfile.set("compexec","powerdata",powerdata )
        else:
          cntrlfile.set("compexec","powerdata","%s/%s" %(os.getcwd(),powerdata))
      except AssertionError:
        timePowerList = map(utilities.ExtractListData,  powerdata.split("@"))
        # write default power file 
        noptsteps   = config.getint( "qoi_0" ,"noptsteps" )
        numideal    = config.getint( "qoi_0" ,"numideal"  )
        MRTI_ntime  = config.getint( "mrti"  ,"ntime"     )
        Maxtime = max(time_window[1] + numideal * (noptsteps-1),MRTI_ntime)
        utilities.write_power_file(Maxtime,timePowerList,
                                   "%s/%s/files/power.ini" % (jobid,namejob)) 
        cntrlfile.remove_option("compexec","powerdata") # clean up
      cntrlfile.set("qoi_0","optimize_w_0"     ,     optimize_w_0      )
      cntrlfile.set("qoi_0","optimize_k_0"     ,     optimize_k_0      )
      cntrlfile.set("qoi_0","optimize_k_1"     ,     optimize_k_1      )
      cntrlfile.set("qoi_0","optimize_k_2"     ,     optimize_k_2      )
      cntrlfile.set("qoi_0","optimize_k_3"     ,     optimize_k_3      )
      cntrlfile.set("qoi_0","optimize_pow"     ,     optimize_pow      )
      cntrlfile.set("qoi_0","optimize_mu_a"    ,     optimize_mu_a     )
      cntrlfile.set("qoi_0","optimize_mu_s"    ,     optimize_mu_s     )
      cntrlfile.set("qoi_0","ideal_nzero_init" , "%d" % time_window[0] )
      cntrlfile.set("qoi_0","ideal_ntime_init" , "%d" % time_window[1] )
      cntrlfile.set("method",     "qoi"        ,       objective       )
      cntrlfile.set("method",     "pde"        ,          pde          )
      joblist.append( [namejob, numproc, " ", cntrlfile, method]   )
   #close script
   fcnvalfile.close; fcnvalfile.flush()
   #make executable
   fcnvalfile=os.chmod("%s/printvalue.txt" % jobid ,00770)

   return joblist

def setupkalman(iniFile):

   # get the thermal image info
   MRTI_nslice = iniFile.getint("mrti","nslice")
   nzeroList  = map(int,iniFile.get("mrti","nzero").split(Delimiter))
   ntimeList  = map(int,iniFile.get("mrti","ntime").split(Delimiter))

   # get model parameters
   x_0 = iniFile.getfloat("source_laser","x_0")
   y_0 = iniFile.getfloat("source_laser","y_0")
   z_0 = iniFile.getfloat("source_laser","z_0")
   listk_0  = map(float,iniFile.get("thermal_conductivity","k_0").split(Delimiter))
   listw_0  = map(float,iniFile.get("perfusion","w_0").split(Delimiter))
   meascovList   = map(float,iniFile.get("kalman","meascov").split(Delimiter))
   statecovList  = map(float,iniFile.get("kalman","statecov").split(Delimiter))
   
   # variations in ROI 
   # assumed of the form [ix,nx,iy,ny];[ix,nx,iy,ny];...
   roiList = map(utilities.ExtractListData,
                   iniFile.get("kalman","roi").split(Delimiter))

   # variations in solution methods
   linalgList = map(int,iniFile.get("kalman","method").split(Delimiter))
   solverList = iniFile.get("compexec","solver").split(Delimiter)
   
   # create list of command line params
   cmdLineParams =[ (nzero,ntime,meascov,statecov,roi,linalg,solver) 
                    for nzero     in   nzeroList
                    for ntime     in   ntimeList
                    for meascov   in   meascovList
                    for statecov  in   statecovList
                    for roi       in   roiList
                    for linalg    in   linalgList
                    for solver    in   solverList
                  ]
   listcmdLine=[]
   for (nzero,ntime,meascov,statecov,roi,linalg,solver)  in cmdLineParams:
      cmdLineOpt  = "-nslice %d -nzero %d -ntime %d" % \
                         (MRTI_nslice,nzero,ntime)
      cmdLineOpt  = cmdLineOpt  + \
                    " -X_0 %f -Y_0 %f -Z_0 %f -meascov %f -statecov %f" % \
                          (x_0,y_0,z_0,meascov,statecov)
      cmdLineOpt  = cmdLineOpt  + \
                    " -ix %d -nx %d -iy %d -ny %d" % \
                           (roi[0],roi[1],roi[2],roi[3])
      cmdLineOpt  = cmdLineOpt  + \
                    " -method %d -solver %s" % (linalg,solver)
      listcmdLine.append( cmdLineOpt  )

   # variations in execution options
   listmethod  = iniFile.get("compexec","method").split(Delimiter)
   listnumproc = map(int,iniFile.get("compexec","numproc").split(Delimiter))

   paramlist =[ (numproc,cmdLine_options,method) 
                    for numproc          in listnumproc 
                    for cmdLine_options  in listcmdLine 
                    for method           in listmethod 
              ]
   id = 0 
   for (numproc,cmdLine_options,method) in paramlist:
      namejob= "%s%02d" % (jobid,id)
      # create directories
      utilities.create_directories(jobid,namejob)
      joblist.append( [namejob, numproc, cmdLine_options, iniFile, method]   )
      id = id + 1 

   # don't run too many
   utilities.verify_job_submission(len(joblist),40)

   return joblist 

#debugging
if __name__ == "__main__":
   print __doc__
   import utilities
   cntrldflt = ConfigParser.ConfigParser()
   cntrldflt.add_section("qoi_0")
   cntrldflt.add_section("compexec")
   cntrldflt.set( "qoi_0" ,"noptsteps"       ,"4") 
   cntrldflt.set( "qoi_0" ,"noffset"         ,"1 ") 
   cntrldflt.set( "qoi_0" ,"numideal"        ,"5") 
   cntrldflt.set( "qoi_0" ,"ideal_nzero_init","5") 
   cntrldflt.set( "qoi_0" ,"ideal_ntime_init","10") 
   cntrldflt.set( "compexec" ,"meshdata","tags") 
   cntrldflt.set( "compexec" ,"powerdata","tags") 
   jobid = "test"
   workdir = os.getcwd()
   print cntrldflt.get( "qoi_0" ,"noptsteps" )
   JOBS = distribute_nopts(cntrldflt)
   print cntrldflt.get( "qoi_0" ,"noptsteps" )
