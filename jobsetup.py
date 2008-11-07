"""various functions to setup jobs for dddas scripts
 
      distribute_nopts

"""
import os  # make os a global variable in this modules namespace
import sys
import utilities
import copy
import ConfigParser 

def setupjob(config):
   print """param_study:
     default setup, may vary the intial guess of parameters 
     into several batch jobs 
   """
   joblist=[]
   workdir = config.get( "compexec" , "workdir" ) 
   jobid   = config.get( "compexec" , "jobid" ) 
   #default lists
   listk_0          = [config.getfloat("thermalcond","k_0")]
   listk_1          = [config.getfloat("thermalcond","k_1")]
   listk_2          = [config.getfloat("thermalcond","k_2")]
   listk_3          = [config.getfloat("thermalcond","k_3")]
   listw_0          = [config.getfloat("perfusivity","w_0")]
   listw_n          = [config.getfloat("perfusivity","w_n")]
   listw_i          = [config.getfloat("perfusivity","w_i")]
   listw_d          = [config.getfloat("perfusivity","w_d")]
   listw_2          = [config.getfloat("perfusivity","w_2")]
   listw_ni         = [config.getfloat("perfusivity","w_ni")]
   listw_id         = [config.getfloat("perfusivity","w_id")]
   listx_0          = [config.getfloat("source_laser","x_0")]
   listy_0          = [config.getfloat("source_laser","y_0")]
   listz_0          = [config.getfloat("source_laser","z_0")]
   listg_flux       = [config.getfloat("neumann_boundary","g_flux")]
   listcoeff_cool   = [config.getfloat("cauchy_boundary","coeff_cool")]
   listmethod       = [" "]
   listoptimize_w_0 = [config.get("qoi_0","optimize_w_0")]
   listoptimize_k_0 = [config.get("qoi_0","optimize_k_0")]
   listoptimize_k_1 = [config.get("qoi_0","optimize_k_1")]
   listoptimize_k_2 = [config.get("qoi_0","optimize_k_2")]
   listoptimize_k_3 = [config.get("qoi_0","optimize_k_3")]
   listoptimize_pow = [config.get("qoi_0","optimize_pow")]
   listoptimize_mu_a= [config.get("qoi_0","optimize_mu_a")]
   listoptimize_mu_s= [config.get("qoi_0","optimize_mu_s")]
   listk_1_ub       = [config.getfloat("thermalcond","k_1_ub")]
   listpde          = [config.get("hp3d","pde")]
   listw_0_field    = [config.get("field","w_0_field")]
   listk_0_field    = [config.get("field","k_0_field")]
   listobjective    = [config.get("qoi_0","objective")]
   # laser params default from control file
   listmu_a         = [config.getfloat("source_laser","mu_a")]
   listmu_s         = [config.getfloat("source_laser","mu_s")]
   listanfact       = [config.getfloat("source_laser","anfact")]

   # built commands for parameters study
   paramstudyid = []; paramstudvar = []

   # setup monte carlo 
   # the database* lists mimic the ordering of the matlab code used 
   # to generate MC grids
   databasemu_s   = [3.0   *100.0,   # agar (as a guess ~ water *100)
                     47.0  *100.0,   # human prostate in-vitro
                     435.0 *100.0]   #  brain adult white matter 
   databasemu_a   = [  9.19*100.0,   # canine_prostate
                       0.60*100.0,   # rat prostate (in-vivo)
                       0.45*100.0,   # initial guess
                       0.30*100.0]   # human brain
   databaseanfact = [.71,       # myocardium
                     .862,      # human prostate in-vitro
                     .97]       # prostate rat tumor
   try: # monte carlo filename(s) found
      mc_filelocation=config.get("source_laser","mc_filelocation")
      # during preprocessing mc_file_name is used as a printf file format
      mc_filefmt=config.get("source_laser","mc_file_name")
   except ConfigParser.NoOptionError:
      #mc_file_name not found
      #mc_filecmd is meant to do nothing
      mc_filecmd= '""echo %%d %s/%%s/files/MC.asc > /dev/null""' % jobid   
      #error checking
      if(config.get("hp3d","pde")=="nonlinpennesmonte"  or
         config.getboolean("compexec","vary_pde")         ):
        raise "\n\n    NO MC GRID FOUND FOR MONTE CARLO RUN!!!!!!!  "
   else:
      if(config.getboolean("compexec","vary_mc_file")):
        # this picks out the files from the MC grid database that we are 
        #  interested in. vary_mc_cmds should be a python command to build
        #  the lists of interest. it should have NO SPACES i.e.
        #  config.get("compexec","vary_mc_cmds") == listmu_s=[3.*100.,47.0*100,435.0*100.];listmu_a=[15.5*100.,1.23*100.,.36*100.,.04*100.];listanfact=[.862]
        exec config.get("compexec","vary_mc_cmds")
        if(len(listmu_s) > 1):
           paramstudyid.append('mu_s=%.1e')
           paramstudvar.append('mu_s')
        if(len(listmu_a) > 1):
           paramstudyid.append('mu_a=%.1e')
           paramstudvar.append('mu_a')
        if(len(listanfact) > 1):
           paramstudyid.append('anfact=%.1e')
           paramstudvar.append('anfact')
      else: # same mc grid for each run
        listmu_a        = [config.getfloat("source_laser","mu_a")]
        listmu_s        = [config.getfloat("source_laser","mu_s")]
        listanfact      = [config.getfloat("source_laser","anfact")]
      # mc_filelocation gives the computation host the path to the mc_file
      # the path is the same for all
      config.set("source_laser","mc_filelocation","files" )
      # set the default MC filename
      config.set("source_laser","mc_file_name","MC.asc")
      # mc_filecmd copies the mc grids to the default location
      mc_filecmd= '""cp %s/%s %s/%%s/files/MC.asc""' % \
                                       (mc_filelocation,mc_filefmt,jobid)   

   #get the host name of the visualization computer
   vishost = config.get("output","vishost")
   viswork=config.get("output","viswork") # get viswork directory
   # get mesh data and power data info
   meshfile =config.get("compexec","meshdata" ) 
   powerfile=config.get("compexec","powerdata") 
   listmeshcmd      = []
   if(config.getint("compexec","vary_mesh")):
     print """
        vary_mesh == true
        setting up to use vary the mesh file in a list for the batch jobs 
           assumes that <mesh_dir> contains the mesh files of the 
                                       form mesh_dir/input_compact_%d
     """
     #get number of meshes
     nmeshes   = config.getint( "compexec" , "vary_mesh" ) 
     for i in range(nmeshes):
       listmeshcmd.append(['""cp %s/input_compact_%d %s/%%s/files/input_compact;       cp %s %s/%%s/files/power.dat""' %  (meshfile,i,jobid,powerfile,jobid),i])
     # meshcmd copies the meshes to the default location
     config.set("compexec","compfilelocation","files" )
     def procmap(i):
         return "%d" % (12*i + 20)
     numproclist   = map(procmap,range(nmeshes))
   elif(config.getint("compexec","vary_power")):
     print """
        vary_power == true
        setting up to vary the initial power file in a list for the batch jobs 
           assumes that <power_dir> contains the power files of the 
                                       form power_dir/power%d.dat
     """
     #get number of meshes
     npower   = config.getint( "compexec" , "vary_power" ) 
     for i in range(npower):
       listmeshcmd.append(['""cp %s %s/%%s/files/input_compact; cp %s/power%d.dat %s/%%s/files/power.dat""' %  (meshfile,jobid,powerfile,i,jobid),0])
     # meshcmd copies the meshes to the default location
     config.set("compexec","compfilelocation","files" )
     numproclist   = [config.get("compexec","numproc" )] 
   else: 
     # default is to use a single pre-registered mesh for each run
     print """
        vary_mesh  == false
        vary_power == false
        meshdata and powerdata should point to a file 
     """
     #the command appended to meshcmd is meant to do nothing
     listmeshcmd.append(['""echo %s/%%s %s/%%s > /dev/null""' % \
                                                     (jobid,jobid) , 0 ])
     numproclist   = [config.get("compexec","numproc" )] 
     # compfilelocation gives the computation host the path to the 
     # mesh and power file, the path is the same for all
     config.set("compexec","compfilelocation","%s/%s" % (workdir,jobid))
     #copy the mesh file and the power file to the working directory
     if(os.system('cp %s %s/input_compact' % (meshfile ,jobid))):
           raise "\nerror copying mesh file %s \n" % meshfile
     if(os.system('cp %s %s/power.dat' % (powerfile ,jobid))):
           raise "\nerror copying power file %s \n" % powerfile

   #get data used for batch submission
   if(config.getboolean("compexec","vary_k_0")):
     listk_0=utilities.variable_range(config,["thermalcond","k_0"],"atan")
     paramstudyid.append('k_0=%.1e')
     paramstudvar.append('k_0')
   if(config.getboolean("compexec","vary_k_1")):
     listk_1=utilities.variable_range(config,["thermalcond","k_1"],"atan")
     paramstudyid.append('k_1=%.1e')
     paramstudvar.append('k_1')
   if(config.getboolean("compexec","vary_k_2")):
     listk_2=utilities.variable_range(config,["thermalcond","k_2"],"atan")
     paramstudyid.append('k_2=%.1e')
     paramstudvar.append('k_2')
   if(config.getboolean("compexec","vary_k_3")):
     listk_3=utilities.variable_range(config,["thermalcond","k_3"],"atan")
     paramstudyid.append('k_3=%.1e')
     paramstudvar.append('k_3')
   if(config.getboolean("compexec","vary_w_0")):
     listw_0=utilities.variable_range(config,["perfusivity","w_0"],"atan")
     paramstudyid.append('w_0=%.1e')
     paramstudvar.append('w_0')
   if(config.getboolean("compexec","vary_w_n")):
     listw_n=utilities.variable_range(config,["perfusivity","w_n"])
   if(config.getboolean("compexec","vary_w_i")):
     listw_i=utilities.variable_range(config,["perfusivity","w_i"])
   if(config.getboolean("compexec","vary_w_d")):
     listw_d=utilities.variable_range(config,["perfusivity","w_d"])
   if(config.getboolean("compexec","vary_w_2")):
     listw_2=utilities.variable_range(config,["perfusivity","w_2"],"atan")
     paramstudyid.append('w_2=%.1e')
     paramstudvar.append('w_2')
   if(config.getboolean("compexec","vary_w_ni")):
     listw_ni=utilities.variable_range(config,["perfusivity","w_nid"],"atan",
                                                                      "lb","md")
     paramstudyid.append('w_ni=%.1e')
     paramstudvar.append('w_ni')
   if(config.getboolean("compexec","vary_w_id")):
     listw_id=utilities.variable_range(config,["perfusivity","w_nid"],"atan",
                                                                      "md","ub")
     paramstudyid.append('w_id=%.1e')
     paramstudvar.append('w_id')
   if(config.getboolean("compexec","vary_g_flux")):
     listg_flux=utilities.variable_range(config,
                                         ["neumann_boundary","g_flux"],"atan")
     paramstudyid.append('g_flux=%.1e')
     paramstudvar.append('g_flux')
   if(config.getboolean("compexec","vary_coeff_cool")):
     listcoeff_cool=utilities.variable_range(config,
                                        ["cauchy_boundary","coeff_cool"],"atan")
     paramstudyid.append('coeff_cool=%.1e')
     paramstudvar.append('coeff_cool')
   if(config.getboolean("compexec","vary_method")):
     #-tao_method tao_nm     has bug in initial release of tao-1.9
     #-tao_method tao_gpcg   requires hessian
     # the unconstrained solver -tao_method tao_lmvm,tao_cg 
     #  require that the ls stepmax be smaller than the default value
     #   so as not to take the penalty term to inf on the first step
     #listmethod=[" -tao_method tao_lmvm -tao_ls_stepmax 0.05",
     #            " -tao_method tao_lmvm -tao_ls_stepmax 0.3",
     #            " -tao_method tao_lmvm -tao_ls_stepmax 0.6",
     #            " -tao_method tao_blmvm" ]
     #listmethod=[" -tao_lmm_delta_max .1 -tao_ls_gtol .97 -tao_ls_ftol 1.0e-5",
     #            " -tao_lmm_delta_max .2 -tao_ls_gtol .97 -tao_ls_ftol 1.0e-5",
     #            " -tao_lmm_delta_max .1 -tao_ls_gtol 1.0 -tao_ls_ftol 1.0e-6",
     #            " -tao_lmm_delta_max .2 -tao_ls_gtol 1.0 -tao_ls_ftol 1.0e-6"]
     listmethod=[" -tao_lmm_delta_max 0.1 ",
                 " -tao_lmm_delta_max 1.0 "]
   if(config.getboolean("compexec","vary_optimize_w_0")):
     listoptimize_w_0=["false","true"]
     paramstudyid.append('opt_w_0=%d')
     paramstudvar.append('optimize_w_0')
   if(config.getboolean("compexec","vary_optimize_k_0")):
     listoptimize_k_0=["false","true"]
     paramstudyid.append('opt_k_0=%d')
     paramstudvar.append('optimize_k_0')
   if(config.getboolean("compexec","vary_optimize_k_1")):
     listoptimize_k_1=["false","true"]
     paramstudyid.append('opt_k_1=%d')
     paramstudvar.append('optimize_k_1')
   if(config.getboolean("compexec","vary_optimize_k_2")):
     listoptimize_k_2=["false","true"]
     paramstudyid.append('opt_k_2=%d')
     paramstudvar.append('optimize_k_2')
   if(config.getboolean("compexec","vary_optimize_k_3")):
     listoptimize_k_3=["false","true"]
     paramstudyid.append('opt_k_3=%d')
     paramstudvar.append('optimize_k_3')
   if(config.getboolean("compexec","vary_optimize_pow")):
     listoptimize_pow=["false","true"]
     paramstudyid.append('opt_pow=%d')
     paramstudvar.append('optimize_pow')
   if(config.getboolean("compexec","vary_optimize_mu_a")):
     listoptimize_mu_a=["false","true"]
     paramstudyid.append('opt_mu_a=%d')
     paramstudvar.append('optimize_mu_a')
   if(config.getboolean("compexec","vary_optimize_mu_s")):
     listoptimize_mu_s=["false","true"]
     paramstudyid.append('opt_mu_s=%d')
     paramstudvar.append('optimize_mu_s')
   if(config.getboolean("compexec","vary_k_1_ub")):
     listk_1_ub=[ 0.31e0 ,0.33e0,0.35e0,0.38e0]
     paramstudyid.append('k_1_ub=%.1e')
     paramstudvar.append('k_1_ub')
   if(config.getboolean("compexec","vary_objective")):
     listobjective = ["temp_control","dam_control_arr"]
   if(config.getboolean("compexec","vary_pde")):
     listpde       = ["nonlinpennesisolaser","nonlinpennesmonte"]
     paramstudyid.append('pde=%8s')
     paramstudvar.append('pde.split("nonlinpennes").pop()')
   if(config.getboolean("compexec","vary_w_0_field")):
     listw_0_field = ["false","true"]
     paramstudyid.append('w_0_fld=%d')
     paramstudvar.append('w_0_field')
   if(config.getboolean("compexec","vary_k_0_field")):
     listk_0_field = ["false","true"]
     paramstudyid.append('k_0_fld=%d')
     paramstudvar.append('k_0_field')

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
   print "listg_flux    "      , listg_flux    
   print "listcoeff_cool"      , listcoeff_cool
   print "listmethod"          , listmethod
   print "listoptimize_w_0"    , listoptimize_w_0
   print "listoptimize_k_0"    , listoptimize_k_0
   print "listoptimize_k_1"    , listoptimize_k_1
   print "listoptimize_k_2"    , listoptimize_k_2
   print "listoptimize_k_3"    , listoptimize_k_3
   print "listoptimize_pow"    , listoptimize_pow
   print "listoptimize_mu_a"   , listoptimize_mu_a
   print "listoptimize_mu_s"   , listoptimize_mu_s
   print "listk_1_ub"          , listk_1_ub
   print "listmeshcmd"         , listmeshcmd
   print "listpde"             , listpde
   print "listobjective"       , listobjective
   print "listw_0_field"       , listw_0_field 
   print "listk_0_field"       , listk_0_field 
   print "mc_filecmd"          , mc_filecmd 
   #use list comprehension to build entire set of parameter list
   paramlist =[ (meshcmditer,
                 k_0,k_1,k_2,k_3,w_0,w_n,w_i,w_d,w_2,w_ni,w_id,x_0,y_0,z_0,
                 mu_s,mu_a,anfact,g_flux,coeff_cool, method, optimize_w_0,
                 optimize_k_0, optimize_k_1, optimize_k_2, optimize_k_3, 
                 optimize_pow, optimize_mu_a, optimize_mu_s, 
                      w_0_field, k_0_field,k_1_ub,objective,pde)
                    for meshcmditer      in listmeshcmd   
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
                    for g_flux           in listg_flux 
                    for coeff_cool       in listcoeff_cool
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
                    for k_1_ub           in listk_1_ub        
                    for objective        in listobjective        
                    for pde              in listpde       ]

   if( len(paramlist)  > 40 ) : 
      print "\n\n    %d jobs > 40 job max at TACC \n\n" % len(paramlist) 
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
   #create batch job list
   id = 0 
   #create script to process the optimization iteration function value outputs
   fcnvalfile=open("%s/printvalue.txt" % jobid ,"w")
   #create gnuplot script to visualize iteration history
   gnuplotfile=open("%s/%s.plt" % (jobid,jobid) ,"w")
   gnuplotfile.write('# Gnuplot script file for "%s/%s??.dat" \n'
                                                               % (jobid,jobid) )
   gnuplotfile.write('# use set term <n> to switch between windows \n' )
   gnuplotfile.write('set   autoscale    # scale axes automatically   \n' )
   gnuplotfile.write('unset log          # remove any log-scaling     \n' )
   gnuplotfile.write('unset label        # remove any previous labels \n' )
   gnuplotfile.write('set xtic auto      # set xtics automatically    \n' )
   gnuplotfile.write('set ytic auto      # set ytics automatically    \n' )
   gnuplotfile.write('set xlabel "iteration #" \n' )
   gnuplotfile.write('#set xr [0:30] \n#set yr [0:100]\n' )
   gnuplotfile.write('set key below\n' )
   gnuplotfunc=[] ; gnuplotiter=[] ; gnuplotmu_s=[] ; gnuplotmu_a=[] ;
   gnuplotw_0=[]  ; gnuplotk_0=[]  ; gnuplotw_n=[]  ; gnuplotw_i=[]
   gnuplotw_d=[]  ; gnuplotw_2=[]  ; gnuplotk_1=[]  ; gnuplotk_2=[] 
   gnuplotk_3=[]  ; gnuplotx_0=[]  ; gnuploty_0=[]  ; gnuplotz_0=[]
   gnuplotw_ni=[] ; gnuplotw_id=[] ;
   for (meshcmditer,
                 k_0,k_1,k_2,k_3,w_0,w_n,w_i,w_d,w_2,w_ni,w_id,x_0,y_0,z_0,
                 mu_s,mu_a,anfact,g_flux,coeff_cool, method, optimize_w_0,
                 optimize_k_0, optimize_k_1, optimize_k_2, optimize_k_3, 
                 optimize_pow, optimize_mu_a, optimize_mu_s, 
                   w_0_field, k_0_field,k_1_ub,objective,pde) in paramlist:
      # create directory hierarchy to store files
      namejob= "%s%02d" % (jobid,id)
      fcnvalfile.write("echo %s using %s \n" % (namejob,pde))
      fcnvalfile.write("grep 'Function value' %s/out.o* | head -1 \n" % namejob)
      fcnvalfile.write("grep  func_eval       %s/out.o* | tail -1 \n" % namejob)
      fcnvalfile.write("grep 'Function value' %s/out.o* | tail -1 \n" % namejob)
      fcnvalfile.write("grep 'Function value' %s/out.o* | cut  -d ' ' -f3,6 | cut -d ',' -f'1 2' --output-delimiter=' ' > %s/iter.dat \n" % (namejob,namejob))
      fcnvalfile.write("grep 'Function evaluation' %s/out.o* | cut  -d ':' -f3 | awk '{print NR $0}' > %s/func.dat \n" % (namejob,namejob))
      constfldparamextract= "grep '%s(   0)' %s/files/qoi_0out_0.o* | cut -d '=' -f2 | awk '{print NR $0}' > %s/%s.dat \n"
      fcnvalfile.write(constfldparamextract% ("W_0",namejob,namejob,"W_0"))
      fcnvalfile.write(constfldparamextract% ("K_0",namejob,namejob,"K_0"))
      paramextract= "grep %s %s/files/qoi_0out_0.o*  | grep -v read_control | grep -v NOTE | awk '{print NR $0}' > %s/%s.dat \n"
      fcnvalfile.write(paramextract % ("K_1=" ,namejob,namejob,"K_1"  ))
      fcnvalfile.write(paramextract % ("K_2=" ,namejob,namejob,"K_2"  ))
      fcnvalfile.write(paramextract % ("K_3=" ,namejob,namejob,"K_3"  ))
      fcnvalfile.write(paramextract % ("W_N=" ,namejob,namejob,"W_N"  ))
      fcnvalfile.write(paramextract % ("W_I=" ,namejob,namejob,"W_I"  ))
      fcnvalfile.write(paramextract % ("W_D=" ,namejob,namejob,"W_D"  ))
      fcnvalfile.write(paramextract % ("W_2=" ,namejob,namejob,"W_2"  ))
      fcnvalfile.write(paramextract % ("W_NI=",namejob,namejob,"W_NID"))
      fcnvalfile.write(paramextract % ("X_0=" ,namejob,namejob,"X_0"  ))
      fcnvalfile.write(paramextract % ("Y_0=" ,namejob,namejob,"Y_0"  ))
      fcnvalfile.write(paramextract % ("Z_0=" ,namejob,namejob,"Z_0"  ))
      fcnvalfile.write(paramextract % ("MU_A=",namejob,namejob,"MU_A" ))
      fcnvalfile.write(paramextract % ("MU_S=",namejob,namejob,"MU_S" ))
      fcnvalfile.write("grep 'POW(' %s/files/qoi_0out_0.o*  | cut -d')' -f1,2 --output-delimiter=' ' | cut -d'=' -f1,2 --output-delimiter=' '  | cut -d'(' -f2,3  > %s/powhist.his \n" % (namejob,namejob) )
      print "creating job %s using %s " % (namejob,pde)
      id = id + 1 # update counter
      # gnuplot data
      exec "legend = '"+",".join(paramstudyid)+"' % (" + ",".join(paramstudvar) + ")"
      gnuplotfunc.append('"%s/func.dat" using 1:2 title "%s %s" w lp lw 2 ' % \
                                                  (namejob,namejob,legend) )
      gnuplotiter.append('"%s/iter.dat" using 1:2 title "%s %s" w lp lw 2 ' % \
                                                  (namejob,namejob,legend) )
      paramplot='"%s/%s.dat" using 1:%d title "%s %s" w lp lw 2 '
      gnuplotw_0.append( paramplot % (namejob,"W_0"  ,2,namejob,legend))
      gnuplotk_0.append( paramplot % (namejob,"K_0"  ,2,namejob,legend))
      gnuplotw_n.append( paramplot % (namejob,"W_N"  ,3,namejob,legend))
      gnuplotw_i.append( paramplot % (namejob,"W_I"  ,3,namejob,legend))
      gnuplotw_d.append( paramplot % (namejob,"W_D"  ,3,namejob,legend))
      gnuplotw_2.append( paramplot % (namejob,"W_2"  ,3,namejob,legend))
      gnuplotk_1.append( paramplot % (namejob,"K_1"  ,3,namejob,legend))
      gnuplotk_2.append( paramplot % (namejob,"K_2"  ,3,namejob,legend))
      gnuplotk_3.append( paramplot % (namejob,"K_3"  ,3,namejob,legend))
      gnuplotx_0.append( paramplot % (namejob,"X_0"  ,3,namejob,legend))
      gnuploty_0.append( paramplot % (namejob,"Y_0"  ,3,namejob,legend))
      gnuplotz_0.append( paramplot % (namejob,"Z_0"  ,3,namejob,legend))
      gnuplotw_ni.append(paramplot % (namejob,"W_NID",3,namejob,legend))
      gnuplotw_id.append(paramplot % (namejob,"W_NID",5,namejob,legend))
      gnuplotmu_s.append(paramplot % (namejob,"MU_S" ,3,namejob,legend))
      gnuplotmu_a.append(paramplot % (namejob,"MU_A" ,3,namejob,legend))
      utilities.create_directories(jobid,namejob)
      # functions are call by reference need to deepcopy
      cntrlfile = copy.deepcopy(config) 
      cntrlfile.set(     "field"      ,"k_0_field" ,  k_0_field        )
      cntrlfile.set(     "field"      ,"w_0_field" ,  w_0_field        )
      cntrlfile.set("thermalcond"     ,"k_0"       , "%f" % k_0        )
      cntrlfile.set("thermalcond"     ,"k_1"       , "%f" % k_1        )
      cntrlfile.set("thermalcond"     ,"k_1_ub"    , "%f" % k_1_ub     )
      cntrlfile.set("thermalcond"     ,"k_2"       , "%f" % k_2        )
      cntrlfile.set("thermalcond"     ,"k_3"       , "%f" % k_3        )
      cntrlfile.set("perfusivity"     ,"w_0"       , "%f" % w_0        )
      cntrlfile.set("perfusivity"     ,"w_n"       , "%f" % w_n        )
      cntrlfile.set("perfusivity"     ,"w_i"       , "%f" % w_i        )
      cntrlfile.set("perfusivity"     ,"w_d"       , "%f" % w_d        )
      cntrlfile.set("perfusivity"     ,"w_2"       , "%f" % w_2        )
      cntrlfile.set("perfusivity"     ,"w_ni"      , "%f" % w_ni       )
      cntrlfile.set("perfusivity"     ,"w_id"      , "%f" % w_id       )
      cntrlfile.set("source_laser"    ,"x_0"       , "%f" % x_0        ) 
      cntrlfile.set("source_laser"    ,"y_0"       , "%f" % y_0        )
      cntrlfile.set("source_laser"    ,"z_0"       , "%f" % z_0        )
      cntrlfile.set("source_laser"    ,"mu_a"      , "%f" % mu_a       )
      cntrlfile.set("source_laser"    ,"mu_s"      , "%f" % mu_s       )
      cntrlfile.set("source_laser"    ,"anfact"    , "%f" % anfact     )
      cntrlfile.set("neumann_boundary","g_flux"    , "%f" % g_flux     )
      cntrlfile.set("cauchy_boundary" ,"coeff_cool", "%f" % coeff_cool )
      run = cntrlfile.get( "compexec","run")
      cntrlfile.set("compexec" ,"run", run + method )
      meshcmd       = meshcmditer[0]
      numproclistid = meshcmditer[1]
      cntrlfile.set("compexec" ,"numproc", numproclist[numproclistid])
      # execute any mesh file and power file commands
      if(os.system(meshcmd % (namejob,namejob))):
         raise "\n\n    error with %s %% (%s,%s)" % (meshcmd,namejob,namejob)
      # execute any mc grid commands, the grid numbering from the database
      #  of mc grids must be used
      indexmu_a   = databasemu_a.index(mu_a)
      indexmu_s   = databasemu_s.index(mu_s)
      indexanfact = databaseanfact.index(anfact)
      ilocMCgrid = indexanfact + indexmu_a *len(databaseanfact) \
                             + indexmu_s*len(databaseanfact)*len(databasemu_a) 
      if(os.system(mc_filecmd % (ilocMCgrid,namejob))):
         raise "\n\n    error with %s %% (%d,%s)" % \
                                        (mc_filecmd,ilocMCgrid,namejob)
      cntrlfile.set("qoi_0","optimize_w_0" , optimize_w_0 )
      cntrlfile.set("qoi_0","optimize_k_0" , optimize_k_0 )
      cntrlfile.set("qoi_0","optimize_k_1" , optimize_k_1 )
      cntrlfile.set("qoi_0","optimize_k_2" , optimize_k_2 )
      cntrlfile.set("qoi_0","optimize_k_3" , optimize_k_3 )
      cntrlfile.set("qoi_0","optimize_pow" , optimize_pow )
      cntrlfile.set("qoi_0","optimize_mu_a", optimize_mu_a)
      cntrlfile.set("qoi_0","optimize_mu_s", optimize_mu_s)
      cntrlfile.set("qoi_0",  "objective"  ,   objective  )
      cntrlfile.set("hp3d" ,     "pde"     ,      pde     )
      joblist.append([namejob,cntrlfile])
   #close script
   fcnvalfile.close; fcnvalfile.flush()
   #append plot commands and close gnuplot file
   gnuplotfile.write('\nset title "optimization history" \n' )
   termid = 3
   gnuplotparams = [ ("W_0 [kg/s/m^3]",  gnuplotw_0 ),
                     ("k_0 [J/s/m/K]" ,  gnuplotk_0 ),
                     ("k_1 [J/s/m/K]" ,  gnuplotk_1 ),
                     ("k_2 [1/K]"     ,  gnuplotk_2 ),
                     ("k_3 [K]"       ,  gnuplotk_3 ),
                     ("w_n [kg/s/m^3]",  gnuplotw_n ),         
                     ("w_i [kg/s/m^3]",  gnuplotw_i ),
                     ("w_d [kg/s/m^3]",  gnuplotw_d ),
                     ("w_2 [kg/s/m^3]",  gnuplotw_2 ),
                     ("x_0 [m]"       ,  gnuplotx_0 ),
                     ("y_0 [m]"       ,  gnuploty_0 ),
                     ("z_0 [m]"       ,  gnuplotz_0 ),
                     ("mu_a [m^-1]"   ,  gnuplotmu_a),
                     ("mu_s [m^-1]"   ,  gnuplotmu_s),
                     ("w_ni [1/K]"    ,  gnuplotw_ni),
                     ("w_id [1/K]"    ,  gnuplotw_id) ]
   for (gnuparam,gnuparamcmds) in gnuplotparams:
      gnuplotfile.write('set ylabel "%s" \n' % gnuparam)
      gnuplotfile.write('set term x11 %d title "%s (%d)"\nplot \\\n' % (termid,gnuparam,termid) )
      gnuplotfile.write(", \\\n".join(gnuparamcmds))
      gnuplotfile.write("\nshow plot add2history\n")
      termid = termid + 1 
   gnuplotfile.write('\nset ylabel "function value [qoi units]" ' )
   gnuplotfile.write('\nset title "optimization function history" \n' )
   gnuplotfile.write('set term x11 2 title "function history" \nplot \\\n' )
   gnuplotfile.write(", \\\n".join(gnuplotfunc))
   gnuplotfile.write("\nshow plot add2history\n")
   gnuplotfile.write('\nset title "optimization iteration history" \n' )
   gnuplotfile.write('set term x11 1 title "iteration history"\nplot \\\n' )
   gnuplotfile.write(", \\\n".join(gnuplotiter))
   gnuplotfile.write("\nshow plot add2history\n")
   gnuplotfile.close; gnuplotfile.flush()
   #write python script for power plotting 
   pwrplt=open("%s/power.py" % jobid,"w")
   pwrplt.write(utilities.plotpower())
   pwrplt.close;pwrplt.flush() # close & ensure entire file written 
   #write gnuplot script for viewing optimizer history
   hisplt=open("%s/history.plt" % jobid,"w")
   hisplt.write(utilities.singlehistory())
   hisplt.close;hisplt.flush() # close & ensure entire file written 
   #make executable
   fcnvalfile=os.chmod("%s/printvalue.txt" % jobid ,00770)

   return joblist

def distribute_nopts(config):
   print """distribute_nopts:
      setting up to distribute each optimization step into several batch jobs 
   """
   joblist=[]
   workdir = config.get( "compexec" , "workdir" ) 
   jobid   = config.get( "compexec" , "jobid" ) 
   #get data used for batch submission
   noptsteps   = config.getint( "qoi_0" , "noptsteps"        ) 
   noffset     = config.getint( "qoi_0" , "noffset"          ) 
   numideal    = config.getint( "qoi_0" , "numideal"         ) 
   ideal_nzero = config.getint( "qoi_0" , "ideal_nzero_init" ) 
   ideal_ntime = config.getint( "qoi_0" , "ideal_ntime_init" ) 
   #create batch job list
   for i in range(noptsteps): # range stops at n - 1
      # create directory hierarchy to store files
      namejob= "%s%d" % (jobid,i)
      utilities.create_directories(jobid,namejob)
      # functions are call by reference need to deepcopy
      cntrlfile = copy.deepcopy(config) 
      cntrlfile.set("compexec","num_qoi"         ,"1"              )
      cntrlfile.set("qoi_0","noptsteps"       ,"1"              )
      cntrlfile.set("qoi_0","ideal_nzero_init","%d" % (ideal_nzero+i*noffset ) )
      cntrlfile.set("qoi_0","ideal_ntime_init","%d" % (ideal_ntime+i*numideal) )
      # compfilelocation gives the computation host the path to the 
      # use the same mesh file and power file for all batch runs
      cntrlfile.set("compexec","compfilelocation","%s/%s" % (workdir,jobid) )
      joblist.append([namejob,cntrlfile])
   #meshdata and powerdata used to point to files
   #copy the mesh file and the power file to the working directory
   meshfile =config.get("compexec","meshdata" ) 
   powerfile=config.get("compexec","powerdata") 
   os.system('cp %s %s/input_compact' % (meshfile ,jobid))
   os.system('cp %s %s/power.dat' % (powerfile ,jobid))
   return [joblist,["",""]]

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
