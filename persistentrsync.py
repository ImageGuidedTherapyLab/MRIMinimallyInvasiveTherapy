import os
import time


# setup command line parser to control execution
from optparse import OptionParser
parser = OptionParser()
parser.add_option( "--socketpath",
                  action="store", dest="socketpath", default=None,
                  help="all connections will go through this file", metavar="FILE")
parser.add_option( "--localdirectory", 
                  action="store", dest="localdirectory", default=None,
                  help="transfer FROM localhost:localdirectory TO  remoteuser@remoteserver:remotedirectory", metavar="PATH")
parser.add_option( "--remoteuser", 
                  action="store", dest="remoteuser", default=None,
                  help="transfer FROM localhost:localdirectory TO  remoteuser@remoteserver:remotedirectory", metavar="USER")
parser.add_option( "--remotedirectory", 
                  action="store", dest="remotedirectory", default=None,
                  help="transfer FROM localhost:localdirectory TO  remoteuser@remoteserver:remotedirectory", metavar="PATH")
parser.add_option( "--remoteserver", 
                  action="store", dest="remoteserver", default=None,
                  help="transfer FROM localhost:localdirectory TO  remoteuser@remoteserver:remotedirectory", metavar="IP")
parser.add_option( "--remotersync", 
                  action="store", dest="remotersync", default=None,
                  help="rsync install on remoteserver", metavar="PATH")
(options, args) = parser.parse_args()

DefaultPort=22
if (options.socketpath   != None and options.localdirectory  != None  
                                 and
    options.remoteserver != None and options.remotedirectory != None ):
  # all connections will go through this file
  socketfile = "%s/%s@%s:%s" % (options.socketpath,options.remoteuser,options.remoteserver,DefaultPort) 
  # create persistent connections if they do not exist...
  #    ssh -O exit to kill a connection
  CheckConnectionCMD = "ssh -O check -S %s %s@%s" % (socketfile,options.remoteuser,options.remoteserver) 
  KillConnectionCMD  = "ssh -O exit  -S %s %s@%s" % (socketfile,options.remoteuser,options.remoteserver) 
  print CheckConnectionCMD 
  if( os.system(CheckConnectionCMD) ):
     print "\n\n   Creating Persisent Connection on %s!!!! \n\n " % options.remoteserver
     CreateConnectionCMD='ssh -MNf -S %s %s@%s ' %  (socketfile,options.remoteuser,options.remoteserver) 
     os.system(CreateConnectionCMD)
  else:
     print "\n\n   Found Persisent Connection on %s!!!! \n\n " % options.remoteserver
  print "Kill a connection: %s " % KillConnectionCMD 
  
  RunRsync=True
  while(RunRsync):
    try:
      # ensure local directory available
      localmkdir= 'mkdir -p %s  ' % (options.localdirectory)
      print localmkdir
      os.system( localmkdir )
      # start rsync
      rsynccmd = 'rsync -e "ssh -S %s"  ' %  (socketfile)
      if(options.remotersync != None):
         rsynccmd = rsynccmd + ' --rsync-path=%s ' % (options.remotersync)
      rsynccmd = rsynccmd + ' -avz %s@%s:%s/ %s/ ' %  (options.remoteuser,options.remoteserver,options.remotedirectory,options.localdirectory)
      print rsynccmd 
      os.system( rsynccmd )
      time.sleep(1)
    except KeyboardInterrupt:
      # kill connection on KeyboardInterrupt
      RunRsync=False
      print KillConnectionCMD
      os.system(KillConnectionCMD)
else:
  parser.print_help()
  print options

