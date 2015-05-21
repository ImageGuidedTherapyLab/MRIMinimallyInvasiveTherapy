import os
import time


# setup command line parser to control execution
from optparse import OptionParser
parser = OptionParser()
parser.add_option( "--socketpath",
                  action="store", dest="socketpath", default=None,
                  help="all connections will go through this file", metavar="FILE")
parser.add_option( "--localserver", 
                  action="store", dest="localserver", default=None,
                  help="transfer FROM localserver:localdirectory TO  remoteuser@remoteserver:remotedirectory", metavar="IP")
parser.add_option( "--localdirectory", 
                  action="store", dest="localdirectory", default=None,
                  help="transfer FROM localserver:localdirectory TO  remoteuser@remoteserver:remotedirectory", metavar="PATH")
parser.add_option( "--remoteuser", 
                  action="store", dest="remoteuser", default=None,
                  help="transfer FROM localserver:localdirectory TO  remoteuser@remoteserver:remotedirectory", metavar="USER")
parser.add_option( "--remotedirectory", 
                  action="store", dest="remotedirectory", default=None,
                  help="transfer FROM localserver:localdirectory TO  remoteuser@remoteserver:remotedirectory", metavar="PATH")
parser.add_option( "--remoteserver", 
                  action="store", dest="remoteserver", default=None,
                  help="transfer FROM localserver:localdirectory TO  remoteuser@remoteserver:remotedirectory", metavar="IP")
(options, args) = parser.parse_args()

DefaultPort=22
if (options.localserver  != None and options.localdirectory  != None
                                 and options.socketpath      != None
    options.remoteserver != None and options.remotedirectory != None ):
  # all connections will go through this file
  socketfile = "%s/%s@%s:%s" % (options.socketpath,options.remoteuser,options.remoteserver,DefaultPort) 
  # create persistent connections if they do not exist...
  #    ssh -O exit to kill a connection
  CheckConnectionCMD = "ssh -O check -S %s %s@%s" % (socketfile,username,remoteserver) 
  KillConnectionCMD  = "ssh -O exit  -S %s %s@%s" % (socketfile,username,remoteserver) 
  print CheckConnectionCMD 
  if( os.system(CheckConnectionCMD) ):
     print "\n\n   Creating Persisent Connection on %s!!!! \n\n " % remoteserver
     CreateConnectionCMD='ssh -MNf -S %s %s@%s ' %  (socketfile,username,remoteserver) 
     os.system(CreateConnectionCMD)
  else:
     print "\n\n   Found Persisent Connection on %s!!!! \n\n " % remoteserver
     print "Kill a connection: %s " % KillConnectionCMD 
  
  RunRsync=True
  while(RunRsync):
    try:
      # ensure local directory available
      localmkdir= 'mkdir -p %s  ' % (localdirectory)
      print localmkdir
      os.system( localmkdir )
      # start rsync
      rsynccmd = 'rsync -e "ssh -S %s" -avz %s@%s:%s/ %s/ ' %  (socketfile,username,remoteserver,remotedirectory,localdirectory)
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

