import os
import time

# location of local directory
localdirectory='datatest/'

# setup remote server
#   transfer FROM remotedirectory to localdirectory
username='sdc'
remoteserver='10.115.24.111'
username='fuentes'
remoteserver='10.115.8.170'
remotedirectory='/home/cjmaclellan/mfgreData/datatest'

# all connections will go through this file
socketfile='/tmp/%r@%h:%p'

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
