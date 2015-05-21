"""LBIE_2_hp3d.py script

   Transform LBIE mesh to hp3d format and apply boundary conditions

      usage: 
        
      python LBIE_2_hp3d.py  polyorder refinementlevel 
                             [filename].rawn [registrationfile].ini  
            where [filename].rawn contains the mesh  
                  polyorder is the polynomial order of mesh  
                  refinementlevel is the mesh refinement level  
  
"""
# Numerica imports intrinsic math functions(zeros,cos,sin,etc...)
import ConfigParser
import sys
import time
try:
  from numpy import *
except:
  from Numeric import *
  def cross(x,y):
      return  array([
                      x[1]*y[2] - x[2]*y[1],
                      x[2]*y[0] - x[0]*y[2],
                      x[0]*y[1] - x[1]*y[0],
                    ])


def elementcheck(oneelm,idelm):
    #compute vector giving orientation of 1-2-3-4 face wrt 5-6-7-8 face
    elem_dir = ( oneelm[5 - 1] + oneelm[6 - 1]  +       \
                 oneelm[7 - 1] + oneelm[8 - 1] )/4.0       \
              -( oneelm[1 - 1] + oneelm[2 - 1]  +       \
                 oneelm[3 - 1] + oneelm[4 - 1] )/4.0
    #compute normal of 1-2-3-4 face and make sure that the dot product
    # of the normal and element orientation is positive
    adum = oneelm[1] - oneelm[0]
    bdum = oneelm[3] - oneelm[0]
    normal = cross(adum,bdum)
    if(dot(normal,elem_dir) < 0.0):
      print "inverted element detected: element %d failed first check %le" % \
                                                  (idelm,dot(normal,elem_dir) )
      print "   element orientation ", elem_dir
      print "   element normal      ", normal
      return True
    # double check w/ possibly different plane
    adum = oneelm[3] - oneelm[2]
    bdum = oneelm[1] - oneelm[2]
    normal = cross(adum,bdum)
    if(dot(normal,elem_dir) < 0.0):
      print "inverted element detected: element %d failed first check %le" % \
                                                  (idelm,dot(normal,elem_dir) )
      print "   element orientation ", elem_dir
      print "   element normal      ", normal
      return True
    return False
    
def linmap(x,xmin,xmax,ymin,ymax):
    delx=xmax-xmin
    y = (xmax - x)/delx*ymin + (x-xmin)/delx*ymax
    return y

xhat=zeros(3,'d')
def register(A,x,b): # A x + b
    xhat[0] =  A[0][0] * x[0] + A[0][1] * x[1] + A[0][2] * x[2]  +  b[0]
    xhat[1] =  A[1][0] * x[0] + A[1][1] * x[1] + A[1][2] * x[2]  +  b[1]
    xhat[2] =  A[2][0] * x[0] + A[2][1] * x[1] + A[2][2] * x[2]  +  b[2]
    return xhat
if (len(sys.argv) < 4 ):
    assert False, __doc__
# polynomial order and refinement data for hp3d and number of processors 

p         = int( sys.argv[1] ) 
ref       = int( sys.argv[2] ) 
meshfile  =      sys.argv[3] 

# echo refinement level and polynomial order
print "meshfile %s " % meshfile
print "polynomial order   = %d " % p
print "element refinement = %d " % ref
print " "
filename = meshfile.split("/").pop()
filebase = filename.split(".")[0]
fileext  = filename.split(".")[1]
#import data
LBIE=open(meshfile,"r")
data = [line.split() for line in LBIE]
LBIE.close()
numnod = int(data[0][0]) #number of nodes
numelm = int(data[0][1]) #number of elements
print "number of nodes =",numnod ,"number of elements =", numelm
#initialize storage for nodes and element and bc data
nodes=zeros((numnod,3),'d')
nodbc=zeros((numnod,1),'i')
elems=zeros((numelm,8),'i')
if(fileext == 'raw'):
  print "RAW FILE DETECTED"
  column_bc_ID = 3
elif(fileext == 'rawn'):
  print "NORMAL DATA DETECTED"
  column_bc_ID = 6
else:
  print "UNKNOWN FILE TYPE "
  raise
#store nodes
for i in range(numnod):
  nodes[i]=[float(data[i+1][0]),float(data[i+1][1]),float(data[i+1][2])]
  nodbc[i]=int(data[i+1][column_bc_ID])
#get bounds on the data
mins=argmin(nodes,0)   ; maxs=argmax(nodes,0)
xmin=nodes[mins[0]][0] ; ymin=nodes[mins[1]][1] ; zmin=nodes[mins[2]][2]
xmax=nodes[maxs[0]][0] ; ymax=nodes[maxs[1]][1] ; zmax=nodes[maxs[2]][2]

#use test parser to get registration data
config = ConfigParser.ConfigParser()
#default the registration transformation to identity
config.add_section("registration_params")
config.set("registration_params",    "umin"  , "%f" % xmin )
config.set("registration_params",    "umax"  , "%f" % xmax )
config.set("registration_params",    "vmin"  , "%f" % ymin )
config.set("registration_params",    "vmax"  , "%f" % ymax )
config.set("registration_params",    "wmin"  , "%f" % zmin )
config.set("registration_params",    "wmax"  , "%f" % zmax )
config.set("registration_params", "RotMat00" ,     "1"     )
config.set("registration_params", "RotMat01" ,     "0"     )
config.set("registration_params", "RotMat02" ,     "0"     )
config.set("registration_params", "RotMat10" ,     "0"     )
config.set("registration_params", "RotMat11" ,     "1"     )
config.set("registration_params", "RotMat12" ,     "0"     )
config.set("registration_params", "RotMat20" ,     "0"     )
config.set("registration_params", "RotMat21" ,     "0"     )
config.set("registration_params", "RotMat22" ,     "1"     )
config.set("registration_params", "OffSet0"  ,     "0"     )
config.set("registration_params", "OffSet1"  ,     "0"     )
config.set("registration_params", "OffSet2"  ,     "0"     )
try: #try to open the control file containing registration data
  registrationfile = sys.argv[4]
  fname = open(registrationfile,"r")
  config.readfp(fname)
  registrationbase = "_%s" % registrationfile.split(".")[0]
except: #try to open the control file containing registration data
  print "NO REGISTRATION FILE FOUND"
  print "NO REGISTRATION FILE FOUND"
  print "NO REGISTRATION FILE FOUND"
  registrationbase = ""
#if (scalecase==1): 
#  print "scale case  1"
#  umin = -.063 ; vmin = .078 ; wmin = -.083 ;
#  umax =  .068 ; vmax = .020 ; wmax =  .095 ;

umin = config.getfloat( "registration_params" ,"umin")
umax = config.getfloat( "registration_params" ,"umax")
vmin = config.getfloat( "registration_params" ,"vmin")
vmax = config.getfloat( "registration_params" ,"vmax")
wmin = config.getfloat( "registration_params" ,"wmin")
wmax = config.getfloat( "registration_params" ,"wmax")
 
print "[umin,umax] =" , [umin,umax]
print "[vmin,vmax] =" , [vmin,vmax]
print "[wmin,wmax] =" , [wmin,wmax]

#scale the mesh
for i in range(numnod):
  nodes[i]=[linmap(float(data[i+1][0]),xmin,xmax,umin,umax),
            linmap(float(data[i+1][1]),ymin,ymax,vmin,vmax),
            linmap(float(data[i+1][2]),zmin,zmax,wmin,wmax)]

##get rotation matrix and offset
RotMat=zeros((3,3),'d')
OffSet=zeros(3,'d')
RotMat[0][0] = config.getfloat("registration_params", "RotMat00" )
RotMat[0][1] = config.getfloat("registration_params", "RotMat01" )
RotMat[0][2] = config.getfloat("registration_params", "RotMat02" )
RotMat[1][0] = config.getfloat("registration_params", "RotMat10" )
RotMat[1][1] = config.getfloat("registration_params", "RotMat11" )
RotMat[1][2] = config.getfloat("registration_params", "RotMat12" )
RotMat[2][0] = config.getfloat("registration_params", "RotMat20" )
RotMat[2][1] = config.getfloat("registration_params", "RotMat21" )
RotMat[2][2] = config.getfloat("registration_params", "RotMat22" )
OffSet[0]    = config.getfloat("registration_params", "OffSet0"  )
OffSet[1]    = config.getfloat("registration_params", "OffSet1"  )
OffSet[2]    = config.getfloat("registration_params", "OffSet2"  )

print "OffSet"   , OffSet
print "RotMat\n" , RotMat

for i in range(numnod):
  nodes[i]=register(RotMat,nodes[i],OffSet)

#store elements
for i in range(numelm):
  j = i+numnod  #jessica file format
  elems[i]=[ int(data[j+1][0]),int(data[j+1][1]),
             int(data[j+1][2]),int(data[j+1][3]),
             int(data[j+1][4]),int(data[j+1][5]),
             int(data[j+1][6]),int(data[j+1][7]) ]

print time.time()," Checking Element Quality "
errordetected = False
#compute element centroid
oneelm=zeros((8,3),'d')
elemscentroid=zeros((numelm,3),'d')
for i in range(numelm):
    for j in range(8):
        oneelm[j]= nodes[elems[i][j]]
    if(elementcheck(oneelm,i)):
          errordetected = True
    elemscentroid[i]=average(oneelm,0)

# domain decomposition embedded in FEM code

mapelm=range(numelm) # debugging


#                        z
#
#                        |
#                        |
#                        |
#                        |
#                        |               
#                      4 ___________________________ 7
#                       /?                         /|
#                      / ?                        / |
#                     /  ?                       /  |
#                    /   ?                      /   |
#                   /    ?                     /    |
#                  /     ?                    /     |
#               5 /______?__________________6/      |
#                |       ?                  |       |
#                |       ?                  |       |
#                |       ?                  |       |
#                |       ?                  |       |
#                |       ?                  |       |
#                |     0 ???????????????????|???????| 3  --------- y
#                |       ?                  |      /
#                |      ?                   |     /
#                |     ?                    |    /
#                |    ?                     |   /    
#                |   ?                      |  /
#                |  ?                       | /
#                | ?                        |/
#                |?_________________________/
#              1 /                           2
#               /
#              /
#             /       see ticam report 99-29 for face/node ordering
#            /
#          x
#temp storage for computing cross product on faces
adum=zeros(3,'d')
bdum=zeros(3,'d')
normal=zeros(3,'d')
#master element faces
face=zeros((6,4),'i')
face[0]=[0,1,2,3]  # face 1
face[1]=[4,5,6,7]  # face 2
face[2]=[0,4,5,1]  # face 3
face[3]=[1,2,6,5]  # face 4
face[4]=[3,7,6,2]  # face 5
face[5]=[0,3,7,4]  # face 6
print time.time()," Applying BC "
#if all four nodes of a face are flagged as bc then the face is a boundary
bc=zeros((numelm),'i') 
#store output
for i in range(numelm):
   k = mapelm[i]
   bc_nat=zeros((6),'i') 
   for j in range(6):
      bcflag=(nodbc[elems[k][face[j][0]]]+nodbc[elems[k][face[j][1]]]+
                nodbc[elems[k][face[j][2]]]+nodbc[elems[k][face[j][3]]] )
      adum= nodes[elems[k][face[j][1]]] - nodes[elems[k][face[j][0]]]
      bdum= nodes[elems[k][face[j][3]]] - nodes[elems[k][face[j][0]]]
      adum= 1./sqrt(dot(adum,adum))*adum
      bdum= 1./sqrt(dot(bdum,bdum))*bdum
      if (bcflag==0): 
         bc_nat[j]=0; # internal face 
      elif (bcflag==1):
         bc_nat[j]=0; # must be an internal face 
      elif (bcflag==2):
         bc_nat[j]=0; # internal face 
      elif (bcflag==3):
         bc_nat[j]=0; # internal face 
      elif (bcflag==4):
         normal[0] = adum[1]*bdum[2] - adum[2]*bdum[1] 
         normal[1] = adum[2]*bdum[0] - adum[0]*bdum[2] 
         normal[2] = adum[0]*bdum[1] - adum[1]*bdum[0] 
         #normalize
         normal= 1./sqrt(dot(normal,normal))*normal
         if (abs(dot(normal,[0.,0.,1.])) <= .98):
           bc_nat[j]=3;
         else: 
           bc_nat[j]=2;
      else:
         print '%d nodes on boundry: elem %d face %d' % (bcflag,i+1,j)
   bc[i]=  bc_nat[0]+bc_nat[1]*10+bc_nat[2]*100+bc_nat[3]*1000+bc_nat[4]*10000+bc_nat[5]*100000

print time.time()," writing files " 
# nodes start from 1 in gmv and hp3d
elems=elems+1
####################writout gmv file
gmv=open("%s%s.%s"% (filebase,registrationbase,"txt"),"w")
gmv.write("%s\n" % 'gmvinput ascii') 
gmv.write("nodev %d\n"%numnod) 
for i in range(numnod):
    gmv.write("%f %f %f \n" % (nodes[i][0] , nodes[i][1] , nodes[i][2] ))
gmv.write("cells %d\n"%numelm) 
for i in range(numelm):
    j = mapelm[i]
    gmv.write("hex 8\n  %d  %d  %d  %d  %d  %d  %d  %d \n" % 
        ( elems[j][0], elems[j][1], elems[j][2], elems[j][3],
          elems[j][4], elems[j][5], elems[j][6], elems[j][7]) )
gmv.write("%s\n" % 'variable') 
gmv.write("node__bc %d\n" % 1) 
for i in range(numnod):
    gmv.write("%d   " % nodbc[i]) 
gmv.write("\n") 
gmv.write("elemntbc %d\n" % 0) 
for i in range(numelm):
    gmv.write("%d   " % bc[i]) 
gmv.write("\n %s\n" % 'endvars') 
gmv.write("%s\n"%'endgmv') 
gmv.close
####################writout avs file
avs=open("%s%s.%s"% (filebase,registrationbase,"inp"),"w")
avs.write("%d %d 1 1 0\n" % (numnod,numelm)) 
#nodes
for i in range(numnod):
    avs.write("  %d %f %f %f \n" % (i+1,nodes[i][0],nodes[i][1],nodes[i][2] ))
#elements
for i in range(numelm):
    j = mapelm[i]
    avs.write("%d 1 hex %d  %d  %d  %d  %d  %d  %d  %d \n" % 
        ( i+1, elems[j][0], elems[j][1], elems[j][2], elems[j][3],
               elems[j][4], elems[j][5], elems[j][6], elems[j][7]) )
avs.write("   1   1\n")  # node variables
avs.write("node__bc, none \n") 
for i in range(numnod):
    avs.write("    %d  %d\n" % (i+1,nodbc[i])) 
avs.write("   1   1\n")  # element variables
avs.write("elemntbc, none \n") 
for i in range(numelm):
    avs.write("    %d  %d\n" % (i+1,bc[i])) 
avs.close
####################writout hp3d file 
## open the input deck file
hp3d=open("input_compact%s" % registrationbase ,"w")
## write out the dimensions and manifold dimensions
hp3d.write(" %d  %d  %s \n\n" % (3,3,'...NDIM,MANDIM') )
## write out the surfaces
hp3d.write(" %d   %s \n\n" % (0,'...NUMBER OF SURFACES') )
#
## write out the GMP hex nodes
hp3d.write(' %d   %s \n\n' % (numnod,'...NUMBER OF POINTS')  )
for i in range(numnod):
    hp3d.write('%s   %s   %d \n' % ('Regular','...TYPE OF POINT',i) )
    hp3d.write("%f %f %f  ...COORDINATES OF POINTS\n\n" %
                                        (nodes[i][0],nodes[i][1],nodes[i][2])  )
##write out connectivites
hp3d.write('\n\n  %d   %s \n\n' % (numelm,'...NUMBER OF HEXAHEDRONS') )
for i in range(len(elems)):
    j = mapelm[i]
    hp3d.write(' %d  %d  %d  %d  %d  %d  %d  %d EIGHT VERTICES OF HEX\n' % 
                        ( elems[j][0], elems[j][1], elems[j][2], elems[j][3],
                          elems[j][4], elems[j][5], elems[j][6], elems[j][7]) )
    hp3d.write(' %d  %d  %d  %d  %d  %d  \n\n' %  (0,0,0,0,0,0) )
##write out polynomial order and boundary conditions
hp3d.write('\n\n  %d   %s \n\n'%(numelm,'...NUMBER OF BLOCKS IN THE MESH') )
for i in range(numelm):
     hp3d.write('%d   %s   %d \n' % (10*(i+1)+2,'...BLOCK NICKNAME',i) )
     hp3d.write('   %d  %d  %d   %d %d %d  ...APPROX and SUBDIVISION NO\n' % 
		                                           (p,p,p,ref,ref,ref) )
     #@@@@  recall BC ARE ORDERED BACKWARDS!!!!!!!!!!!!!!!!!!!!!!!
     hp3d.write('   %d  ...BOUNDARY CONDITIONS\n\n' % bc[i] )
#space for direct solver (not used)
hp3d.write('\n %d   %s\n' % (100,'...SPACE FOR FRONTAL SOLVER(not used in petsc)'))
#estimate number of elements, vertices, and higher order nodes needed
NRELIB= numelm*pow(ref,3)
NRVERB= 8*numelm*pow(ref,3)
NRNODB= 19*numelm*pow(ref,3)
hp3d.write('\n %d %d %d %d  %s \n' % \
                (NRELIB,NRVERB,NRNODB,1,'...MAXELIB,MAXVERB,MAXNODB,NREQNB') )
hp3d.write('\n %d   %s \n' % (0,'     NELNEW') )

# close the input deck file
hp3d.close
