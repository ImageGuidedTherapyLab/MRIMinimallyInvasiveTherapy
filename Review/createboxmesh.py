"""createboxmesh.py script

create a mesh of a box from 2 points
                                                        y 
          ____________                                 |
         |            |?                               |
         |            |  ?                             |
         |            |    ?                           |
         |            |      ?                         |__________  x
         |            |        ?                        \\
         |            |          ?                         \\
         |            |            ?                          \\
         |            |              ?                           \\  z
         |____________|                ?                 
       p1?                               ?             
           ?                                ?          
             ?                                 ?       
               ?                    ____________  p2   
                 ?                 |            |         
                   ?               |            |         
                     ?             |            |         
                       ?           |            |        
                         ?         |            |        
                           ?       |            |     
                             ?     |            |     
                               ?   |            |     
                                 ? |____________|     
                                                    
  p = uniform polynomial order       refine = refinement level of the GMP hexes
       p1  =  first point   ;   p2  = second point
       nx,ny,nz  =  # of GMP hexes in x,y,z direction

  [usage]
         python createboxmesh.py <inifile>

           if used inifile is expected to contain a [backupmesh] section
   
             with either   
                 boundingbox= bx=[?,?];by=[?,?];bz=[?,?]  (NO SPACES)
                        -or-
                 points= p1=[?,?,?];p2=[?,?,?]            (NO SPACES)

"""
print __doc__
import sys
import ConfigParser
try:
  from numpy import *
except:
  from Numeric import *
controlfile = ConfigParser.ConfigParser()
#set defaults
controlfile.add_section("backupmesh")
controlfile.set("backupmesh" ,  "p" ,  "1" )
controlfile.set("backupmesh" , "ref",  "1" )
controlfile.set("backupmesh" , "nx" , "65" )
controlfile.set("backupmesh" , "ny" , "65" )
controlfile.set("backupmesh" , "nz" , "15" )
try: # get data from an ini file 
   controlfile.readfp( open("%s" % sys.argv[1], "r") )
   X0 =  controlfile.getfloat("planning_images", "X0" )
   Y0 =  controlfile.getfloat("planning_images", "Y0" )
   Z0 =  controlfile.getfloat("planning_images", "Z0" )
   DX =  controlfile.getfloat("planning_images", "DX" )
   DY =  controlfile.getfloat("planning_images", "DY" )
   DZ =  controlfile.getfloat("planning_images", "DZ" )
   print controlfile.get("planning_images","boundingbox")
   exec  controlfile.get("planning_images","boundingbox")
   p1 = [X0 + bx[0]* DX, Y0 + by[0]* DY, Z0 + bz[0]* DZ ]
   p2 = [X0 + bx[1]* DX, Y0 + by[1]* DY, Z0 + bz[1]* DZ ]
   nz = bz[1]-bz[0] + 1 # ensure sufficiently resolved for thermal images
except ConfigParser.NoOptionError: 
   print "bounding box not found trying to read points directly..."
   print "expected of the form: points = p1=[p1x,p1y,p1z];p2=[p2x,p2y,p2z]"
   print controlfile.get("backupmesh","points")
   exec  controlfile.get("backupmesh","points")
   nz  =  controlfile.getint("backupmesh","nz" )
except IndexError: 
   print "control file not found..."
   p1 = [0.0,0.0,0.0]; p2 = [1.0,1.0,1.0]
   nz  =  controlfile.getint("backupmesh","nz" )
except: 
   print "error inputing bounding box..."
   p1 = [0.0,0.0,0.0]; p2 = [1.0,1.0,1.0]
   nz  =  controlfile.getint("backupmesh","nz" )
p   =  controlfile.getint("backupmesh", "p" )
ref =  controlfile.getint("backupmesh","ref")
nx  =  controlfile.getint("backupmesh","nx" )
ny  =  controlfile.getint("backupmesh","ny" )

# echo the data
print "polynomial order %d, hp3d refinement level %d"%(p,ref)
print "nx %d, ny %d, nz %d" % (nx,ny,nz)
print "point 1 ", p1
print "point 2 ", p2

# open the input deck file
hp3d=open("input_compact_box","w")
# write out the dimensions and manifold dimensions
hp3d.write(" %d  %d  %s \n\n" % (3,3,'...NDIM,MANDIM') )
# write out the surfaces
hp3d.write(" %d   %s \n\n" % (0,'...NUMBER OF SURFACES') )

# define nodes
nodes =  [];
nodbc=zeros( ( (nx+1)*(ny+1)*(nz+1) ,1 ) , 'i' )
icount = 0 

for k in range(nz+1):
   zcoord = float(nz-k)/float(nz) * p1[2] + float(k)/float(nz) * p2[2]
   for j in range(ny+1):
      ycoord = float(ny-j)/float(ny) * p1[1] + float(j)/float(ny) * p2[1]
      for i in range(nx+1):
          xcoord = float(nx-i)/float(nx) * p1[0] + float(i)/float(nx) * p2[0]
          nodes.append((xcoord,ycoord,zcoord))
          if(i== 0):
             nodbc[icount] = 3 
          if(i==nx):
             nodbc[icount] = 3 
          if(j== 0):
             nodbc[icount] = 3 
          if(j==ny):
             nodbc[icount] = 3 
          if(k== 0):
             nodbc[icount] = 2 
          if(k==nz):
             nodbc[icount] = 2 
          icount = icount + 1


# write out the GMP hex nodes
hp3d.write(' %d   %s \n\n' % (len(nodes),'...NUMBER OF POINTS')  )
for i in range(len(nodes)):
    hp3d.write('%s   %s   %d \n' % ('Regular','...TYPE OF POINT',i+1) )
    hp3d.write("%f %f %f  ...COORDINATES OF POINTS\n\n" % nodes[i] )
#
#                        z
#
#                        |
#                        |
#                        |
#                        |
#                        |           a16
#                     a5 _____________x____________  a8
#                       /?                         /|
#                      / ?                        / |
#                     /  ?                       /  |
#                 a13x   ?      a22x         a15x   |
#                   /    ?                     /    |
#                  /     ?    a14             /     |
#              a6 /______?_____x___________a7/      |
#                |    a17x         a26x     |    a20x
#                |       ?                  |       |
#                |       ?                  |       |
#                |a23x   ?     a27x         | a25x  |
#                |       ?           a12    |       |
#                |    a1 ????????????x??????|???????| a4 --------- y
#             a18x       ? a24x          a19x      /
#                |      ?                   |     /
#                |     ?                    |    /
#                |  a9x           x         |   x a11
#                |   ?           a21        |  /
#                |  ?                       | /             for BC and surfaces
#                | ?                        |/             --------------------
#                |?___________x_____________/                 a21 = face 1
#               a2           a10           a3                 a22 = face 2
#               /                                             a23 = face 3
#              /                                              a24 = face 4
#             /                                               a25 = face 5
#            /                                                a26 = face 6
#          x
# counter clockwise on nodes gives outward normal on all faces
# define element connectivities
elems = []; bc = [];                                        #  start numbering
for k in range(nz):                                         #   |    at ONE
   for j in range(ny):                                      #   |
      for i in range(nx):                                   #  \|/
         elems.append(  (  k  *(nx+1)*(ny+1)+  j  *(nx+1)+i   + 1,
                           k  *(nx+1)*(ny+1)+  j  *(nx+1)+i+1 + 1,
                           k  *(nx+1)*(ny+1)+(j+1)*(nx+1)+i+1 + 1,
                           k  *(nx+1)*(ny+1)+(j+1)*(nx+1)+i   + 1,
                         (k+1)*(nx+1)*(ny+1)+  j  *(nx+1)+i   + 1,
                         (k+1)*(nx+1)*(ny+1)+  j  *(nx+1)+i+1 + 1,
                         (k+1)*(nx+1)*(ny+1)+(j+1)*(nx+1)+i+1 + 1, 
                         (k+1)*(nx+1)*(ny+1)+(j+1)*(nx+1)+i   + 1 )  )
         bcx0=0; bcxN=0; bcy0=0; bcyN=0; bcz0=0; bczN=0;
         if(i==0):
            bcx0=3;
         if(i==nx-1):
            bcxN=3;
         if(j==0):
            bcy0=3;
         if(j==ny-1):
            bcyN=3;
         if(k==0):
            bcz0=2;
         if(k==nz-1):
            bczN=2;
         bc.append((bcz0,bczN,bcy0,bcxN,bcyN,bcx0))
elmsurf=(0,0,0,0,0,0)
#write out connectivites
hp3d.write('\n\n  %d   %s \n\n' % (len(elems),'...NUMBER OF HEXAHEDRONS') )
for i in range(len(elems)):
    hp3d.write(' %d  %d  %d  %d  %d  %d  %d  %d EIGHT VERTICES OF HEX\n' % elems[i])
    hp3d.write(' %d  %d  %d  %d  %d  %d  \n\n' %  elmsurf ) 
#write out polynomial order and boundary conditions
hp3d.write('\n\n  %d   %s \n\n'%(len(elems),'...NUMBER OF BLOCKS IN THE MESH') )
for i in range(len(elems)):
     hp3d.write('%d   %s   %d \n' % (10*(i+1)+2,'...BLOCK NICKNAME',i+1) )
     hp3d.write('   %d  %d  %d   %d %d %d  ...APPROX and SUBDIVISION NO\n' % 
		                                           (p,p,p,ref,ref,ref) )
     #@@@@  recall BC ARE ORDERED BACKWARDS!!!!!!!!!!!!!!!!!!!!!!!
     hp3d.write('   %d%d%d%d%d%d  ...BOUNDARY CONDITIONS\n\n' % 
		       (bc[i][5],bc[i][4],bc[i][3],bc[i][2],bc[i][1],bc[i][0]) )
#space for direct solver (not used)
hp3d.write('\n %d   %s\n' % (100,'...SPACE FOR FRONTAL SOLVER(not used in petsc)'))
#number of elements, vertices, and higher order nodes
# damn hp3d bugs
MAXELIB = max(len(elems)+10,100)
MAXVERB = max(8*len(elems)*pow(ref,3),1000)
MAXNODB = max(19*len(elems)*pow(ref,3),1000)
hp3d.write('\n %d %d %d %d  %s \n' %  (MAXELIB, MAXVERB, MAXNODB,1,
                                        '...MAXELIB,MAXVERB,MAXNODB,NREQNB') )
hp3d.write('\n %d   %s \n' % (0,'     NELNEW') )

# close the input deck file
hp3d.close


####################writout avs file
avs=open("box.inp","w")
avs.write("%d %d 1 1 0\n" % (len(nodes),len(elems))) 
#nodes
for i in range(len(nodes)):
    avs.write("  %d %f %f %f \n" % (i+1,nodes[i][0],nodes[i][1],nodes[i][2] ))
#elements
for i in range(len(elems)):
    avs.write("%d 1 hex %d  %d  %d  %d  %d  %d  %d  %d \n" % 
        ( i+1, elems[i][0], elems[i][1], elems[i][2], elems[i][3],
               elems[i][4], elems[i][5], elems[i][6], elems[i][7]) )
avs.write("   1   1\n")  # node variables
avs.write("node__bc, none \n") 
for i in range(len(nodes)):
    avs.write("    %d  %d\n" % (i+1,nodbc[i])) 
avs.write("   1   1\n")  # element variables
avs.write("elemntbc, none \n") 
for i in range(len(elems)):
    avs.write("    %d  %d%d%d%d%d%d\n" % 
	         (i+1,bc[i][5],bc[i][4],bc[i][3],bc[i][2],bc[i][1],bc[i][0]) )
avs.close

