import itkUtilities # custom python module for access to itk library
import sys

#dataFile = "/mnt/data/scratch/fuentes/data/dir/lung_may10/Cases/u17441t/u17441t_sample150_Source_xyz.txt"
#deformationField = "/mnt/data/scratch/fuentes/data/dir/lung_may10/Cases/u17441t/DeformableRegistration15/deformationField.mha"
#outputFile = "output.txt"
# get file names from command line
dataFile         = sys.argv[1]
deformationField = sys.argv[2]
outputFile       = sys.argv[3]

#import data
dirData=open(dataFile,"r")
landMarkData = [map(int,line.split()) for line in dirData]
dirData.close()

#check extension
if (deformationField.split(".").pop() != "mha"):
   raise IOError("\n\n    only works w/ mha files?? " )
# read from mha file and transform indices
transformedfile=open( outputFile ,"w")
for landMark in landMarkData:
    offset= itkUtilities.GetPixelValue( deformationField  , landMark) 
    transformedfile.write("%d %d %d\n" % ( landMark[0]+offset[0],
                                           landMark[1]+offset[1],
                                           landMark[2]+offset[2])   )

transformedfile.close; transformedfile.flush() 
