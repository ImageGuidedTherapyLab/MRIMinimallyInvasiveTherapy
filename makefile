include ${PETSC_DIR}/conf/base

.PHONY: tags

all: Makefile
	make -f Makefile
Makefile: CMakeLists.txt
	cmake -DITK_DIR=$(ITK_DIR) .
tags:
	ctags --langmap=c++:+.txx --languages=c,c++ -R $(PETSC_DIR) $(ITK_SOURCE) CImg.h *

ITK_CXXFLAGS=-ftemplate-depth-50 
# taken from cmake with verbose makefile
ITK_INCLUDE= -I$(ITK_HOME)/include/InsightToolkit/gdcm/src \
             -I$(ITK_HOME)/include/InsightToolkit/gdcm \
             -I$(ITK_HOME)/include/InsightToolkit/Utilities \
             -I$(ITK_HOME)/include/InsightToolkit/Utilities/vxl/core \
             -I$(ITK_HOME)/include/InsightToolkit/Utilities/vxl/vcl \
             -I$(ITK_HOME)/include/InsightToolkit/Utilities/vxl/v3p/netlib \
             -I$(ITK_HOME)/include/InsightToolkit/Utilities/itkExtHdrs \
             -I$(ITK_HOME)/include/InsightToolkit/Utilities/nifti/znzlib \
             -I$(ITK_HOME)/include/InsightToolkit/Utilities/nifti/niftilib \
             -I$(ITK_HOME)/include/InsightToolkit/Utilities/expat \
             -I$(ITK_HOME)/include/InsightToolkit/Utilities/DICOMParser \
             -I$(ITK_HOME)/include/InsightToolkit/Utilities/NrrdIO \
             -I$(ITK_HOME)/include/InsightToolkit/Utilities/MetaIO \
             -I$(ITK_HOME)/include/InsightToolkit/SpatialObject \
             -I$(ITK_HOME)/include/InsightToolkit/Numerics/NeuralNetworks \
             -I$(ITK_HOME)/include/InsightToolkit/Numerics/Statistics \
             -I$(ITK_HOME)/include/InsightToolkit/Numerics/FEM \
             -I$(ITK_HOME)/include/InsightToolkit/IO \
             -I$(ITK_HOME)/include/InsightToolkit/Numerics \
             -I$(ITK_HOME)/include/InsightToolkit/Common \
             -I$(ITK_HOME)/include/InsightToolkit/BasicFilters \
             -I$(ITK_HOME)/include/InsightToolkit/Algorithms \
             -I$(ITK_HOME)/include/InsightToolkit  
ITK_LIB= -rdynamic -L$(ITK_DIR) -lITKFEM -lITKIO -lITKStatistics -lITKNumerics -lITKBasicFilters -lITKNrrdIO -litkgdcm -litkjpeg12 -litkjpeg16 -litkopenjpeg -luuid -litkpng -litktiff -litkjpeg8 -lITKSpatialObject -lITKMetaIO -lITKDICOMParser -lITKEXPAT -lITKniftiio -lITKznz -litkzlib -lITKCommon -litkvnl_inst -litkvnl_algo -litkv3p_netlib -litkvnl -litkvcl -lm -litksys -lpthread -ldl -lm -Wl,-rpath,$(ITK_DIR)
#
INCLUDE= $(PETSC_INCLUDE) $(ITK_INCLUDE) 

%.o: %.cxx
	@echo Compiling $<
	$(PCC) $(PREC_FLAGS) $(PCC_FLAGS) $(ITK_CXXFLAGS) $(INCLUDE) -o $@ -c $<

DicomReadImagRealWriteTemp: DicomReadImagRealWriteTemp.o chkopts
	-${CLINKER} -o DicomReadImagRealWriteTemp DicomReadImagRealWriteTemp.o  ${PETSC_DM_LIB} ${ITK_LIB} 
