# This project is designed to be built outside the Insight source tree.
PROJECT(Registration)

# Echo Compiler Options 
SET(CMAKE_VERBOSE_MAKEFILE ON)

# Find ITK.
FIND_PACKAGE(ITK REQUIRED)
IF(ITK_FOUND)
  INCLUDE(${ITK_USE_FILE})
ENDIF(ITK_FOUND)

# Find X11 libraries to link with
FIND_PACKAGE(X11 REQUIRED)
IF(NOT X11_FOUND)
   MESSAGE("Could not find X11 package...")
ENDIF(NOT X11_FOUND)

ADD_EXECUTABLE(RigidRegistration RigidRegistration.cxx )
TARGET_LINK_LIBRARIES(RigidRegistration ITKNumerics ITKIO)

ADD_EXECUTABLE(DicomSeriesConvert DicomSeriesConvert.cxx )
TARGET_LINK_LIBRARIES(DicomSeriesConvert ITKCommon ITKIO X11)

ADD_EXECUTABLE(DicomCopyHeader DicomCopyHeader.cxx)
TARGET_LINK_LIBRARIES(DicomCopyHeader ITKCommon ITKIO)

ADD_EXECUTABLE(DicomSeriesCopyHeader DicomSeriesCopyHeader.cxx)
TARGET_LINK_LIBRARIES(DicomSeriesCopyHeader ITKCommon ITKIO)

ADD_EXECUTABLE(SubsampleNoFilter SubsampleNoFilter.cxx )
TARGET_LINK_LIBRARIES(SubsampleNoFilter  ITKBasicFilters ITKIO)
