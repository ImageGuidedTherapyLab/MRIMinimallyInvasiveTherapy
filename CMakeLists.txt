cmake_minimum_required(VERSION 2.4)

# This project is designed to be built outside the Insight source tree.
PROJECT(RealTime)

# Echo Compiler Options 
SET(CMAKE_VERBOSE_MAKEFILE ON)

# Find ITK.
FIND_PACKAGE(ITK REQUIRED)
IF(ITK_FOUND)
  set(ITK_BUILD_SETTINGS_FILE)
  INCLUDE(${ITK_USE_FILE})
ENDIF(ITK_FOUND)

# Set Compile Method
STRING(COMPARE EQUAL "$ENV{METHOD}" "dbg" DBGMODE)
STRING(COMPARE EQUAL "$ENV{METHOD}" "opt" OPTMODE)
IF(DBGMODE)
   message("dbg mode")
   SET(CMAKE_BUILD_TYPE Debug)
ENDIF(DBGMODE)
IF(OPTMODE)
   message("opt mode")
   SET(CMAKE_BUILD_TYPE Release)
ENDIF(OPTMODE)

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

ADD_EXECUTABLE(DicomReadImagRealWriteTemp DicomReadImagRealWriteTemp.cxx )
TARGET_LINK_LIBRARIES(DicomReadImagRealWriteTemp ITKIO)
