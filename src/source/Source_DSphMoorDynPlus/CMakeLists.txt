#--------------------------------------------------------------
# CMAKE file to compile MoorDynPlus
#--------------------------------------------------------------

cmake_minimum_required(VERSION 3.0)

#--------------------------------------------------------------
# 1. Project name
#--------------------------------------------------------------
PROJECT(DSphMoorDyn)

#--------------------------------------------------------------
# 2. Configures library on execution time
#--------------------------------------------------------------
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MT")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MTd")

#--------------------------------------------------------------
# 3. Configures Build types for VS
#--------------------------------------------------------------
set(CMAKE_CONFIGURATION_TYPES "")
if(MSVC)
  set_property(CACHE CMAKE_CONFIGURATION_TYPES PROPERTY STRINGS "Release" "Debug" ) 
  # Cathes the build type selected by user
  set(CXX_FLAGS_MoorDyn "")
  foreach(CONFIG_TYPES CMAKE_CONFIGURATION_TYPES)
    if(${CONFIG_TYPES} MATCHES "Release")
      set(CXX_FLAGS_MoorDyn "${CMAKE_CXX_FLAGS_RELEASE}")
      set(CMAKE_CONFIGURATION_TYPES "Release" CACHE STRING "" FORCE)
    elseif(${CONFIG_TYPES} MATCHES "Debug")
      set(CXX_FLAGS_MoorDyn "${CMAKE_CXX_FLAGS_DEBUG}")
      set(CMAKE_CONFIGURATION_TYPES "Debug" CACHE STRING "" FORCE)
    endif()
  endforeach(CONFIG_TYPES)
endif()

#--------------------------------------------------------------
# 4. Source files
#--------------------------------------------------------------
# MoorDynPlus source files
set(SRC
  IEnvironment.cpp  
  IQSlines.cpp       
  IBody.cpp         
  ILine.cpp         
  kiss_fft.c         
  DSphMoorDynPlus.cpp
  Misc.cpp
  IConnection.cpp
  IOutput.cpp
  FunMoorDynPlus.cpp
  MoorDynPlus.cpp
)

# MoorDynPlus header files
set(HEADER 
  FunMoorDynPlus.h
  MoorDynPlus.h
  IApi.h
  IEnvironment.h
  IQSlines.h
  TypesMoorDynPlus.h
  IBody.h
  ILine.h
  _kiss_fft_guts.h 
  kiss_fft.h
  Misc.h
  IConnection.h
  IOutput.h
  kissfft.hh
)

# DualSPHysics header files
set(DSPH_HEADERS 
  Tools/Functions.h 
  Tools/FunGeo3d.h
  Tools/JDataArrays.h
  Tools/JException.h 
  Tools/JLog2.h 
  Tools/JObject.h 
  Tools/JRangeFilter.h
  Tools/JSaveCsv2.h 
  Tools/JXml.h 
  Tools/tinystr.h 
  Tools/tinyxml.h 
  Tools/TypesDef.h
)

# DualSPHysics source files
set(DSPH_SRC 
  Tools/Functions.cpp
  Tools/FunGeo3d.cpp
  Tools/JDataArrays.cpp
  Tools/JException.cpp
  Tools/JLog2.cpp
  Tools/JObject.cpp
  Tools/JRangeFilter.cpp
  Tools/JSaveCsv2.cpp
  Tools/JXml.cpp
  Tools/tinystr.cpp
  Tools/tinyxml.cpp
  Tools/tinyxmlerror.cpp
  Tools/tinyxmlparser.cpp
)

# Set the group files
source_group(core FILES ${SRC} ${HEADER})
source_group(Tools FILES ${DSPH_HEADERS} ${DSPH_SRC})

# Add the groups to project files
set(LibDSphMoorDyn_FILES ${SRC} ${HEADER} ${DSPH_HEADERS} ${DSPH_SRC})

# Include directories
set(CXX_INCLUDES "${CMAKE_CURRENT_SOURCE_DIR}")
set(CXX_INCLUDES  "${CXX_INCLUDES}/Tools " )
include_directories(${CXX_INCLUDES})

#--------------------------------------------------------------
# 5. CMake flags
#--------------------------------------------------------------
set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DDISABLE_DSPH")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DDISABLE_DSPH")

#--------------------------------------------------------------
# 6. Don't use C++11 ABI for compatibility between GCC v4.9 and >=5
#--------------------------------------------------------------
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  message(STATUS "GCC version:  ${CMAKE_CXX_COMPILER_VERSION}")
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -c -O3 -static -static-libgcc -DLINUX -fPIC ")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x -c -O3 -static -static-libgcc -std=gnu++0x -DLINUX -fPIC ")
  if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 4.9.9)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D_GLIBCXX_USE_CXX11_ABI=0 ")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_GLIBCXX_USE_CXX11_ABI=0 -static-libstdc++ ")
  endif()
endif()

#--------------------------------------------------------------
# 7. Linker flags
#--------------------------------------------------------------
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
set(LibDSphMoorDyn "dsphmoordyn_64")
  set(LIBNAME "${LibDSphMoorDyn}_${CMAKE_VS_PLATFORM_NAME}")
  add_library(${LibDSphMoorDyn} STATIC ${LibDSphMoorDyn_FILES})
elseif(MSVC)
    # MSVC 2015
    set(LibDSphMoorDyn "LibDSphMoorDyn")
    set(LIBNAME "${LibDSphMoorDyn}_${CMAKE_VS_PLATFORM_NAME}_${CMAKE_VS_PLATFORM_TOOLSET}_${CMAKE_CONFIGURATION_TYPES}")
    add_library(${LIBNAME} STATIC ${LibDSphMoorDyn_FILES})
  foreach(CPP_FLAGS CMAKE_CXX_FLAGS CXX_FLAGS_MoorDyn)
    if(${CPP_FLAGS} MATCHES "/MD")
      string(REGEX REPLACE "/MD" "/MT" ${CPP_FLAGS} "${CPP_FLAGS}")
    endif()
  endforeach(CPP_FLAGS)
endif()
