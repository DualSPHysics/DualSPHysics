cmake_minimum_required(VERSION 3.15)

project(LibDSphMoorDynPlus LANGUAGES CXX)

# Define build type defaults solo si no existe 
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE) 
endif() 
# Forzar uso de /MT o /MTd si estás en Visual Studio 
if(MSVC) 
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MT") 
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MTd") 
endif() 

# --------------------------
# 1. Source files
# --------------------------
set(SRC DSphMoorDynPlus.cpp FunMoorDynPlus.cpp IBody.cpp IConnection.cpp IEnvironment.cpp ILine.cpp IQSlines.cpp IOutput.cpp Misc.cpp MoorDynPlus.cpp )
set(HEADER FunMoorDynPlus.h IApi.h IBody.h IConnection.h IEnvironment.h ILine.h IQSlines.h IOutput.h TypesMoorDynPlus.h Misc.h MoorDynPlus.h )
set(SRC_EXTRA kiss_fft.c)
set(HEADER_EXTRA _kiss_fft_guts.h kiss_fft.h kissfft.hh)
set(DSPH_SRC ../Functions.cpp ../FunGeo3d.cpp ../JDataArrays.cpp ../JException.cpp ../JLog2.cpp ../JObject.cpp ../JSaveCsv2.cpp ../JSpVtkShape.cpp ../JRangeFilter.cpp ../JXml.cpp ../tinystr.cpp ../tinyxml.cpp ../tinyxmlerror.cpp ../tinyxmlparser.cpp )
set(DSPH_HEADERS ../DualSphDef.h ../Functions.h ../FunGeo3d.h ../FunctionsMath.h ../JDataArrays.h ../JException.h ../JLog2.h ../JMatrix4.h ../JObject.h ../JSaveCsv2.h ../JRangeFilter.h ../JSpVtkShape.h ../JXml.h ../JXmlDef.h ../RunExceptionDef.h ../tinystr.h ../tinyxml.h )

# --------------------------
# 2. Folder organisation files
# --------------------------
source_group(source FILES ${SRC} ${HEADER} ${SRC_EXTRA} ${HEADER_EXTRA})
source_group(tools FILES ${DSPH_HEADERS} ${DSPH_SRC})
set(LibDSphMoorDynPlus_FILES ${SRC} ${HEADER} ${SRC_EXTRA} ${HEADER_EXTRA} ${DSPH_HEADERS} ${DSPH_SRC} )

# --------------------------
# 3. Include directories
# --------------------------
include_directories( ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/.. )

#--------------------------------------------------------------
# 4. Don't use C++11 ABI for compatibility between GCC v4.9 and >=5
#--------------------------------------------------------------
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  message(STATUS "GCC version: ${CMAKE_CXX_COMPILER_VERSION}")
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -c -O3 -static -static-libgcc -DLINUX ")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x -c -O3 -static -static-libgcc -std=gnu++0x -DLINUX ")
  if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 4.9.9) 
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D_GLIBCXX_USE_CXX11_ABI=0 ")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_GLIBCXX_USE_CXX11_ABI=0 -static-libstdc++ ")
  endif()
endif()

#--------------------------------------------------------------
# 5. Linker flags
#--------------------------------------------------------------
set(LIBNAME LibDSphMoorDynPlus)
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  add_library(${LIBNAME} STATIC ${LibDSphMoorDynPlus_FILES})
  set_target_properties(${LIBNAME} PROPERTIES 
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/../../lib/linux_gcc"
    OUTPUT_NAME "dsphmoordynplus_64")
elseif(MSVC)
  add_library(${LIBNAME} STATIC ${LibDSphMoorDynPlus_FILES}) 
  target_compile_definitions(${LIBNAME} PRIVATE DISABLE_DSPH)
  foreach(CONFIG_TYPE Debug Release RelWithDebInfo MinSizeRel)
    string(TOUPPER "${CONFIG_TYPE}" CONFIG_TYPE_UPPER)
    set_target_properties(${LIBNAME} PROPERTIES 
      ARCHIVE_OUTPUT_DIRECTORY_${CONFIG_TYPE_UPPER} "${CMAKE_CURRENT_SOURCE_DIR}/../../lib/vs2022" 
      OUTPUT_NAME "${LIBNAME}_${CMAKE_VS_PLATFORM_NAME}_${CMAKE_VS_PLATFORM_TOOLSET}_$<CONFIG>"
    )
  endforeach() 
  foreach(CPP_FLAGS CMAKE_CXX_FLAGS CXX_FLAGS_MoorDynPlus) 
    if(${CPP_FLAGS} MATCHES "/MD")
      string(REGEX REPLACE "/MD" "/MT" ${CPP_FLAGS} "${CPP_FLAGS}")
    endif()
  endforeach(CPP_FLAGS)
endif()